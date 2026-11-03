#!/usr/bin/env python3
"""Negative controls for the fail-closed checkers (standard library only).

For each checker this script
  1. confirms statically (via ast) that it contains no ``assert`` statement;
  2. runs the unmodified checker under ``python -I -B`` and ``python -I -B -O``
     and requires exit status 0 with ``all_checks_passed: true``;
  3. writes deliberately mutated copies (one textual mutation each, which must
     match the source exactly once) into a fresh temporary work directory,
     runs each under both modes, and requires exit status 1, a JSON report with
     ``all_checks_passed: false`` and at least one failed check label.

The temporary work directory (under the system temp location) is never
deleted by this script. Results are written to
results/negative_controls.json next to this script (or --output).
"""

from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent            # checks/
WORK = Path(tempfile.mkdtemp(prefix="ent_concav_negative_controls_"))

SCRIPTS = {
    "coefficient_checks.py": "coefficient_checks.json",
    "local_positive_checks.py": "local_positive_checks.json",
    "identity_checks.py": "identity_checks.json",
}

MUTANTS = [
    # coefficient_checks.py
    ("coef_target_minus7_gamma", "coefficient_checks.py", "wrong curvature target -7 gamma instead of -8 gamma",
     "target_curvature_per_n=-8 * gamma", "target_curvature_per_n=-7 * gamma"),
    ("coef_b0_identity_minus7n", "coefficient_checks.py", "exact b=0 identity tested as C''=-7nC",
     "c2 == -8 * n * c0", "c2 == -7 * n * c0"),
    ("coef_gamma_closed_form", "coefficient_checks.py", "gamma closed form 20 sqrt5 q/(49q-4)",
     "gamma_closed = 21 * math.sqrt(5)", "gamma_closed = 20 * math.sqrt(5)"),
    ("coef_norm_prefactor", "coefficient_checks.py", "norm prefactor (M23) 4.4/(pi Delta sqrt det H0)",
     "prefactor = 4 / (math.pi * delta", "prefactor = 4.4 / (math.pi * delta"),
    ("coef_matrix_second_jet", "coefficient_checks.py", "B''(0) entry -3e instead of -4e",
     "[z, -4 * e, z, -4 * c],", "[z, -3 * e, z, -4 * c],"),
    # local_positive_checks.py
    ("lp_928_to_929", "local_positive_checks.py", "half-signal constant 928 -> 929 (auxiliary group)",
     "sum_b_squared == 928", "sum_b_squared == 929"),
    ("lp_133_to_134", "local_positive_checks.py", "half-signal constant 133 -> 134 (auxiliary group)",
     "sum_n_kernel == 133", "sum_n_kernel == 134"),
    ("lp_table_s_digit", "local_positive_checks.py", "Section I.7 table: s(0.50) last digit 5 -> 6",
     's="0.5202847075"', 's="0.5202847076"'),
    ("lp_table_radius_digit", "local_positive_checks.py", "Section I.7 table: radius(0.50) 1.7134e-5 -> 1.7135e-5",
     'radius="1.7134e-5"', 'radius="1.7135e-5"'),
    ("lp_J_exponent", "local_positive_checks.py", "(P10) J with (pi/eps)^0.26",
     'J = c(s) / direct_c_b * direct_beta_b / beta(s) * (mp.pi / eps) ** mp.mpf("0.25")',
     'J = c(s) / direct_c_b * direct_beta_b / beta(s) * (mp.pi / eps) ** mp.mpf("0.26")'),
    ("lp_A_constant", "local_positive_checks.py", "(P21) A_q with 30 instead of 60",
     "A = K * J / 8 + 60 + J / (8 * eta**2)", "A = K * J / 8 + 30 + J / (8 * eta**2)"),
    # identity_checks.py
    ("id_wrong_Delta", "identity_checks.py", "(M20) Delta with (1+a+a/1000) prefactor",
     "delta = (1 + A) * (", "delta = (1 + A + A / 1000) * ("),
    ("id_wrong_M21_entry", "identity_checks.py", "(M21) entry B_33 = d+3eu instead of d+2eu",
     "c * w, d + 2 * e * u, e * w]", "c * w, d + 3 * e * u, e * w]"),
    ("id_wrong_lemma4_block", "identity_checks.py", "Lemma II.5 off-diagonal block beta v_t v_t^T",
     "m[i][2 + j] = bt * vt[i] * vs[j]", "m[i][2 + j] = bt * vt[i] * vt[j]"),
    ("id_wrong_4nfact2", "identity_checks.py", "(M20) prefactor 3(n!)^2 instead of 4(n!)^2",
     '4 * FACT[n] ** 2 / p["Delta"] * coefficient_cn(bm, n)', '3 * FACT[n] ** 2 / p["Delta"] * coefficient_cn(bm, n)'),
    ("id_wrong_generating_function", "identity_checks.py", "double generating function w^2 coefficient perturbed",
     "d_ = -A / (2 * (1 + A))", 'd_ = -A / (2 * (1 + A)) * D("1.000001")'),
]

MODES = {"normal": ["-I", "-B"], "optimized_O": ["-I", "-B", "-O"]}


def asserts_in(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return [node.lineno for node in ast.walk(tree) if isinstance(node, ast.Assert)]


def run_script(script_path, out_dir, mode, json_name):
    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / json_name
    if report.exists():
        report.write_text("", encoding="utf-8")   # never read a stale report
    start = time.perf_counter()
    proc = subprocess.run([sys.executable, *MODES[mode], str(script_path), "--output-dir", str(out_dir)],
                          cwd=str(out_dir), capture_output=True, text=True, timeout=600)
    elapsed = time.perf_counter() - start
    try:
        data = json.loads(report.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        data = None
    return dict(mode=mode, exit_code=proc.returncode, seconds=round(elapsed, 3),
                report_found=data is not None,
                all_checks_passed=None if data is None else data.get("all_checks_passed"),
                python_optimize_flag=None if data is None else data.get("python_optimize_flag"),
                number_of_checks=None if data is None else data.get("number_of_checks"),
                failed_labels=[] if data is None else data.get("failed_labels", [])[:8],
                stderr_tail=proc.stderr[-400:])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "results" / "negative_controls.json")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    start = time.perf_counter()
    problems = []
    result = {"python": sys.version, "work_directory": "system temporary directory (not recorded)", "static_no_assert": {},
              "originals": [], "mutants": []}

    for name in SCRIPTS:
        lines = asserts_in(HERE / name)
        result["static_no_assert"][name] = lines
        if lines:
            problems.append(f"{name}: assert statements at lines {lines}")

    tasks = []
    for name, json_name in SCRIPTS.items():
        for mode in MODES:
            tasks.append(("original", name, None, HERE / name, WORK / "original" / f"{name[:-3]}_{mode}", mode, json_name))
    for mid, name, desc, old, new in MUTANTS:
        src = (HERE / name).read_text(encoding="utf-8")
        count = src.count(old)
        entry = dict(id=mid, script=name, description=desc, mutation_old=old, mutation_new=new,
                     occurrences=count, runs=[])
        result["mutants"].append(entry)
        if count != 1:
            problems.append(f"mutant {mid}: mutation target occurs {count} times (expected exactly 1)")
            continue
        mdir = WORK / mid / "checks"
        mdir.mkdir(parents=True, exist_ok=True)
        mpath = mdir / name
        mpath.write_text(src.replace(old, new), encoding="utf-8")
        for mode in MODES:
            tasks.append(("mutant", name, entry, mpath, WORK / mid / f"out_{mode}", mode, SCRIPTS[name]))

    def work(task):
        kind, name, entry, path, out_dir, mode, json_name = task
        return task, run_script(path, out_dir, mode, json_name)

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        outcomes = list(pool.map(work, tasks))

    for (kind, name, entry, path, out_dir, mode, json_name), run in outcomes:
        expected_flag = 1 if mode == "optimized_O" else 0
        if run["python_optimize_flag"] != expected_flag and run["report_found"]:
            problems.append(f"{kind} {name} {mode}: optimize flag {run['python_optimize_flag']} != {expected_flag}")
        if kind == "original":
            ok = run["exit_code"] == 0 and run["all_checks_passed"] is True
            run.update(script=name, expected="exit 0, all_checks_passed true", ok=ok)
            result["originals"].append(run)
            if not ok:
                problems.append(f"original {name} {mode}: exit {run['exit_code']}, passed={run['all_checks_passed']}")
        else:
            ok = (run["exit_code"] == 1 and run["all_checks_passed"] is False and bool(run["failed_labels"]))
            run.update(expected="exit 1, all_checks_passed false, >=1 failed label", ok=ok)
            entry["runs"].append(run)
            if not ok:
                problems.append(f"mutant {entry['id']} {mode}: exit {run['exit_code']}, passed={run['all_checks_passed']}")

    for entry in result["mutants"]:
        entry["caught_in_both_modes"] = len(entry["runs"]) == 2 and all(r["ok"] for r in entry["runs"])
    result["number_of_mutants"] = len(MUTANTS)
    result["mutants_caught_in_both_modes"] = sum(e["caught_in_both_modes"] for e in result["mutants"])
    result["originals_pass_in_both_modes"] = (len(result["originals"]) == 2 * len(SCRIPTS)
                                              and all(r["ok"] for r in result["originals"]))
    result["problems"] = problems
    result["all_controls_passed"] = not problems and result["originals_pass_in_both_modes"] \
        and result["mutants_caught_in_both_modes"] == len(MUTANTS)
    result["elapsed_seconds"] = round(time.perf_counter() - start, 3)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for r in result["originals"]:
        print(f"original {r['script']:26s} {r['mode']:12s} exit={r['exit_code']} checks={r['number_of_checks']} "
              f"{r['seconds']:.2f}s {'OK' if r['ok'] else 'PROBLEM'}")
    for e in result["mutants"]:
        modes = " ".join(f"{r['mode']}:exit={r['exit_code']}" for r in e["runs"])
        print(f"mutant {e['id']:30s} {modes:34s} {'CAUGHT' if e['caught_in_both_modes'] else 'NOT CAUGHT'}")
    for p in problems:
        print("PROBLEM:", p)
    print(f"negative controls: {'ALL PASSED' if result['all_controls_passed'] else 'FAILED'}; "
          f"{result['mutants_caught_in_both_modes']}/{len(MUTANTS)} mutants caught in both modes; "
          f"{result['elapsed_seconds']:.1f} s")
    sys.exit(0 if result["all_controls_passed"] else 1)


if __name__ == "__main__":
    main()

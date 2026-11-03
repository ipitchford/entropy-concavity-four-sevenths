#!/usr/bin/env python3
"""Fail-closed coefficient and curvature checks for Part II of the entropy manuscript.

Standard library only. Every check goes through ``Checker.check``; the pass flag
is computed from the recorded results and any failure (or uncaught exception)
ends the run with exit status 1. There is no ``assert`` statement, so running
under ``python -O`` changes nothing.

What is checked (manuscript Part II):

* closed forms (M4) b, (M13)-(M15) d, e, c, kappa, -alpha, d(4/7), (M16) the
  exponential comparison and its factorisation, (M22)/(M3) gamma(q)
  = (kappa+b)/(kappa-b) = 21 sqrt5 q/(49q-4), (M24) det H_0, and the
  eigenvector/jet identities used in the proof of Lemma II.4;
* thirty floating-point cases (q in {4/7, 0.575, 0.6, 0.75, 1},
  n in {2,4,8,12,16,24}) of the Gaussian-moment recurrence: C_n'(0)=0,
  |C''/(nC) - (-8 gamma)| <= C_CURV/n, observed O(1/n) decay, and
  |normalised norm / (M23) prefactor - 1| <= C_NORM/n;
* exact rational surrogate: recurrence equals an independent multigraph
  enumeration (n=2,4,6), and the exact b=0 identity C''/C = -8n (Section II.7).

These finite computations are sanity checks; they are not a proof of the
asymptotic statements.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
import time
import traceback
from fractions import Fraction as F
from pathlib import Path

# Tolerance constants for the O(1/n) asymptotic comparisons (see CHANGES.md).
# Observed max over q in {4/7,0.575,0.6,0.75,1} and n in {2,...,24}:
#   n * |C''/(nC) + 8 gamma|          <= 1.77  (q = 4/7, n = 24)
#   n * |norm/(M23) prefactor - 1|    <= 0.12  (q = 1,   n = 24)
# The constants below give a margin of about 1.13x and 2x respectively, so a
# wrong target (e.g. -7 gamma, which is off by gamma >= 1.04 at every n) or a
# wrong prefactor fails, while the true O(1/n) error passes.
C_CURV = 2.0
C_NORM = 0.25
# The ratio err(24)/err(12) for an exact 1/n law is 0.5.
DECAY_WINDOW = (0.45, 0.55)


class Checker:
    def __init__(self):
        self.results = []

    def check(self, cond, label, detail="", group="manuscript", ref=None):
        ok = bool(cond)
        self.results.append(dict(label=label, passed=ok, detail=str(detail),
                                 group=group, manuscript_ref=ref))
        if not ok:
            print(f"FAIL [{group}] {label}: {detail}", flush=True)
        return ok

    def summary(self):
        failed = [r["label"] for r in self.results if not r["passed"]]
        manuscript = [r for r in self.results if r["group"] == "manuscript"]
        return dict(
            number_of_checks=len(self.results),
            number_failed=len(failed),
            failed_labels=failed,
            manuscript_claim_checks=len(manuscript),
            manuscript_claim_checks_passed=all(r["passed"] for r in manuscript) and bool(manuscript),
            all_checks_passed=(not failed) and bool(self.results),
        )


def rel_close(x, y, tol):
    return abs(x - y) <= tol * max(abs(x), abs(y))


def matrix_jets(d, e, c):
    """B(0), B'(0), B''(0) for the matrix (M21); w = sqrt(1-4u^2), w'(0)=0, w''(0)=-4."""
    z = d * 0
    b0 = [[d, e, c, c], [e, d, c, c], [c, c, d, e], [c, c, e, d]]
    b1 = [
        [2 * e, z, 2 * c, z],
        [z, -2 * e, z, -2 * c],
        [2 * c, z, 2 * e, z],
        [z, -2 * c, z, -2 * e],
    ]
    b2 = [
        [z, -4 * e, z, -4 * c],
        [-4 * e, z, -4 * c, z],
        [z, -4 * c, z, -4 * e],
        [-4 * c, z, -4 * e, z],
    ]
    return b0, b1, b2


def diagonal_jets(d, e, c, ns):
    """Coefficient C_n(0) and its first two u-derivatives via the recurrence
    alpha_i c_alpha = sum_j B_ij c_(alpha-e_i-e_j), differentiated twice."""
    maximum = max(ns)
    base = maximum + 1
    strides = (base**3, base**2, base, 1)
    length = base**4
    zero = d * 0
    one = zero + 1
    values = [[zero] * length for _ in range(3)]
    values[0][0] = one
    b0, b1, b2 = matrix_jets(d, e, c)
    c0, c1, c2 = values
    wanted = set(ns)
    result = {}
    for a0 in range(base):
        for a1 in range(base):
            for a2 in range(base):
                prefix = a0 * strides[0] + a1 * strides[1] + a2 * strides[2]
                for a3 in range((a0 + a1 + a2) % 2, base, 2):
                    alpha = (a0, a1, a2, a3)
                    if alpha == (0, 0, 0, 0):
                        continue
                    i = next(j for j, count in enumerate(alpha) if count)
                    index = prefix + a3
                    v0 = v1 = v2 = zero
                    for j in range(i, 4):
                        if alpha[j] <= (j == i):
                            continue
                        previous = index - strides[i] - strides[j]
                        p0, p1, p2 = c0[previous], c1[previous], c2[previous]
                        a, b, h = b0[i][j], b1[i][j], b2[i][j]
                        v0 += a * p0
                        v1 += a * p1 + b * p0
                        v2 += a * p2 + 2 * b * p1 + h * p0
                    c0[index] = v0 / alpha[i]
                    c1[index] = v1 / alpha[i]
                    c2[index] = v2 / alpha[i]
                    if a0 == a1 == a2 == a3 and a0 in wanted:
                        result[a0] = (c0[index], c1[index], c2[index])
    return result


def enumerate_edge_count_jets(d: F, e: F, c: F, n: int):
    """Independent exact expansion over six nonloop and four loop counts."""
    b0, b1, b2 = matrix_jets(d, e, c)
    off_diagonal = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    total = [F(0), F(0), F(0)]
    configurations = 0

    def visit(which, degrees, weight, score, log_curvature):
        nonlocal configurations
        if which == len(off_diagonal):
            if any(degree % 2 for degree in degrees):
                return
            for i, remaining in enumerate(degrees):
                k = remaining // 2
                a, b, h = b0[i][i] / 2, b1[i][i] / 2, b2[i][i] / 2
                weight *= a**k / math.factorial(k)
                score += k * b / a
                log_curvature += k * (h / a - (b / a) ** 2)
            total[0] += weight
            total[1] += weight * score
            total[2] += weight * (score**2 + log_curvature)
            configurations += 1
            return
        i, j = off_diagonal[which]
        a, b, h = b0[i][j], b1[i][j], b2[i][j]
        for k in range(min(degrees[i], degrees[j]) + 1):
            remainder = list(degrees)
            remainder[i] -= k
            remainder[j] -= k
            visit(which + 1, tuple(remainder), weight * a**k / math.factorial(k),
                  score + k * b / a, log_curvature + k * (h / a - (b / a) ** 2))

    visit(0, (n, n, n, n), F(1), F(0), F(0))
    return tuple(total), configurations


def det4(m):
    a = [row[:] for row in m]
    det = 1.0
    for col in range(4):
        piv = max(range(col, 4), key=lambda r: abs(a[r][col]))
        if a[piv][col] == 0:
            return 0.0
        if piv != col:
            a[col], a[piv] = a[piv], a[col]
            det = -det
        det *= a[col][col]
        for r in range(col + 1, 4):
            f = a[r][col] / a[col][col]
            for k in range(col, 4):
                a[r][k] -= f * a[col][k]
    return det


def matvec(m, v):
    return [sum(m[i][j] * v[j] for j in range(4)) for i in range(4)]


def ratio_string(value: F):
    return f"{value.numerator}/{value.denominator}"


def parameters(q, ck, label):
    ref = "Part II (M4),(M13)-(M16),(M22),(M24),(M3)"
    a = 7 / (3 * math.sqrt(5))
    b = (a - 1) / (a + 1)
    lp = (a * (1 + q) - 1) / (a * (1 + q) + 1)
    lm = (a * (1 - q) - 1) / (a * (1 - q) + 1)
    alpha = (lp + lm) / 2 - b
    beta = (lp - lm) / 2
    d, e, c = -b - alpha / 2, -alpha / 2, beta / 2
    kappa = -lm
    delta = (1 + a) * math.sqrt((1 + a) ** 2 - a * a * q * q)
    hdet = 2 * (1 - lp / kappa) * (1 - b / kappa) ** 2
    gamma = (kappa + b) / (kappa - b)
    gamma_closed = 21 * math.sqrt(5) * q / (49 * q - 4)
    prefactor = 4 / (math.pi * delta * math.sqrt(hdet))
    b0, b1, b2 = matrix_jets(d, e, c)

    ck.check(min(d, e, c) > 0, f"q={label}: d,e,c>0 (M14)", (d, e, c), ref=ref)
    ck.check(rel_close(d + e + 2 * c, kappa, 1e-13), f"q={label}: d+e+2c=kappa (M15)",
             (d + e + 2 * c, kappa), ref=ref)
    ck.check(abs((e - d) - b) <= 1e-15, f"q={label}: e-d=b (M15)", (e - d, b), ref=ref)
    ck.check(rel_close(-alpha, 2 * a * a * q * q / ((1 + a) * ((1 + a) ** 2 - a * a * q * q)), 1e-12),
             f"q={label}: closed form of -alpha (M14)", -alpha, ref=ref)
    ck.check(rel_close(gamma, gamma_closed, 1e-12),
             f"q={label}: gamma=(kappa+b)/(kappa-b)=21sqrt5 q/(49q-4) (M22),(M3)",
             (gamma, gamma_closed), ref=ref)
    h0 = [[(1.0 if i == j else 0.0) + b0[i][j] / kappa for j in range(4)] for i in range(4)]
    ck.check(rel_close(det4(h0), hdet, 1e-12), f"q={label}: det(I+B(0)/kappa) equals (M24)",
             (det4(h0), hdet), ref=ref)
    one = [1.0] * 4
    v = [1.0, -1.0, 1.0, -1.0]
    w1 = [1.0, 1.0, -1.0, -1.0]
    w2 = [1.0, -1.0, -1.0, 1.0]
    U = kappa + b
    eig_ok = (all(abs(x - kappa) < 1e-14 for x in matvec(b0, one))
              and all(abs(x + b * y) < 1e-14 for x, y in zip(matvec(b0, v), v))
              and all(abs(x + lp * y) < 1e-14 for x, y in zip(matvec(b0, w1), w1))
              and all(abs(x + b * y) < 1e-14 for x, y in zip(matvec(b0, w2), w2)))
    ck.check(eig_ok, f"q={label}: B(0) eigenpairs kappa,-lambda_+,-b,-b (Lemma II.4)", "", ref=ref)
    jet_ok = (all(abs(x - U * y) < 1e-14 for x, y in zip(matvec(b1, one), v))
              and abs(sum(matvec(b2, one)) + 8 * U) < 1e-13)
    ck.check(jet_ok, f"q={label}: B'(0)1=Uv and 1^T B''(0) 1=-8U (Lemma II.4)", "", ref=ref)
    tau_half = (1 - a + a * q / 2) / (1 + a - a * q / 2)
    factor = -(a / 45) * (7 * q - 4) * (7 * q * q - 24 * q + 2) / ((1 + a - a * q) ** 2 * (1 + a - a * q / 2))
    ck.check(abs((kappa**2 - tau_half) - factor) <= 1e-14,
             f"q={label}: kappa^2 - tau(q/2) factorisation (M16)", (kappa**2 - tau_half, factor), ref=ref)
    if label == "4/7":
        ck.check(abs(kappa**2 - tau_half) < 1e-14, "q=4/7: equality kappa^2 = tau(2/7) (M16)",
                 kappa**2 - tau_half, ref=ref)
    else:
        ck.check(kappa**2 - tau_half > 0, f"q={label}: strict tau(q/2) < kappa^2 (M16)",
                 kappa**2 - tau_half, ref=ref)
    ck.check(math.sqrt(b) < kappa, f"q={label}: sqrt(b) < kappa (M16 second margin)", (math.sqrt(b), kappa), ref=ref)
    return dict(a=a, b=b, d=d, e=e, c=c, kappa=kappa, delta=delta,
                lambda_plus=lp, hessian_determinant=hdet, gamma=gamma,
                gamma_closed_form=gamma_closed,
                target_curvature_per_n=-8 * gamma, target_norm_prefactor=prefactor)


def run(ck, results, record):
    a = 7 / (3 * math.sqrt(5))
    b = (a - 1) / (a + 1)
    ck.check(abs(b - (47 - 21 * math.sqrt(5)) / 2) < 1e-15 and 0 < b < 1 / 40,
             "b=(a-1)/(a+1)=(47-21sqrt5)/2 in (0,1/40) (M4)", b, ref="Part II (M4)")
    q47 = 4 / 7
    lp = (a * (1 + q47) - 1) / (a * (1 + q47) + 1)
    lm = (a * (1 - q47) - 1) / (a * (1 - q47) + 1)
    d47 = -b - ((lp + lm) / 2 - b) / 2
    ck.check(abs(d47 - (203 * math.sqrt(5) - 453) / 38) < 1e-14 and 5 * 203**2 - 453**2 == 836,
             "d(4/7)=(203sqrt5-453)/38>0 (M14)", d47, ref="Part II (M14)")

    ns = [2, 4, 8, 12, 16, 24]
    cases = [("4/7", 4 / 7), ("0.575", .575), ("0.6", .6), ("0.75", .75), ("1", 1.0)]
    numeric_rows = 0
    for label, q in cases:
        scope = ("coefficient_only: Lemma II.4 applies to the coefficient at q=4/7 (entries positive by (M14)), "
                 "but Theorem B and (M16) do not; no entropy conclusion") if label == "4/7" else "q>4/7"
        p = parameters(q, ck, label)
        jets = diagonal_jets(p["d"], p["e"], p["c"], ns)
        case = {"q_label": label, "q": q, "scope": scope, "parameters": p, "rows": []}
        record(f"q={label}; target C''/(n C)={p['target_curvature_per_n']:.12g}; "
               f"target norm prefactor={p['target_norm_prefactor']:.12g}; scope={'coefficient-only' if label == '4/7' else 'q>4/7'}")
        errs = {}
        for n in ns:
            c0, c1, c2 = jets[n]
            ratio = c2 / (n * c0)
            norm_squared = 4 * math.factorial(n) ** 2 * c0 / p["delta"]
            norm_prefactor = n * norm_squared / p["kappa"] ** (2 * n)
            first_relative = c1 / c0
            curv_err = ratio - p["target_curvature_per_n"]
            norm_err = norm_prefactor / p["target_norm_prefactor"] - 1
            errs[n] = curv_err
            ref22 = "Part II (M22)"
            ck.check(c0 > 0 and math.isfinite(c0), f"q={label} n={n}: C_n(0)>0 finite", c0, ref=ref22)
            ck.check(abs(first_relative) < 1e-10 * max(1.0, n), f"q={label} n={n}: C_n'(0)=0 (evenness)",
                     first_relative, ref=ref22)
            ck.check(ratio < 0, f"q={label} n={n}: C_n''(0)<0", ratio, ref=ref22)
            ck.check(abs(curv_err) <= C_CURV / n,
                     f"q={label} n={n}: |C''/(nC) + 8 gamma| <= {C_CURV}/n",
                     f"error={curv_err:.6g}, bound={C_CURV / n:.6g}", ref=ref22)
            ck.check(abs(norm_err) <= C_NORM / n,
                     f"q={label} n={n}: |n||M_1/2||^2 kappa^(-2n) / (M23) prefactor - 1| <= {C_NORM}/n",
                     f"error={norm_err:.6g}, bound={C_NORM / n:.6g}", ref="Part II (M23),(M24)")
            row = dict(n=n, coefficient=c0, coefficient_first_derivative=c1,
                       coefficient_second_derivative=c2, relative_first_derivative=first_relative,
                       curvature_per_n=ratio, curvature_error=curv_err,
                       curvature_error_times_n=curv_err * n,
                       norm_squared=norm_squared, normalised_norm_prefactor=norm_prefactor,
                       norm_prefactor_relative_error=norm_err,
                       norm_prefactor_relative_error_times_n=norm_err * n)
            case["rows"].append(row)
            numeric_rows += 1
            record(f"  n={n:2d}  C''/(n C)={ratio: .12g}  n*err={curv_err * n:+.4f}  "
                   f"normalised norm={norm_prefactor:.12g}  n*rel.norm err={norm_err * n:+.4f}")
        decay = errs[24] / errs[12]
        ck.check(DECAY_WINDOW[0] <= decay <= DECAY_WINDOW[1],
                 f"q={label}: curvature error decays like 1/n (err24/err12 in {DECAY_WINDOW})",
                 decay, ref="Part II (M22): O(1) remainder")
        ck.check(all(errs[m] > errs[n] > 0 for m, n in zip(ns, ns[1:])),
                 f"q={label}: curvature error positive and decreasing in n", "", ref="Part II (M22)")
        case["curvature_error_ratio_n24_over_n12"] = decay
        results["numerical"].append(case)
    ck.check(numeric_rows == 30, "thirty numerical coefficient cases evaluated", numeric_rows,
             ref="Section II.7")
    q06 = next(cse for cse in results["numerical"] if cse["q_label"] == "0.6")
    r24 = next(r for r in q06["rows"] if r["n"] == 24)["curvature_per_n"]
    ck.check(f"{r24:.12g}" == "-8.8098805868" and abs(r24 + 8.80988058680) < 5e-12,
             "Section II.7 quoted value C''/(nC)=-8.80988058680 at q=0.6, n=24", r24, ref="Section II.7")
    ck.check(abs(q06["parameters"]["target_curvature_per_n"] + 8.87384457213) < 5e-12,
             "Section II.7 quoted limit -8.87384457213 at q=0.6",
             q06["parameters"]["target_curvature_per_n"], ref="Section II.7")

    d, e, c = F(1, 20), F(1, 10), F(7, 40)
    kappa, bb = d + e + 2 * c, e - d
    target = -8 * (kappa + bb) / (kappa - bb)
    ck.check(target == F(-88, 9), "exact surrogate target -8(kappa+b)/(kappa-b) = -88/9", target,
             ref="Part II (M22)")
    jets = diagonal_jets(d, e, c, [2, 4, 6])
    record("Exact rational surrogate d=1/20, e=1/10, c=7/40; target curvature per n=-88/9")
    previous_error = None
    for n in [2, 4, 6]:
        direct, configurations = enumerate_edge_count_jets(d, e, c, n)
        ck.check(jets[n] == direct, f"surrogate n={n}: recurrence equals multigraph enumeration",
                 configurations, ref="Part II (M20)-(M21) coefficient C_n")
        c0, c1, c2 = direct
        ck.check(c1 == 0 and c2 < 0, f"surrogate n={n}: C'=0 and C''<0 exactly", "", ref="Part II (M22)")
        ratio = c2 / (n * c0)
        err = ratio - target
        ck.check(err > 0 and (previous_error is None or err < previous_error),
                 f"surrogate n={n}: exact curvature error positive and decreasing toward -88/9",
                 float(err), ref="Part II (M22)")
        previous_error = err
        results["exact_surrogate"].append(dict(
            n=n, coefficient=ratio_string(c0), first_derivative=ratio_string(c1),
            second_derivative=ratio_string(c2), curvature_per_n=ratio_string(ratio),
            curvature_per_n_decimal=float(ratio),
            independent_enumeration_configurations=configurations,
            recurrence_equals_enumeration=jets[n] == direct))
        record(f"  n={n}: recurrence vs {configurations} edge-count configurations; "
               f"C''/(n C)={ratio_string(ratio)} = {float(ratio):.12g}")

    jets = diagonal_jets(F(1, 10), F(1, 10), F(1, 5), [2, 4, 6])
    record("Exact b=0 identity check: C''/C = -8n")
    for n in [2, 4, 6]:
        c0, c1, c2 = jets[n]
        ok = ck.check(c1 == 0 and c2 == -8 * n * c0, f"b=0 n={n}: C''/C = -8n exactly",
                      ratio_string(c2 / c0), ref="Section II.7 (b=0 identity)")
        results["exact_b_zero"].append(dict(n=n, exact_curvature=ratio_string(c2 / c0), passed=ok))
        record(f"  n={n}: C''/C={ratio_string(c2 / c0)}; {'PASS' if ok else 'FAIL'}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    log = []

    def record(line):
        print(line, flush=True)
        log.append(line)

    ck = Checker()
    results = {
        "purpose": "Finite exact and floating-point checks; not a proof of asymptotics.",
        "python_version": platform.python_version(),
        "python_optimize_flag": sys.flags.optimize,
        "recurrence": "alpha_i*c_alpha = sum_j B_ij*c_(alpha-e_i-e_j), differentiated twice",
        "tolerances": dict(C_CURV=C_CURV, C_NORM=C_NORM, decay_window=DECAY_WINDOW),
        "numerical": [], "exact_surrogate": [], "exact_b_zero": [],
    }
    try:
        run(ck, results, record)
    except Exception as exc:  # fail closed on any crash
        ck.check(False, "uncaught exception", repr(exc) + "\n" + traceback.format_exc())
    results.update(ck.summary())
    results["checks"] = ck.results
    results["elapsed_seconds"] = time.perf_counter() - start
    status = "All checks passed" if results["all_checks_passed"] else f"{results['number_failed']} CHECK(S) FAILED"
    record(f"{status}: {results['number_of_checks']} checks in {results['elapsed_seconds']:.3f} seconds.")
    (args.output_dir / "coefficient_checks.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    (args.output_dir / "coefficient_checks.log").write_text("\n".join(log) + "\n", encoding="utf-8")
    sys.exit(0 if results["all_checks_passed"] else 1)


if __name__ == "__main__":
    main()

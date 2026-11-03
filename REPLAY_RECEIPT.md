# Replay receipt (release 0.1.0-candidate, 3 November 2026)

Clean replay before freezing, from a fresh copy of the approved package (not the working directory), CPython 3.12.11
(macOS arm64, standard library only), `PYTHONDONTWRITEBYTECODE=1`, run inside `checks/`:

| Command | Normal (`-I -B`) | Optimized (`-I -B -O`) |
|---|---|---|
| `coefficient_checks.py --output-dir ...` | exit 0, 233 checks, 2 s | exit 0, 2 s |
| `local_positive_checks.py --output-dir ...` | exit 0, 165 checks, <1 s | exit 0, <1 s |
| `identity_checks.py --output-dir ...` | exit 0, 178 checks, 3 s | exit 0, 4 s |
| `extension_checks.py --output-dir ...` | exit 0, 315 checks, 2 s | exit 0, 2 s |
| `negative_controls.py --output ...` (`-I -B`; runs both modes itself) | exit 0, 16/16 mutants caught in both modes, 21 s | |

All nine JSON outputs were identical to the shipped files in `checks/results/` apart from runtime and interpreter
fields (`elapsed_seconds`, `runtime_s`, `python_version`, `python`, `work_directory`); the text logs were identical apart
from their final timing line. The only later change to the package was text (DOI insertion in README.md and the
manuscript status section, rebuilt PDF and generated LaTeX, and the release metadata files); the checker code and
results are byte-identical to the replayed copy. CI on GitHub reruns every checker on Linux.

This receipt records execution only. It is not a proof of the analytic statements, which rest on the written proofs.

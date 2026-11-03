# Verification record

Version 0.1.0-candidate, 3 November 2026. This file says what was checked, how, and what each check does not establish.

## Analytic proofs

All theorems rest on the written proofs in [manuscript.md](manuscript.md). No proof assistant was used.

| Result | Location | Internal checking |
|---|---|---|
| Theorem A, Corollary (B) | Part I (Sections I.1-I.7), Section 2.3 | Line-by-line review by a methodology reviewer agent (all closed forms recomputed, Section I.7 table reproduced); different-model adversarial review (no gap found; endpoints, uniformity and logarithmic remainder checked) |
| Theorem B | Part II (Sections II.1-II.7) | Same two reviews; gamma(q), (M20)-(M21), Lemma II.5 matrix independently recomputed to 50+ digits |
| Proposition S, Corollary S.5 | Section III.1 | Proof audit by a reviewer agent (sympy/mpmath recomputation of (S4), (S8), (S9), (S10)); different-model audit. A gap in Corollary S.5(4) (lower bound on the linear rate for all a) was found and closed with an integral representation |
| Proposition E | Section III.2 | Same two audits; constants (E5), (E8), 136/15 recomputed exactly; size bookkeeping checked; scope wording corrected (negativity only near t = 1/2) |
| Remark R | Section III.3 | Both audits; C_delta = 5.92873e-3 and exponent 9/4 recomputed independently |

## Shipped computations (`checks/`, standard library only)

| Script | What it checks | Checks |
|---|---|---|
| coefficient_checks.py | Part II coefficient asymptotics: 30 cases q in {4/7, 0.575, 0.6, 0.75, 1}, n in {2..24}; C''/(nC) against -8 gamma(q) within 2/n and the (M23) prefactor (q > 4/7); exact recurrence vs multigraph enumeration; exact b = 0 identity; closed forms (M14)-(M16), (M22)-(M24) | 233 |
| local_positive_checks.py | Part I recipe (P7), (P9), (P10), (P21) in 100-digit arithmetic; Section I.7 table to displayed digits; kernel square root vs Mehler parametrisation. Constants of a superseded half-signal note are labelled auxiliary | 165 |
| identity_checks.py | (M11), (M20)-(M21) against a direct Hermite expansion (27 cases) and the Lemma II.5 two-weight inner products (36 cases), 80-digit Decimal, agreement about 1e-78 | 178 |
| extension_checks.py | Proposition S envelope s*(r) = b(r) and failure above it; (S9) coincidences at r = 2/7; exact (S8); Proposition E finite-n diagnostics (labelled diagnostics); Remark R exponent and constant; 6 internal negative controls | 315 |
| negative_controls.py | No `assert` in the checkers; originals pass under normal and `-O`; 16 deliberately broken copies fail in both modes | 16 mutants |

Results of the release replay (CPython 3.12, `python3 -I -B` and `-I -B -O`): all checks passed in both modes
([checks/results/](checks/results/)); 16/16 mutants caught.

**Limits.** Floating-point and finite-precision evaluations are not interval certificates. Finite-n coefficient values
support, but do not prove, the asymptotic statements. The q = 4/7 coefficient rows concern the generating coefficient
only, not entropy concavity there. No numerical integration of entropy curves is claimed.

## Replay history

- The first internal version's two checkers were replayed from a fresh copy on 3 November 2026; their outputs
  matched the shipped JSON, but they gated on `assert` and passed under `-O` even when broken. They were replaced by
  the fail-closed checkers above.

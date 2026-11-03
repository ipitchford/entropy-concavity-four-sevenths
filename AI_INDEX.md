# AI index: the critical signal fraction for local entropy concavity (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) and the
paper [entropy_concavity_threshold.pdf](entropy_concavity_threshold.pdf) (source [manuscript.md](manuscript.md)).

## Normalisation

X, Y i.i.d. with even density f, mean 0, variance 1; G standard Gaussian; signal fraction q in (0,1];
F_{q,f}(t) = h(sqrt(qt) X + sqrt(q(1-t)) Y + sqrt(1-q) G), natural logarithms. f = phi (1+g) with g admissible:
real, even, bounded, int g phi = int x^2 g phi = 0. T_rho g(x) = E g(rho x + sqrt(1-rho^2) G); N e_j = j e_j.
C = { q : some delta > 0 makes F_{q, phi(1+g)} concave on [0,1] for all admissible ||g||_inf <= delta }.
SNR = q/(1-q); q = 4/7 is SNR 4/3.

## Claims and where they are proved

| ID | Claim | Location in [manuscript.md](manuscript.md) | Evidence |
|---|---|---|---|
| C1 | Theorem A: for Q < 4/7, explicit delta_Q with F'' <= -(1/16) ||N(T_sqrt(qt) g + T_sqrt(q(1-t)) g)||^2 for all q <= Q, t in [0,1] | Section 2.3; Part I | proof |
| C2 | Lemma I.2: ||T_sqrt(s) g||^2 <= J ||g||_inf ||T_sqrt(r) g|| for r < s < b(r) = sqrt(2r(1+r)) - r | Section I.2 | proof (extends Luo F.1 to bounded inputs) |
| C3 | Corollary: (0,4/7) subset C subset (0,4/7], sup C = 4/7; same in relative C_b^m neighbourhoods and log-concave moment-matched classes | Section 2.3 | proof + Luo (2026) Thm 1.2 |
| C4 | Theorem B: Luo's moment-corrected sequence has F''(1/2) > 0 eventually for every q in (4/7,1]; F''(1/2) ~ 4 n gamma(q) eps_n^4 ||M||^2, gamma = 21 sqrt5 q/(49q-4); profile 1 - exp(-4 gamma s^2) on scale n^(-1/2) | Section 2.3; Part II | proof (sign on (4/7, 2-sqrt2] also proposed earlier in [8]) |
| C5 | Proposition S: for s in (b(r),1] the interpolation quotient is unbounded (also for admissible g); envelope (S4); extremal Hermite parameter a = (1-r^2)^(-1/2), = 7/(3 sqrt5) at r = 2/7 | Section III.1 | proof |
| C6 | Corollary S.5(4): for q <= 4/7 no Hermite family h_{n,a} has the exponential gap (M16) | Section III.1.4 | proof (rates) |
| C7 | Proposition E: at q = 4/7, Luo's family with eps_n -> 0, eps_n n^(1/4) -> 0 has F'' <= -c eps^2 n^2 ||L_t||^2 < 0 on a neighbourhood of t = 1/2; midpoint constant 136/(15 sqrt(2(1+a))) | Section III.2 | proof |
| C8 | Remark R: delta_Q ~ C_delta (4/7-Q)^(9/4), C_delta ~ 5.92873e-3 | Section III.3 | first-order expansion; numerical confirmation |
| C9 | Blind-source-separation reading (no interior minimum below SNR 4/3; spurious 45-degree minimum above) | Discussion D.3 | corollary of C1, C4 |

Not claimed: membership of 4/7 in C; any global (non-local) concavity result; asymmetric inputs (they fail for every
q > 0, Luo Thm 1.1); optimality of delta_Q or of 1/16; behaviour of Luo's family away from t = 1/2 at q = 4/7.

## Evidence files

- [checks/coefficient_checks.py](checks/coefficient_checks.py), [checks/local_positive_checks.py](checks/local_positive_checks.py),
  [checks/identity_checks.py](checks/identity_checks.py), [checks/extension_checks.py](checks/extension_checks.py),
  [checks/negative_controls.py](checks/negative_controls.py); outputs in [checks/results/](checks/results/)
  (normal and `-O` runs). Fail-closed; standard library only. See [VERIFICATION.md](VERIFICATION.md).
- [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md): internal AI-assisted review objections and dispositions.
- [RESEARCH_GATES.json](RESEARCH_GATES.json): contribution map, sources read, limitation triage.
- [boundary_status.md](boundary_status.md): what is known at q = 4/7.
- [refs.bib](refs.bib): bibliography copied from registry records.
- [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md): release replay; [MANIFEST.sha256](MANIFEST.sha256): SHA-256 of every file;
  [CITATION.cff](CITATION.cff) and [.zenodo.json](.zenodo.json): citation metadata (DOI 10.5281/zenodo.23910629).
- [history/manuscript_v1_2026-11-02.md](history/manuscript_v1_2026-11-02.md): superseded first internal version.

## Limits

No human peer review, independent human reproduction or formal verification. Numerical evaluations are not interval
certificates. The cited 2026 sources (Luo, Schubr, Abu-Ajamieh, the GitHub proposal) are unrefereed.

# Review summary (public, rating-free)

Internal AI-assisted review, 3 November 2026. Not human peer review. Reviewers: five independent Claude Code reviewer
agents (editor, methodology/proof, domain, cross-disciplinary perspective, devil's advocate); one adversarial review
and one proof audit by a different model (Astra via the Codex companion); a separate proof-audit agent and a
re-review agent for the revision; and a clean replay of every shipped script by the editor.

## First round: verdicts

No reviewer found a mathematical error or a counterexample. All load-bearing constants were recomputed
independently: the condition b(q/2) > q iff q < 4/7, the radius table, gamma(q), the Hessian determinant, (M20)-(M21)
and the Lemma II.5 matrix. Requested changes concerned positioning, presentation and reproducibility.

## Main objections and dispositions

| Objection | Disposition |
|---|---|
| An earlier public proposal for Luo's Problem 10.3 (GitHub issue the-omega-institute/trureturing #25275, 1 Nov 2026) used the same construction and a differentiated saddle mechanism, and was not cited | Cited as [8] throughout. Theorem B's novelty restated: the sign on (4/7, 2 - sqrt2] was proposed there; new here are the exact constant, the local profile, one sequence for all q, a single argument covering q = 1, and a complete written proof |
| What is new versus Luo (Appendices E, F) and Schubr was unclear; the new range of Theorem B's sign is only about 0.014 wide | New table in Section 1.2; explanation of why Luo's interpolation does not already give Theorem A; regimes table separates Luo's range q > 2 - sqrt2 |
| Why 4/7? The two thresholds looked coincidental | New Section 1.3 and Proposition S: the smoothing inequality is sharp at s = b(r); 4/7 is the fixed point q = b(q/2); within the Hermite test family the extremal member is Luo's counterexample family |
| "Critical" depends on the topology; Abu-Ajamieh's 5/8 for a fixed direction looked inconsistent | The neighbourhood is named everywhere; Section 1.3 and Discussion D.1 reconcile the fixed He_4 direction |
| Motivation and context missing for non-specialists; thin bibliography | New abstract and introduction (conjecture, Lieb's inequality, Gaussian-mixture results, Luo's refutation); 18 references copied from Crossref, arXiv or GitHub records |
| Two stapled parts with restarted numbering and duplicated definitions | Single setup section; unique section labels; roadmaps; figure added |
| Several compressed steps (chain rule before (P20); smoothing margin; regularity in Lemma II.5; (M29); (M33)) | Written out |
| Shipped scripts relied on `assert`, hard-coded a pass flag and never compared key values with their targets | Rewritten fail-closed with target comparisons; new identity and extension checkers; negative controls (16/16 caught) |
| The boundary q = 4/7 | New Proposition E: Luo's family is strictly concave near equal weighting at q = 4/7; the boundary remains open |
| Practical size; interpretation outside the subfield | Leading-order dip sizes stated; Discussion D.3 gives a blind-source-separation reading with caveats |
| Suggested but not done: Rényi-entropy and k-copy versions, degree-resolved thresholds, neighbourhood-topology dependence, a better radius | Listed as open problems (D.5) |

## Verification round

The re-review confirmed that the fixes were made and found no mathematical error. The proof audits found no critical
issue. They found two overstatements (the scope of Proposition E, which holds only near t = 1/2, and "extremal inputs",
which holds only within the Hermite test family), a missing lower bound in Corollary S.5(4), and several local slips
(a margin constant, two exponents or identities in Lemmas S.1 and S.3, an asymmetric bilinear form treated as
symmetric in Lemma E.5, a midpoint identity in E.4). All were corrected; none changed a stated result. A missing blank
line that hid the Part III heading in the PDF was fixed, and the shipped checks folder was replaced by the final
fail-closed set.

#!/usr/bin/env python3
"""Fail-closed diagnostics for Part I (local entropy concavity below 4/7).

Standard library only (decimal at 100 digits, fractions). Every check goes
through ``Checker.check``; the pass flag is computed from the recorded results
and any failure or uncaught exception ends the run with exit status 1. There is
no ``assert`` statement, so ``python -O`` changes nothing.

Manuscript checks (group "manuscript"):
* Lemma 2.1 closed-form recipe (P5), (P7), (P9), (P10) and the radius (P21) at
  q = 0.50, 0.53, 0.56, 0.57, including positivity of E_1, E_2, epsilon, the
  strict kernel-comparison gap (P2), P > beta_B, 2 alpha^eps > beta, the
  ordering r < q < s < b(r) < 1, the direct square root (P9) against the
  Mehler parametrisation and direct Gaussian composition, and agreement of
  s, epsilon, gap, J and the radius with the Part I Section I.7 table to the
  displayed digits;
* the threshold identity (3q/2)^2 - q(1+q/2) = q(7q-4)/4 (Part I §4).

Groups excluded from the manuscript claim list (still checked, fail-closed):
* "auxiliary_not_used_by_manuscript": endpoint constants of the superseded
  half-signal note (928, 133, 8126464, 1572/30000), and the alternative
  dyadic epsilon recipe;
* "diagnostic_only": the cosine interpolation ratios (non-adversarial) and the
  auxiliary L2 bilinear sharpness table (unbounded e_n, not admissible).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from decimal import Decimal, getcontext
from fractions import Fraction
from math import comb
from pathlib import Path

getcontext().prec = 110

AUX = "auxiliary_not_used_by_manuscript"
DIAG = "diagnostic_only"

# Part I, Section I.7 table, as displayed in manuscript.md (closed-form recipe (P7)).
SECTION7_TABLE = {
    "0.50": dict(s="0.5202847075", eps="0.1375827353", gap="0.0473458334", J="1.4320544662", radius="1.7134e-5"),
    "0.53": dict(s="0.5419050561", eps="0.0870079716", gap="0.0296201831", J="1.5902847695", radius="4.8520e-6"),
    "0.56": dict(s="0.5633202098", eps="0.0262722883", gap="0.0088148014", J="2.1257516652", radius="2.5717e-7"),
    "0.57": dict(s="0.5704164638", eps="0.0033875871", gap="0.0011301518", J="3.5370963892", radius="2.3569e-9"),
}


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
            excluded_from_manuscript_claims=sorted({r["group"] for r in self.results} - {"manuscript"}),
            all_checks_passed=(not failed) and bool(self.results),
        )


def half_ulp(displayed):
    """Half a unit in the last displayed digit of a decimal string."""
    d = Decimal(displayed)
    return Decimal(1).scaleb(d.as_tuple().exponent) / 2


def matches_display(value, displayed):
    return abs(value - Decimal(displayed)) <= half_ulp(displayed)


def arctangent_reciprocal(n):
    x = Decimal(1) / n
    x_squared = x * x
    power = x
    total = x
    sign = -1
    for k in range(1, 10000):
        power *= x_squared
        term = power / (2 * k + 1)
        total += sign * term
        if abs(term) < Decimal("1e-108"):
            return total
        sign = -sign
    raise ArithmeticError("Machin series failed to converge")


class DecimalMath:
    pi = 16 * arctangent_reciprocal(5) - 4 * arctangent_reciprocal(239)
    mpf = staticmethod(lambda x: Decimal(str(x)))
    sqrt = staticmethod(lambda x: Decimal(x).sqrt())
    exp = staticmethod(lambda x: Decimal(x).exp())
    power = staticmethod(lambda x, y: Decimal(x) ** y)
    nstr = staticmethod(lambda x, digits: format(x, f".{digits}g"))


mp = DecimalMath()


def geometric_moment(k, r):
    """Sum_d d**k r**d using Stirling numbers of the second kind."""
    stirling = [[0] * (k + 1) for _ in range(k + 1)]
    stirling[0][0] = 1
    for n in range(1, k + 1):
        for j in range(1, n + 1):
            stirling[n][j] = j * stirling[n - 1][j] + stirling[n - 1][j - 1]
    ans = Fraction(0)
    factorial = 1
    for j in range(k + 1):
        if j > 0:
            factorial *= j
        ans += stirling[k][j] * factorial * r**j / (1 - r) ** (j + 1)
    return ans


def endpoint_exact_checks(ck):
    """Superseded half-signal (q=1/2) note; NOT used by the current manuscript."""
    ref = "sources/half_signal_supplement.md (superseded; not in manuscript)"
    moments = [geometric_moment(k, Fraction(3, 4)) for k in range(5)]
    sum_b_squared = (moments[4] + 2 * moments[3] + 9 * moments[2] + 8 * moments[1] + 16 * moments[0]) / 16
    sum_n_kernel = 16 * moments[0] + 4 * moments[1] + moments[2] / 4
    sum_tt_kernel = 8192 * (16 * moments[0] + sum_b_squared)
    ck.check(moments == [4, 12, 84, 876, 12180], "half-signal geometric moments 4,12,84,876,12180",
             [str(m) for m in moments], group=AUX, ref=ref)
    ck.check(sum_b_squared == 928, "half-signal sum B_d^2 = 928", sum_b_squared, group=AUX, ref=ref)
    ck.check(sum_n_kernel == 133, "half-signal number-operator kernel sum = 133", sum_n_kernel, group=AUX, ref=ref)
    ck.check(sum_tt_kernel == 8126464, "half-signal second-derivative kernel sum = 8126464", sum_tt_kernel,
             group=AUX, ref=ref)
    ck.check(sum_n_kernel < 12**2, "half-signal 133 < 12^2", "", group=AUX, ref=ref)
    ck.check(sum_tt_kernel < 3000**2, "half-signal 8126464 < 3000^2", "", group=AUX, ref=ref)
    ck.check(Fraction(1572, 30000) < Fraction(1, 16), "half-signal 1572/30000 < 1/16", "", group=AUX, ref=ref)
    return {
        "scope": "auxiliary_not_used_by_manuscript: constants of the superseded half-signal note; "
                 "excluded from the manuscript claim list.",
        "geometric_moments": [str(x) for x in moments],
        "sum_B_d_squared": str(sum_b_squared),
        "number_operator_kernel_sum": str(sum_n_kernel),
        "second_derivative_kernel_sum": str(sum_tt_kernel),
        "endpoint_radius": "1/30000",
        "curvature_error_coefficient": "1572/30000 < 1/16",
    }


def positive_parameters(ck, q_string, recipe):
    manuscript = recipe == "closed_form"
    group = "manuscript" if manuscript else AUX
    tag = f"q={q_string} [{recipe}]"
    q = mp.mpf(q_string)
    r = q / 2
    boundary = mp.sqrt(2 * r * (1 + r)) - r
    s = (q + boundary) / 2
    alpha = lambda u: 1 / (2 * (1 - u * u))
    beta = lambda u: u / (1 - u * u)
    c = lambda u: 1 / (2 * mp.pi * mp.sqrt(1 - u * u))
    if manuscript:
        ck.check(0 < r < q < s < boundary < 1, f"{tag}: 0<r<q<s<b(r)<1 (P5), Part I §4",
                 (r, s, boundary), ref="Part I (P5), §4")
        ck.check(abs((3 * q / 2) ** 2 - q * (1 + q / 2) - q * (7 * q - 4) / 4) < Decimal("1e-95") and q * (7 * q - 4) < 0,
                 f"{tag}: b(q/2)>q iff q<4/7 via (3q/2)^2-q(1+q/2)=q(7q-4)/4", "", ref="Part I §4")
        P0 = 1 / (1 - r)
        B1 = 1 / (1 - s) - P0
        B2 = P0 - 1 / (1 + s)
        E1 = P0 - B1 * B1 / (2 * beta(r))
        E2 = P0 - B2 * B2 / (2 * beta(r))
        ck.check(B1 > 0 and B2 > 0, f"{tag}: B_1,B_2>0", (B1, B2), ref="Part I §2")
        ck.check(E1 > 0 and E2 > 0, f"{tag}: E_1,E_2>0 (positivity of (P7))", (E1, E2), ref="Part I (P7)")
        cc = mp.sqrt(2 * r / (1 + r))
        ck.check(abs(mp.sqrt(2 * beta(r) * P0) - P0 * cc) < Decimal("1e-95"),
                 f"{tag}: sqrt(2 beta P_0) = P_0 c", "", ref="Part I §2")
        ck.check(abs((r + cc) / (1 + cc) - boundary) < Decimal("1e-95"),
                 f"{tag}: (r+c)/(1+c) = sqrt(2r(1+r)) - r", "", ref="Part I §2")
        epsilons = [min(1 / (1 + r), 1 / (1 + s), E1, E2) / 2]
    else:
        eps0 = min(1 / (2 * (1 + r)), 1 / (2 * (1 + s)))
        epsilons = [eps0 * mp.power(2, -m) for m in range(500)]
    found = False
    for m, eps in enumerate(epsilons):
        alpha_r = alpha(r) - eps / 2
        beta_r = beta(r)
        a = mp.sqrt(4 * alpha_r**2 - beta_r**2)
        lam = beta_r / (2 * alpha_r + a)
        v = mp.sqrt(lam)
        eig0 = c(r) * mp.sqrt(2 * mp.pi / (2 * alpha_r + a))
        alpha_b = a * (1 + v * v) / (2 * (1 - v * v))
        beta_b = 2 * a * v / (1 - v * v)
        c_b = mp.sqrt(eig0) * mp.sqrt(a / mp.pi) / mp.sqrt(1 - v * v)
        gap = beta_b - beta(s) - 2 * abs(alpha(s) - eps / 2 - alpha_b)
        if gap > 0:
            found = True
            break
    ck.check(found and eps > 0, f"{tag}: epsilon>0 with strict kernel gap (P2) > 0", gap, group=group,
             ref="Part I (P2),(P7)")
    alpha_composed = alpha_b - beta_b**2 / (8 * alpha_b)
    beta_composed = beta_b**2 / (4 * alpha_b)
    c_composed = c_b**2 * mp.sqrt(mp.pi / (2 * alpha_b))
    composition_error = max(abs(alpha_composed / alpha_r - 1), abs(beta_composed / beta_r - 1),
                            abs(c_composed / c(r) - 1))
    ck.check(composition_error < mp.mpf("1e-85"), f"{tag}: Mehler kernel composition B*B = A_r^eps",
             composition_error, group=group, ref="Part I (P9)")
    P = 1 / (1 - r) - eps
    direct_alpha_b = P / 2
    direct_beta_b = mp.sqrt(2 * beta(r) * P)
    direct_c_b = mp.sqrt(c(r)) * (P / mp.pi) ** mp.mpf("0.25")
    sq_err = max(abs(direct_alpha_b / alpha_b - 1), abs(direct_beta_b / beta_b - 1), abs(direct_c_b / c_b - 1))
    ck.check(sq_err < mp.mpf("1e-85"), f"{tag}: direct square root (P9) equals Mehler parametrisation",
             sq_err, group=group, ref="Part I (P9)")
    ck.check(P > direct_beta_b > 0 and eps < 1 / (1 + r), f"{tag}: P > beta_B > 0", (P, direct_beta_b),
             group=group, ref="Part I (P9)")
    ck.check(all(2 * (alpha(u) - eps / 2) > beta(u) > 0 for u in (r, s)), f"{tag}: 2 alpha_u^eps > beta_u > 0, u=r,s",
             "", group=group, ref="Part I §2")
    if manuscript:
        Pb, Bb = P, direct_beta_b
        ck.check(B1 < Bb and B2 < Bb, f"{tag}: B_1,B_2 < beta_B (strict (P2))", (B1, B2, Bb), ref="Part I §2")
        lhs1 = 2 * (alpha(s) - eps / 2) + beta(s) < Pb + Bb
        lhs2 = Pb - Bb < 2 * (alpha(s) - eps / 2) - beta(s)
        ck.check(lhs1 and lhs2, f"{tag}: (P2) in the form 2a_s^eps +- b_s vs P +- beta_B", "", ref="Part I §2")
    J = c(s) / direct_c_b * direct_beta_b / beta(s) * (mp.pi / eps) ** mp.mpf("0.25")
    rho = q / s
    eta = (1 / rho - 1) / 4
    R = (1 + rho) / 2
    K = mp.sqrt(R) / (1 - mp.sqrt(R)) ** 2
    A = K * J / 8 + 60 + J / (8 * eta**2)
    terms = [mp.mpf(1) / 12, 1 / J, 8 / (K * J), 1 / (16 * A)]
    delta = min(terms)
    ck.check(delta > 0 and 0 < rho < 1 and R < 1, f"{tag}: delta_q>0 with rho=q/s<1 and R_rho<1 (P21)", delta,
             group=group, ref="Part I (P21)")
    if manuscript:
        ck.check(terms.index(delta) == 3, f"{tag}: binding term of (P21) is 1/(16 A_q)", terms.index(delta),
                 ref="Part I (P21)")
        row = SECTION7_TABLE[q_string]
        for key, value in (("s", s), ("eps", eps), ("gap", gap), ("J", J), ("radius", delta)):
            ck.check(matches_display(value, row[key]),
                     f"{tag}: Section I.7 table {key} = {row[key]} to displayed digits",
                     f"computed {mp.nstr(value, 15)}", ref="Part I §7 table")
    # Cosine family g = cos(kx): exact smoothing norm; non-adversarial (diagnostic only).
    cosine_ratios = []
    for frequency in [0, 1, 2, 5, 10, 20]:
        kk = mp.mpf(frequency) ** 2
        norm_s_squared = mp.exp(-(1 - s) * kk) * (1 + mp.exp(-2 * s * kk)) / 2
        norm_r_squared = mp.exp(-(1 - r) * kk) * (1 + mp.exp(-2 * r * kk)) / 2
        ratio = norm_s_squared / (J * mp.sqrt(norm_r_squared))
        ck.check(ratio <= 1, f"{tag}: cosine k={frequency} consistent with (P6)", ratio, group=DIAG,
                 ref="Part I (P6) (non-adversarial sanity test)")
        cosine_ratios.append(mp.nstr(ratio, 20))
    return {
        "q": q_string, "recipe": recipe,
        "scope": "manuscript (P7) recipe; Section I.7 table" if manuscript else
                 "auxiliary_not_used_by_manuscript: alternative dyadic epsilon search",
        "r": mp.nstr(r, 30), "s": mp.nstr(s, 30), "epsilon": mp.nstr(eps, 30),
        "dyadic_epsilon_index": None if manuscript else m,
        "strict_kernel_gap": mp.nstr(gap, 30),
        "Mehler_kernel_composition_relative_error": mp.nstr(composition_error, 8),
        "direct_square_root_vs_Mehler_relative_error": mp.nstr(sq_err, 8),
        "J": mp.nstr(J, 30), "K": mp.nstr(K, 30), "eta": mp.nstr(eta, 30), "A_q": mp.nstr(A, 30),
        "sufficient_radius": mp.nstr(delta, 30),
        "cosine_interpolation_ratios_diagnostic_only": cosine_ratios,
    }


def auxiliary_bilinear_sharpness(ck):
    """Diagnostic only: L2 operator ratios for unbounded e_n (not admissible g)."""
    degrees = [4, 10, 20, 50, 100, 200, 1000]
    rows = []
    for label in ["0.49", "0.50", "0.51", "0.56", "4/7", "0.58"]:
        q = Decimal(4) / 7 if label == "4/7" else Decimal(label)
        ratios = [2 * Decimal(comb(2 * n, n)).sqrt() * (q / 2) ** (n // 2) for n in degrees]
        rows.append({"q": label, "ratios": [mp.nstr(x, 20) for x in ratios]})
        # Internal consistency of the stated formula: ~ (2q)^(n/2), threshold q = 1/2.
        if q < Decimal("0.5"):
            ck.check(ratios[-1] < 1, f"bilinear sharpness q={label}: ratio(1000)<1 below 1/2", ratios[-1], group=DIAG)
        elif q > Decimal("0.5"):
            ck.check(ratios[-1] > 1, f"bilinear sharpness q={label}: ratio(1000)>1 above 1/2", ratios[-1], group=DIAG)
    return {
        "scope": "diagnostic_only: auxiliary L2 operator sharpness for the unbounded e_n, which is not an "
                 "admissible bounded relative-density perturbation; not used by the manuscript.",
        "formula": "2 * sqrt(binomial(2*n,n)) * (q/2)**(n/2)",
        "even_degrees": degrees,
        "values": rows,
    }


def run(ck, result):
    result["weighted_kernel_checks_closed_form_P7"] = [
        positive_parameters(ck, q, "closed_form") for q in ["0.50", "0.53", "0.56", "0.57"]]
    ck.check(len(result["weighted_kernel_checks_closed_form_P7"]) == 4, "four Section I.7 table rows evaluated",
             ref="Part I §7 table")
    result["auxiliary_not_used_by_manuscript"] = {
        "endpoint_checks_half_signal_note": endpoint_exact_checks(ck),
        "dyadic_epsilon_recipe": [positive_parameters(ck, q, "dyadic") for q in ["0.50", "0.53", "0.56", "0.57"]],
    }
    result["diagnostic_only"] = {"auxiliary_bilinear_sharpness": auxiliary_bilinear_sharpness(ck)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("local_positive_output"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    getcontext().prec = 100
    start = time.perf_counter()
    ck = Checker()
    result = {
        "scope": "Exact rational constants and 100-decimal diagnostic checks; no claim of formal verification.",
        "python_optimize_flag": sys.flags.optimize,
    }
    try:
        run(ck, result)
    except Exception as exc:
        ck.check(False, "uncaught exception", repr(exc) + "\n" + traceback.format_exc())
    result.update(ck.summary())
    result["checks"] = ck.results
    result["elapsed_seconds"] = time.perf_counter() - start
    text = json.dumps(result, indent=2) + "\n"
    (args.output_dir / "local_positive_checks.json").write_text(text, encoding="utf-8")
    lines = [f"{'PASS' if r['passed'] else 'FAIL'} [{r['group']}] {r['label']}" for r in ck.results]
    status = "All checks passed" if result["all_checks_passed"] else f"{result['number_failed']} CHECK(S) FAILED"
    lines.append(f"{status}: {result['number_of_checks']} checks "
                 f"({result['manuscript_claim_checks']} manuscript-claim checks) in {result['elapsed_seconds']:.3f} s.")
    (args.output_dir / "local_positive_checks.log").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(lines[-1], flush=True)
    sys.exit(0 if result["all_checks_passed"] else 1)


if __name__ == "__main__":
    main()

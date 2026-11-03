#!/usr/bin/env python3
"""Fail-closed numerical checks for Part III of manuscript.md (Propositions S, E, Remark R).

Python standard library only (decimal, fractions, math, json, argparse).
No assert statements: every check goes through check(), failures are recorded,
all_passed is computed, JSON is written to --output-dir, and the process exits
with status 1 on any failure (also under python -O).

Everything here is a numerical diagnostic or an exact-identity test at finite
size. Nothing in this file is part of the proofs in Part III of manuscript.md.
"""
import argparse
import json
import math
import os
import sys
import time
from decimal import Decimal as D, getcontext
from fractions import Fraction

getcontext().prec = 60

RESULTS = []
CONTROLS = []


def fmt(x):
    if isinstance(x, D):
        return format(x, '.12e')
    return str(x)


def check(cond, label, **data):
    ok = bool(cond)
    RESULTS.append({"label": label, "passed": ok,
                    "data": {k: fmt(v) for k, v in data.items()}})
    return ok


def negative_control(detected, label, **data):
    """A deliberately wrong target. 'detected' must be True (the wrong
    target was caught). Being caught does not fail the run; NOT catching
    it is recorded as a failure of the checker itself."""
    caught = bool(detected)
    CONTROLS.append({"label": label, "caught": caught,
                     "data": {k: fmt(v) for k, v in data.items()}})
    RESULTS.append({"label": "negative control caught: " + label,
                    "passed": caught, "data": {}})
    return caught


ONE = D(1)
TWO = D(2)
SQ5 = D(5).sqrt()


def dpow(x, y):
    return (D(y) * D(x).ln()).exp()


def fact(n):
    return D(math.factorial(n))


# ---------------------------------------------------------------- rates ----
def tau(a, u):
    return (1 - a + a * u) / (1 + a - a * u)


def nu(a, u):
    return (a - 1 + a * u) / (1 + a + a * u)


def rate(a, u):
    return max(abs(tau(a, u)), abs(nu(a, u)))


def b_of(r):
    return (2 * r * (1 + r)).sqrt() - r


def b_closed(r):
    lam = r / (1 + (1 - r * r).sqrt())
    w = lam.sqrt()
    return 2 * w * (1 - w + w * w) / (1 + w ** 4)


def phi_w(w):
    return w / (1 + w + w * w)


# -------------------------------------------- Hermite coefficients of h_n --
def hermite_coeffs(n, a, kmax):
    """Standard (probabilists') orthonormal Hermite coefficients of
    h_{n,a}(x)=sqrt2 exp(-a x^2/2) He_n(sqrt(2a)x)/sqrt(n!), k=0..kmax,
    from the generating function
    sqrt2 (1+a)^(-1/2) exp(b z^2/2 + g z w - c w^2/2)."""
    bb = (a - 1) / (a + 1)
    g = (2 * a).sqrt() / (1 + a)
    c = a / (1 + a)
    pre = TWO.sqrt() / (1 + a).sqrt()
    zc = {}
    for p in range(0, n // 2 + 1):
        zc[2 * p] = (bb / 2) ** p / fact(p)
    out = []
    for k in range(kmax + 1):
        s = D(0)
        for m in range(0, min(n, k) + 1):
            if (n - m) % 2 or (k - m) % 2:
                continue
            p2 = (k - m) // 2
            wc = (-c / 2) ** p2 / fact(p2)
            s += g ** m / fact(m) * zc[n - m] * wc
        out.append(pre * (fact(n) * fact(k)).sqrt() * s)
    return out


def cbin(j):
    return D(math.comb(2 * j, j)) / D(4) ** j


def P_n(a, u, n):
    t, v = tau(a, u), nu(a, u)
    return sum(cbin(j) * cbin(n - j) * t ** j * v ** (n - j) for j in range(n + 1))


def Pi_(a, u):
    return (1 + a) ** 2 - a * a * u * u


# ------------------------------------------------------------ radius (P21) -
def radius(Q):
    q = D(Q)
    r = q / 2
    br = b_of(r)
    s = (q + br) / 2
    P0 = 1 / (1 - r)
    beta = r / (1 - r * r)
    B1 = 1 / (1 - s) - P0
    B2 = P0 - 1 / (1 + s)
    E1 = P0 - B1 ** 2 / (2 * beta)
    E2 = P0 - B2 ** 2 / (2 * beta)
    eps = min(1 / (1 + r), 1 / (1 + s), E1, E2) / 2
    pi = D("3.14159265358979323846264338327950288419716939937510582097494")

    def c(u):
        return 1 / (2 * pi * (1 - u * u).sqrt())
    P = P0 - eps
    bB = (2 * beta * P).sqrt()
    cB = c(r).sqrt() * dpow(P / pi, D("0.25"))
    bs = s / (1 - s * s)
    J = c(s) / cB * bB / bs * dpow(pi / eps, D("0.25"))
    rho = q / s
    eta = (1 / rho - 1) / 4
    R = (1 + rho) / 2
    K = R.sqrt() / (1 - R.sqrt()) ** 2
    A = K * J / 8 + 60 + J / (8 * eta ** 2)
    d = min(ONE / 12, 1 / J, 8 / (K * J), 1 / (16 * A))
    return {"s": s, "eps": eps, "J": J, "delta": d, "pi": pi}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default="./extension_check_output")
    args = ap.parse_args()
    t0 = time.time()
    summary = {}

    # ======================= Prop S: envelope =============================
    env = []
    A_grid = [D(k) / D(500) for k in range(-4000, 1501)]   # log a in [-8,3]
    a_grid = [A.exp() for A in A_grid]

    def maxratio(r, s):
        best, ab = None, None
        for a in a_grid:
            v = rate(a, s) ** 2 / rate(a, r)
            if best is None or v > best:
                best, ab = v, a
        A = ab.ln()
        h = D(1) / D(500)
        for _ in range(6):
            h = h / 20
            for k in range(-40, 41):
                a = (A + k * h).exp()
                v = rate(a, s) ** 2 / rate(a, r)
                if v > best:
                    best, ab = v, a
            A = ab.ln()
        ar = 1 / (1 - r * r).sqrt()
        v = rate(ar, s) ** 2 / rate(ar, r)
        if v > best:
            best, ab = v, ar
        return best, ab

    for rs in ["0.05", "0.1", "0.2", "2/7", "0.4", "0.6", "0.9"]:
        r = D(Fraction(rs).numerator) / D(Fraction(rs).denominator)
        br = b_of(r)
        ar = 1 / (1 - r * r).sqrt()
        check(abs(br - b_closed(r)) < D("1e-50"),
              "S: b(r) equals closed form 2w(1-w+w^2)/(1+w^4), r=" + rs)
        lam = r / (1 + (1 - r * r).sqrt())
        check(abs(tau(ar, r) - lam) < D("1e-50") and abs(nu(ar, r) - lam) < D("1e-50"),
              "S: tau=nu=lambda at a_r, r=" + rs)
        # region-I criterion at a_r equals the envelope equation
        crit = (br - r) / (2 - r - br) - phi_w(lam.sqrt())
        check(abs(crit) < D("1e-50"), "S: (s-r)/(2-r-s)=phi(sqrt lambda) at s=b(r), r=" + rs)
        m0, a0 = maxratio(r, br)
        mlo, _ = maxratio(r, br - D("0.001"))
        mhi, ahi = maxratio(r, br + D("0.001"))
        env.append({"r": rs, "b": fmt(br), "max_at_b": fmt(m0), "argmax_at_b": fmt(a0),
                    "a_r": fmt(ar), "max_below": fmt(mlo), "max_above": fmt(mhi),
                    "argmax_above": fmt(ahi)})
        check(m0 <= 1 + D("1e-40") and m0 >= 1 - D("1e-40"),
              "S: sup_a rate_s^2/rate_r = 1 at s=b(r) (grid+refine), r=" + rs, m=m0)
        check(abs(a0 - ar) < D("1e-6"), "S: maximiser at s=b(r) is a_r=(1-r^2)^(-1/2), r=" + rs, a0=a0, ar=ar)
        check(mlo < 1, "S: sup_a ratio < 1 for s=b(r)-1e-3, r=" + rs, m=mlo)
        check(mhi > 1, "S: sup_a ratio > 1 for s=b(r)+1e-3 (divergence), r=" + rs, m=mhi)
        # s=1 case
        check(rate(ar, ONE) ** 2 / rate(ar, r) > 1, "S: s=1 diverges, r=" + rs)
        # sanity of region analysis: criterion values at sample a
        for A in [D(-3), D(-1), ar.ln() / 2, ar.ln(), ar.ln() * 3 / 2 + D("0.01"), D(1), D(2)]:
            a = A.exp()
            for s in [br - D("0.01"), br, br + D("0.01")]:
                if s >= 1:
                    continue
                ar_s = 1 / (1 - s * s).sqrt()
                Rr, Rs = rate(a, r), rate(a, s)
                direct = Rs ** 2 > Rr
                if a <= ar:
                    td = (s - r) / (2 - r - s)
                    crit2 = td > phi_w(Rr.sqrt())
                    check(direct == crit2, "S: region I criterion matches direct, r=%s a=%s" % (rs, fmt(a)))
                elif a >= ar_s:
                    td = (s - r) / (2 + r + s)
                    crit2 = td > phi_w(Rr.sqrt())
                    check(direct == crit2, "S: region III criterion matches direct, r=%s a=%s" % (rs, fmt(a)))
    summary["envelope"] = env

    # ======================= r=2/7 coincidences ===========================
    a7 = D(7) / (3 * SQ5)
    r7 = D(2) / D(7)
    q7 = D(4) / D(7)
    kap = (3 - SQ5) / 2
    check(abs(1 / (1 - r7 * r7).sqrt() - a7) < D("1e-55"), "C: a_{2/7}=7/(3 sqrt5)")
    check(abs(tau(a7, r7) - kap ** 2) < D("1e-55"), "C: tau_a(2/7)=kappa^2")
    check(abs(nu(a7, r7) - kap ** 2) < D("1e-55"), "C: nu_a(2/7)=kappa^2")
    check(abs(kap ** 2 - (7 - 3 * SQ5) / 2) < D("1e-55"), "C: kappa^2=(7-3sqrt5)/2")
    check(abs(tau(a7, q7) - kap) < D("1e-55"), "C: tau_a(4/7)=kappa")
    lam_m = (a7 * (1 - q7) - 1) / (a7 * (1 - q7) + 1)
    check(abs(-lam_m - kap) < D("1e-55"), "C: kappa=-lambda_-(4/7) as in (M13)")
    check(abs(b_of(r7) - q7) < D("1e-55"), "C: b(2/7)=4/7 (fixed point q=b(q/2))")
    for qq in ["0.3", "0.5", "0.57"]:
        qd = D(qq)
        check(b_of(qd / 2) > qd, "C: b(q/2)>q for q=" + qq)
    for qq in ["0.58", "0.7", "1"]:
        qd = D(qq)
        check(b_of(qd / 2) < qd, "C: b(q/2)<q for q=" + qq)
    check(abs(Pi_(a7, r7) - 2 * (1 + a7)) < D("1e-55"), "C: Pi_a(2/7)=2(1+a)")
    check(abs((a7 - 1) / (a7 + 1) - kap ** 4) < D("1e-55"), "C: b=kappa^4 (so sqrt b=kappa^2<kappa)")
    tp = 2 * a7 / (1 + a7 - a7 * r7) ** 2
    vp = 2 * a7 / (1 + a7 + a7 * r7) ** 2
    al = r7 * tp / kap ** 2
    be = r7 * vp / kap ** 2
    G0 = (3 * al * al + 3 * be * be + 2 * al * be) / 8
    check(abs(al - 3 / SQ5) < D("1e-55") and abs(be - SQ5 / 3) < D("1e-55"), "E: alpha=3/sqrt5, beta=sqrt5/3")
    check(abs((al + be) / 2 - a7) < D("1e-55"), "E: mu_1=(alpha+beta)/2=a")
    check(abs(G0 - D(17) / D(15)) < D("1e-55"), "E: Gamma_0=17/15")
    lead = 8 * G0 / Pi_(a7, r7).sqrt()
    summary["E_leading_constant_8G0_over_sqrtPi"] = fmt(lead)

    # ======================= exact identity (M11) ==========================
    KMAX = 260
    ident = []
    for n in range(0, 13):
        hc = hermite_coeffs(n, a7, KMAX)
        direct = sum(hc[k] ** 2 * r7 ** k for k in range(KMAX + 1))
        closed = 2 / Pi_(a7, r7).sqrt() * kap ** (2 * n)
        rel = abs(direct - closed) / closed
        tail = hc[KMAX] ** 2 * r7 ** KMAX / closed
        ident.append({"n": n, "direct": fmt(direct), "closed": fmt(closed), "rel": fmt(rel)})
        check(rel < D("1e-40") and tail < D("1e-60"), "M11: ||T_sqrt(2/7) h_n||^2 = 2 Pi^-1/2 kappa^(2n), n=%d" % n, rel=rel)
        if n == 6:
            wrong = closed * (1 + D(1) / (n + 1))
            negative_control(abs(direct - wrong) / wrong > D("1e-10"),
                             "wrong closed form 2Pi^-1/2 kappa^(2n)(1+1/(n+1)) rejected at n=6")
    summary["identity_M11_r27"] = ident
    # general a,u including nu<0 and tau<0
    for (aa, uu) in [("0.5", "0.3"), ("0.7", "0.5"), ("1.6", "0.2"), ("3", "0.1")]:
        a, u = D(aa), D(uu)
        for n in [3, 7, 10]:
            hc = hermite_coeffs(n, a, 320)
            direct = sum(hc[k] ** 2 * u ** k for k in range(321))
            closed = 2 / Pi_(a, u).sqrt() * P_n(a, u, n)
            check(abs(direct - closed) < D("1e-35") * (abs(closed) + D("1e-30")),
                  "M11 general: a=%s u=%s n=%d" % (aa, uu, n), direct=direct, closed=closed)
    # sup bound B_n and the n^{-1/4} lower bound on sup at x=0
    for n in [2, 6, 10, 20, 40]:
        # sample |h_n| on a grid (diagnostic for B_n)
        best = D(0)
        for i in range(0, 801):
            y = D(i) / 40   # y = sqrt(2a) x in [0,20]
            # He_n(y) by recurrence
            h0, h1 = ONE, y
            if n == 0:
                he = h0
            else:
                for m in range(1, n):
                    h0, h1 = h1, y * h1 - m * h0
                he = h1
            val = abs(TWO.sqrt() * (-y * y / 4).exp() * he / fact(n).sqrt())
            if val > best:
                best = val
        Bn = (2 * ONE.exp()).sqrt() * dpow(D(n + 1), D("0.25"))
        check(best <= Bn, "S: sampled sup|h_n| <= B_n=sqrt(2e)(n+1)^(1/4), n=%d" % n, sup=best, Bn=Bn)
        check(best <= TWO.sqrt() + D("1e-30"), "S: sampled sup|h_n| <= sqrt2 (consistent with Cramer/Indritz remark), n=%d" % n)

    # ======================= Prop E: <L,L''> closed form ===================
    q = q7
    r = r7
    a = a7
    Pi = Pi_(a, r)
    DA = 1 - a * a * (1 - r) ** 2
    DB = a * a * (1 + r) ** 2 - 1
    Aval = 2 * a * r / DA
    Bval = 2 * a * r / DB
    thA = Aval * (1 - 2 * a * a * r * (1 - r) / DA)
    thB = Bval * (1 - 2 * a * a * r * (1 + r) / DB)
    check(abs(Aval - al) < D("1e-50") and abs(Bval - be) < D("1e-50"), "E: A(r)=alpha, B(r)=beta")
    pi1 = a * a * r * r / Pi
    pi2 = 2 * a * a * r * r / Pi + 3 * a ** 4 * r ** 4 / Pi ** 2
    check(abs(thA) < D("1e-50") and abs(thB) < D("1e-50"), "E: theta A = theta B = 0 at r=2/7")
    G1 = (al * al + be * be) / 8 - al * be / 4 + (thA + thB) / 2 + pi1 * (al + be) - (al + be)
    G2 = pi2 - 2 * pi1
    check(abs(G1 - (D(17) / 15 - 7 * SQ5 / 5)) < D("1e-50"), "E: Gamma_1 = 17/15 - 7 sqrt5/5")
    check(abs(G2 - (D(47) / 30 - 7 * SQ5 / 10)) < D("1e-50"), "E: Gamma_2 = 47/30 - 7 sqrt5/10")
    negative_control(abs(G1 - (D(17) / 15 - 7 * SQ5 / 6)) > D("1e-10"), "wrong Gamma_1 (7sqrt5/6) rejected")

    def LLpp_closed(n):
        p1 = n * (al + be) / 2
        p2 = (al * al + be * be) * (3 * n * n + n) / 8 + al * be * (n * n - n) / 4 + D(n) / 2 * (thA + thB)
        S = 2 / Pi.sqrt() * kap ** (2 * n)
        return 4 * S * (p2 + 2 * pi1 * p1 + pi2 - 2 * pi1 - 2 * p1), S, p1 + pi1

    diag = []
    KM = 300
    for n in [2, 4, 6, 8, 10, 12, 16, 20, 24]:
        hc = hermite_coeffs(n, a, KM)
        S = sum(hc[k] ** 2 * r ** k for k in range(KM + 1))
        mean = sum(hc[k] ** 2 * k * r ** k for k in range(KM + 1)) / S
        LLpp = 4 * sum(hc[k] ** 2 * k * (k - 2) * r ** k for k in range(KM + 1))
        cf, Sc, meanc = LLpp_closed(n)
        check(abs(LLpp - cf) < D("1e-35") * cf, "E: <L,L''>(1/2) direct = closed form, n=%d" % n, d=LLpp, c=cf)
        poly_cf = 8 / Pi.sqrt() * kap ** (2 * n) * (G0 * n * n + G1 * n + G2)
        check(abs(LLpp - poly_cf) < D("1e-35") * cf, "E: <L,L''>(1/2)=8Pi^-1/2 kappa^2n (G0 n^2+G1 n+G2), n=%d" % n)
        check(abs(mean - meanc) < D("1e-35") * meanc, "E: exact mean degree n(alpha+beta)/2+a^2r^2/Pi, n=%d" % n)
        check(mean >= n * a, "E: mean degree >= mu_1 n, n=%d" % n)
        # derivative-power margin: ||N^2 L|| <= C n^2 ||L|| diagnostic
        N2 = (4 * sum(hc[k] ** 2 * D(k) ** 4 * r ** k for k in range(KM + 1))).sqrt()
        Ln = (4 * S).sqrt()
        # bilinear term M_{1/2}=B_q(1/2;h,h) via (P14), and its t-derivatives
        Mc, M1c, M2c = [], [], []
        sb = {}
        for N in range(0, KM + 1, 2):
            s0 = s1 = s2 = D(0)
            for j in range(0, N + 1, 2):
                k = N - j
                if j > KM or k > KM:
                    continue
                key = (N, j)
                w = D(math.comb(N, j)).sqrt()
                term = hc[j] * hc[k] * w
                s0 += term
                s1 += term * (j - k)
                s2 += term * ((j - k) ** 2 - 2 * N)
            fac = dpow(q, D(N) / 2) * dpow(D("0.5"), D(N) / 2)
            Mc.append(fac * s0)
            M1c.append(fac * s1)
            M2c.append(fac * s2)
        Mn = sum(x * x for x in Mc).sqrt()
        Mpp = sum(x * x for x in M2c).sqrt()
        Mp = sum(x * x for x in M1c).sqrt()
        # inner products in the even Hermite basis (index N=2i)
        Lcoef = [2 * hc[N] * dpow(r, D(N) / 2) if N <= KM else D(0) for N in range(0, KM + 1, 2)]
        Lppcoef = [Lcoef[i] * (2 * i) * (2 * i - 2) for i in range(len(Lcoef))]
        LMpp = sum(Lcoef[i] * M2c[i] for i in range(len(Lcoef)))
        MLpp = sum(Mc[i] * Lppcoef[i] for i in range(len(Lcoef)))
        MMpp = sum(Mc[i] * M2c[i] for i in range(len(Lcoef)))
        # (M23) asymptotic for ||M_{1/2}||^2
        lp = (a * (1 + q) - 1) / (a * (1 + q) + 1)
        bb = (a - 1) / (a + 1)
        detH0 = 2 * (1 - lp / kap) * (1 - bb / kap) ** 2
        Delta = (1 + a) * ((1 + a) ** 2 - a * a * q * q).sqrt()
        pi = radius("0.5")["pi"]
        M23 = 4 / (pi * Delta * detH0.sqrt()) * kap ** (2 * n) / n
        row = {"n": n, "S(2/7)": fmt(S), "<L,L''>": fmt(LLpp),
               "<L,L''>/(n^2 S)": fmt(LLpp / (n * n * S)),
               "||N^2L||/(n^2||L||)": fmt(N2 / (n * n * Ln)),
               "||M||^2/M23": fmt(Mn ** 2 / M23),
               "||M''||/(n^2||M||)": fmt(Mpp / (n * n * Mn)),
               "||M'||/(n||M||)": fmt(Mp / (n * Mn))}
        for name, eps in [("exp(-sqrt n)", (-D(n).sqrt()).exp()), ("n^(-1/2)", 1 / D(n).sqrt())]:
            T2 = eps ** 2 * LLpp
            T3 = eps ** 3 * (abs(LMpp) + abs(MLpp))
            T4 = eps ** 4 * abs(MMpp)
            Bn = (2 * ONE.exp()).sqrt() * dpow(D(n + 1), D("0.25"))
            etan = 2 * eps * Bn + eps * eps * Bn * Bn
            Rn = eps * Ln + eps * eps * Mn
            Rpp = eps * (4 * sum(hc[k] ** 2 * (k * (k - 2)) ** 2 * r ** k for k in range(KM + 1))).sqrt() + eps * eps * Mpp
            Tpsi = etan * Rn * Rpp
            poly = -(T2) - eps ** 3 * (LMpp + MLpp) - eps ** 4 * MMpp
            row["eps=" + name] = {"T3/T2": fmt(T3 / T2), "T4/T2": fmt(T4 / T2),
                                  "Tpsi/T2": fmt(Tpsi / T2), "-<R,R''> sign": "neg" if poly < 0 else "nonneg"}
            if n >= 8:
                check(poly < 0, "E diag: polynomial part -<R,R''>(1/2)<0, n=%d eps=%s" % (n, name))
        diag.append(row)
    summary["E_diagnostics_q47_t_half"] = diag
    summary["E_diagnostics_note"] = ("Finite-n diagnostics for h_n WITHOUT moment correction; "
                                     "T3,T4 are exact inner products, Tpsi is the proof's crude bound "
                                     "eta_n ||R|| ||R''|| with eta_n from B_n. Not a proof.")
    # trend checks (diagnostic): ratio <L,L''>/(n^2 S) -> 4*Gamma0
    lastratio = D(diag[-1]["<L,L''>/(n^2 S)"])
    check(abs(lastratio - 4 * G0) / (4 * G0) < D("0.15"), "E diag: <L,L''>/(n^2 S) near 4 Gamma_0 at n=24", v=lastratio)

    # ======================= Remark R: radius decay ========================
    rows = []
    for Qs, expect in [("0.50", D("1.7134e-5")), ("0.53", D("4.8520e-6")),
                       ("0.56", D("2.5717e-7")), ("0.57", D("2.3569e-9"))]:
        d = radius(Qs)["delta"]
        check(abs(d - expect) / expect < D("5e-4"), "R: reproduces Part I table radius at Q=" + Qs, d=d)
    pi = radius("0.5")["pi"]
    r = D(2) / D(7)
    s = D(4) / D(7)

    def c(u):
        return 1 / (2 * pi * (1 - u * u).sqrt())
    cB = c(r).sqrt() * dpow(D(7) / (5 * pi), D("0.25"))
    CJ = c(s) / cB * (D(14) / D(15)) / (D(28) / D(33)) * dpow(pi, D("0.25"))
    Cconst = (D(49) / D(96)) ** 2 * dpow(D(343) / D(144), D("0.25")) / (64 * CJ)
    prev = None
    for k in range(3, 15):
        x = D(10) ** (-k)
        Q = D(4) / D(7) - x
        d = radius(Q)["delta"]
        ratio = d / dpow(x, D("2.25"))
        rows.append({"x": "1e-%d" % k, "delta/x^(9/4)": fmt(ratio), "ratio/C": fmt(ratio / Cconst)})
        if k >= 6:
            check(abs(ratio / Cconst - 1) < D(10) ** (-(k - 3)), "R: delta_Q/x^(9/4) -> C, x=1e-%d" % k)
        # local exponent from consecutive decades
        if prev is not None:
            pe = (d / prev[1]).ln() / (x / prev[0]).ln()
            rows[-1]["local_exponent"] = fmt(pe)
            if k >= 6:
                check(abs(pe - D("2.25")) < D(10) ** (-(k - 4)), "R: local exponent -> 9/4 at x=1e-%d" % k, pe=pe)
            if k == 10:
                negative_control(abs(pe - 2) > D("0.1"), "wrong exponent p=2 rejected", pe=pe)
                negative_control(abs(pe - D("2.5")) > D("0.1"), "wrong exponent p=5/2 rejected", pe=pe)
        prev = (x, d)
    summary["radius_decay"] = {"C_closed": fmt(Cconst), "C_J": fmt(CJ), "rows": rows}

    # ======================= further negative controls =====================
    # wrong envelope: claim s* = b(r)+0.01 at r=0.2 -> s=b+0.005 must already diverge
    r = D("0.2")
    m, _ = maxratio(r, b_of(r) + D("0.005"))
    negative_control(m > 1, "wrong envelope s*=b(r)+0.01 rejected (ratio>1 at b+0.005)", m=m)
    # wrong extremal a: Luo's a with r=2/7 replaced by 1.0 must not give equality tau=nu
    negative_control(abs(tau(ONE, D(2) / 7) - nu(ONE, D(2) / 7)) > D("0.01"),
                     "wrong extremal a=1 at r=2/7 rejected (tau!=nu)")

    all_passed = all(x["passed"] for x in RESULTS)
    os.makedirs(args.output_dir, exist_ok=True)
    out = {"all_passed": all_passed, "n_checks": len(RESULTS),
           "n_failed": sum(1 for x in RESULTS if not x["passed"]),
           "failed": [x for x in RESULTS if not x["passed"]],
           "negative_controls": CONTROLS, "summary": summary,
           "python_optimize_flag": sys.flags.optimize,
           "runtime_s": round(time.time() - t0, 2),
           "checks": RESULTS}
    with open(os.path.join(args.output_dir, "extension_checks.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("checks:", len(RESULTS), "failed:", out["n_failed"],
          "negative controls caught:", sum(1 for c_ in CONTROLS if c_["caught"]), "/", len(CONTROLS),
          "all_passed:", all_passed, "runtime_s:", out["runtime_s"])
    for x in out["failed"]:
        print("FAILED:", x["label"], x["data"])
    return all_passed


if __name__ == "__main__":
    try:
        ok = main()
    except Exception as exc:  # fail closed
        print("EXCEPTION:", repr(exc))
        ok = False
        try:
            od = "./extension_check_output"
            for i, arg in enumerate(sys.argv):
                if arg == "--output-dir" and i + 1 < len(sys.argv):
                    od = sys.argv[i + 1]
                elif arg.startswith("--output-dir="):
                    od = arg.split("=", 1)[1]
            os.makedirs(od, exist_ok=True)
            with open(os.path.join(od, "extension_checks.json"), "w") as fh:
                json.dump({"all_passed": False, "exception": repr(exc), "checks": RESULTS}, fh, indent=1)
        except Exception:
            pass
    sys.exit(0 if ok else 1)

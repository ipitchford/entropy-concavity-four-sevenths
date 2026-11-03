#!/usr/bin/env python3
"""Fail-closed high-precision checks of the Part II covariance identities.

Standard library only (decimal at 80 significant digits, fractions, math).

1. Hermite coefficients of h_n (M4) in the orthonormal basis He_j/sqrt(j!) are
   computed from the double generating function
       sqrt2 (1+a)^(-1/2) exp(b z^2/2 + sqrt(2a) z w/(1+a) - a w^2/(2(1+a))),
   and cross-checked against direct Gaussian-moment integration of
   h_n(x) He_j(x) for small j, and against the linear generating function
   (M11) at r = 0.3, 0.6, 1.
2. M_t is built directly from (P14) with rho = q, z = t (degree-m coefficient
   sum_j beta_j beta_(m-j) binom(m,j)^(1/2) (qt)^(j/2) (q(1-t))^((m-j)/2)),
   and ||M_t||^2 is compared with (M20): 4 (n!)^2 C_n(t-1/2) / Delta, where
   C_n is obtained by independent multigraph enumeration with the matrix (M21).
   q in {4/7, 0.6, 1}, n in {2,4,6}, t in {0.45, 0.5, 0.6}.
3. Lemma II.5: <M_t, M_s> is compared with 4 (n!)^2 C_n(G(t,s)) / Delta for the
   stated two-weight block matrix G, and every entry of G is checked to be
   strictly positive on a grid near (1/2, 1/2). G(t,t) is checked to equal
   B(t-1/2).

Every comparison is relative with tolerance REL_TOL = 1e-40. Every result goes
through Checker.check; any failure or exception gives exit status 1. No assert.
This is a port of the methodology reviewer's mpmath scripts m20.py/lemma4.py.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
import traceback
from decimal import Decimal as D, getcontext
from pathlib import Path

getcontext().prec = 80
REL_TOL = D("1e-40")
TAIL_TOL = D("1e-60")
JMAX = 360


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
        return dict(number_of_checks=len(self.results), number_failed=len(failed),
                    failed_labels=failed, all_checks_passed=(not failed) and bool(self.results))


def fmt(x, digits=20):
    return format(x, f".{digits}g")


def rel_err(x, y):
    scale = max(abs(x), abs(y))
    return abs(x - y) / scale if scale else abs(x - y)


S5 = D(5).sqrt()
A = D(7) / (3 * S5)
B_ = (A - 1) / (A + 1)
FACT = [D(math.factorial(k)) for k in range(2 * JMAX + 10)]
SQRT_FACT = [f.sqrt() for f in FACT]


def params(q):
    lp = (A * (1 + q) - 1) / (A * (1 + q) + 1)
    lm = (A * (1 - q) - 1) / (A * (1 - q) + 1)
    al = (lp + lm) / 2 - B_
    bt = (lp - lm) / 2
    delta = (1 + A) * ((1 + A) ** 2 - A * A * q * q).sqrt()   # (M20)
    return dict(d=-B_ - al / 2, e=-al / 2, c=bt / 2, alpha=al, beta=bt, Delta=delta, kappa=-lm)


def beta_coeffs(n, jmax):
    """Coefficients of h_n on He_j/sqrt(j!) from the double generating function."""
    a_ = B_ / 2
    c_ = (2 * A).sqrt() / (1 + A)
    d_ = -A / (2 * (1 + A))
    pre = D(2).sqrt() / (1 + A).sqrt()
    out = []
    for j in range(jmax + 1):
        tot = D(0)
        for m in range(0, min(n, j) + 1):
            if (n - m) % 2 or (j - m) % 2:
                continue
            tot += c_**m / FACT[m] * a_ ** ((n - m) // 2) / FACT[(n - m) // 2] * d_ ** ((j - m) // 2) / FACT[(j - m) // 2]
        out.append(pre * SQRT_FACT[n] * SQRT_FACT[j] * tot)
    return out


def hermite_poly(n):
    """Integer coefficients of the probabilists' Hermite polynomial He_n."""
    p0, p1 = [1], [0, 1]
    if n == 0:
        return p0
    for k in range(1, n):
        nxt = [0] + p1
        for i, coef in enumerate(p0):
            nxt[i] -= k * coef
        p0, p1 = p1, nxt
    return p1


def beta_direct(n, j):
    """<h_n, He_j/sqrt(j!)>_phi by Gaussian-moment integration of the polynomial product."""
    hn = hermite_poly(n)
    hj = hermite_poly(j)
    s2a = (2 * A).sqrt()
    total = D(0)
    for k, ck in enumerate(hn):
        if ck == 0:
            continue
        for l, cl in enumerate(hj):
            if cl == 0 or (k + l) % 2:
                continue
            m = k + l
            dfact = math.prod(range(m - 1, 0, -2)) if m > 0 else 1
            moment = D(dfact) / (1 + A) ** (D(m + 1) / 2)       # int x^m e^{-a x^2/2} phi
            total += D(ck) * s2a**k * D(cl) * moment
    return D(2).sqrt() * total / (SQRT_FACT[n] * SQRT_FACT[j])


def m11_rhs(n, r):
    tau = (1 - A + A * r) / (1 + A - A * r)
    nu = (A - 1 + A * r) / (1 + A + A * r)
    pre = 2 / ((1 + A) ** 2 - A * A * r * r).sqrt()
    ck = [D(math.comb(2 * k, k)) / D(4) ** k for k in range(n + 1)]
    return pre * sum(ck[k] * ck[n - k] * tau**k * nu ** (n - k) for k in range(n + 1))


def m_coeffs(be, q, t):
    """Hermite coefficients of M_t from (P14) with rho=q, z=t."""
    x = (q * t).sqrt()
    y = (q * (1 - t)).sqrt()
    jmax = len(be) - 1
    xp = [D(1)]
    yp = [D(1)]
    for _ in range(jmax):
        xp.append(xp[-1] * x)
        yp.append(yp[-1] * y)
    out = []
    for m in range(jmax + 1):
        tot = D(0)
        for j in range(m + 1):
            if be[j] == 0 or be[m - j] == 0:
                continue
            tot += be[j] * be[m - j] * D(math.comb(m, j)).sqrt() * xp[j] * yp[m - j]
        out.append(tot)
    return out


def b_matrix(q, u):
    p = params(q)
    d, e, c = p["d"], p["e"], p["c"]
    w = (1 - 4 * u * u).sqrt()
    return [[d + 2 * e * u, e * w, c * (1 + 2 * u), c * w],
            [e * w, d - 2 * e * u, c * w, c * (1 - 2 * u)],
            [c * (1 + 2 * u), c * w, d + 2 * e * u, e * w],
            [c * w, c * (1 - 2 * u), e * w, d - 2 * e * u]]          # (M21)


def g_matrix(q, t, s):
    """Lemma II.5 phase-changed two-weight matrix."""
    p = params(q)
    al, bt = p["alpha"], p["beta"]
    vt = [t.sqrt(), (1 - t).sqrt()]
    vs = [s.sqrt(), (1 - s).sqrt()]
    m = [[D(0)] * 4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            m[i][j] = (-B_ if i == j else 0) - al * vt[i] * vt[j]
            m[2 + i][2 + j] = (-B_ if i == j else 0) - al * vs[i] * vs[j]
            m[i][2 + j] = bt * vt[i] * vs[j]
            m[2 + j][i] = m[i][2 + j]
    return m


def coefficient_cn(bm, n):
    """[z1^n z2^n z3^n z4^n] exp(z^T B z / 2) by multigraph (edge-count) enumeration."""
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    total = [D(0)]

    def rec(k, deg, wt):
        if k == len(pairs):
            if any(x % 2 for x in deg):
                return
            for i in range(4):
                kk = deg[i] // 2
                wt *= (bm[i][i] / 2) ** kk / FACT[kk]
            total[0] += wt
            return
        i, j = pairs[k]
        for mm in range(min(deg[i], deg[j]) + 1):
            d2 = list(deg)
            d2[i] -= mm
            d2[j] -= mm
            rec(k + 1, d2, wt * bm[i][j] ** mm / FACT[mm])

    rec(0, [n] * 4, D(1))
    return total[0]


def run(ck, out):
    qs = [("4/7", D(4) / 7), ("0.6", D("0.6")), ("1", D(1))]
    ns = [2, 4, 6]
    ts = [D("0.45"), D("0.5"), D("0.6")]
    betas = {n: beta_coeffs(n, JMAX) for n in ns}

    # 1. Hermite coefficients: generating function vs direct integration; (M11).
    for n in ns:
        worst = D(0)
        for j in range(0, 15):
            bd = beta_direct(n, j)
            bg = betas[n][j]
            err = abs(bd - bg) / max(abs(bd), abs(bg), D("1e-70"))
            worst = max(worst, err)
        ck.check(worst <= REL_TOL, f"h_n n={n}: generating-function coefficients equal direct integration, j<=14",
                 fmt(worst, 3), ref="Part II (M4)")
        out["hermite_coefficients"].append(dict(n=n, max_rel_error_j_le_14=fmt(worst, 3)))
        for r in [D("0.3"), D("0.6"), D(1)]:
            lhs = sum(bj * bj * r**j for j, bj in enumerate(betas[n]))
            rhs = m11_rhs(n, r)
            e = rel_err(lhs, rhs)
            ck.check(e <= REL_TOL, f"(M11) n={n} r={r}: sum_j beta_j^2 r^j equals generating-function coefficient",
                     fmt(e, 3), ref="Part II (M11)")
            out["M11"].append(dict(n=n, r=str(r), hermite=fmt(lhs), formula=fmt(rhs), rel_error=fmt(e, 3)))

    # 2. (M20)-(M21).
    cache = {}
    for qlab, q in qs:
        p = params(q)
        for n in ns:
            for t in ts:
                mc = m_coeffs(betas[n], q, t)
                cache[(qlab, n, t)] = mc
                norm2 = sum(x * x for x in mc)
                tail = sum(x * x for x in mc[-40:])
                ck.check(tail / norm2 <= TAIL_TOL, f"q={qlab} n={n} t={t}: truncation tail of ||M_t||^2 negligible",
                         fmt(tail / norm2, 3), ref="numerical truncation")
                bm = b_matrix(q, t - D("0.5"))
                ck.check(min(min(row) for row in bm) > 0, f"q={qlab} t={t}: entries of B(t-1/2) positive",
                         fmt(min(min(row) for row in bm), 6), ref="Part II (M21), Lemma II.3 hypothesis")
                rhs = 4 * FACT[n] ** 2 / p["Delta"] * coefficient_cn(bm, n)
                e = rel_err(norm2, rhs)
                ck.check(e <= REL_TOL, f"(M20)-(M21) q={qlab} n={n} t={t}: ||M_t||^2 = 4(n!)^2 C_n(t-1/2)/Delta",
                         fmt(e, 3), ref="Part II (M20),(M21)")
                out["M20_M21"].append(dict(q=qlab, n=n, t=str(t), hermite_norm_squared=fmt(norm2, 25),
                                           formula=fmt(rhs, 25), rel_error=fmt(e, 3), tail_ratio=fmt(tail / norm2, 3)))

    # 3. Lemma II.5 two-weight block matrix.
    pairs = [(D("0.5"), D("0.6")), (D("0.45"), D("0.6")), (D("0.45"), D("0.5")), (D("0.6"), D("0.45"))]
    for qlab, q in qs:
        p = params(q)
        for t in ts:
            g = g_matrix(q, t, t)
            bm = b_matrix(q, t - D("0.5"))
            e = max(abs(g[i][j] - bm[i][j]) for i in range(4) for j in range(4))
            ck.check(e <= D("1e-70"), f"Lemma II.5 q={qlab} t={t}: G(t,t) equals B(t-1/2) of (M21)", fmt(e, 3),
                     ref="Part II Lemma II.5, (M21)")
        grid = [D("0.42") + D("0.02") * k for k in range(9)]   # 0.42 ... 0.58
        mins = min(min(min(row) for row in g_matrix(q, t, s)) for t in grid for s in grid)
        ck.check(mins > 0, f"Lemma II.5 q={qlab}: every entry of G(t,s) > 0 for t,s in [0.42,0.58]", fmt(mins, 6),
                 ref="Part II Lemma II.5")
        out["lemma4_positivity"].append(dict(q=qlab, grid="t,s in {0.42,0.44,...,0.58}", min_entry=fmt(mins, 10)))
        for n in ns:
            for t, s in pairs:
                mt = cache[(qlab, n, t)]
                ms = cache[(qlab, n, s)]
                lhs = sum(x * y for x, y in zip(mt, ms))
                g = g_matrix(q, t, s)
                rhs = 4 * FACT[n] ** 2 / p["Delta"] * coefficient_cn(g, n)
                e = rel_err(lhs, rhs)
                ck.check(e <= REL_TOL, f"Lemma II.5 q={qlab} n={n} t={t} s={s}: <M_t,M_s> = 4(n!)^2 C_n(G(t,s))/Delta",
                         fmt(e, 3), ref="Part II Lemma II.5")
                ck.check(min(min(row) for row in g) > 0, f"Lemma II.5 q={qlab} t={t} s={s}: G entries positive",
                         "", ref="Part II Lemma II.5")
                out["lemma4"].append(dict(q=qlab, n=n, t=str(t), s=str(s), hermite_inner_product=fmt(lhs, 25),
                                          formula=fmt(rhs, 25), rel_error=fmt(e, 3),
                                          min_entry=fmt(min(min(row) for row in g), 8)))
    ck.check(len(out["M20_M21"]) == 27 and len(out["lemma4"]) == 36, "expected number of identity cases evaluated",
             (len(out["M20_M21"]), len(out["lemma4"])))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("identity_output"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    ck = Checker()
    out = {"scope": "High-precision (80-digit Decimal) checks of (M11), (M20)-(M21) and the Lemma II.5 two-weight "
                    "matrix; finite cases, not a proof.",
           "decimal_precision": getcontext().prec, "relative_tolerance": str(REL_TOL), "hermite_truncation": JMAX,
           "python_optimize_flag": sys.flags.optimize,
           "hermite_coefficients": [], "M11": [], "M20_M21": [], "lemma4": [], "lemma4_positivity": []}
    try:
        run(ck, out)
    except Exception as exc:
        ck.check(False, "uncaught exception", repr(exc) + "\n" + traceback.format_exc())
    out.update(ck.summary())
    out["checks"] = ck.results
    out["elapsed_seconds"] = time.perf_counter() - start
    (args.output_dir / "identity_checks.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    lines = [f"{'PASS' if r['passed'] else 'FAIL'} {r['label']}  [{r['detail'][:60]}]" for r in ck.results]
    status = "All checks passed" if out["all_checks_passed"] else f"{out['number_failed']} CHECK(S) FAILED"
    lines.append(f"{status}: {out['number_of_checks']} checks in {out['elapsed_seconds']:.3f} s.")
    (args.output_dir / "identity_checks.log").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(lines[-1], flush=True)
    sys.exit(0 if out["all_checks_passed"] else 1)


if __name__ == "__main__":
    main()

---
title: "The critical signal fraction for local entropy concavity"
subtitle: "A Gaussian neighbourhood theorem, midpoint curvature, and the remaining boundary point"
date: "2 November 2026"
lang: en-GB
fontsize: 10pt
geometry:
  - margin=23mm
colorlinks: true
linkcolor: MidnightBlue
urlcolor: MidnightBlue
---

# Overview and main results

## Abstract

We determine the critical retained-signal fraction for local mixing-weight entropy concavity near the standard Gaussian. For every \(Q<4/7\), we construct a positive relative-density radius such that every symmetric variance-one input in that neighbourhood has a strictly concave entropy curve at every signal fraction \(0<q\le Q\), unless the input is Gaussian. The neighbourhood theorem does not assume log-concavity. Together with Luo's existing counterexamples above \(4/7\), this identifies the supremal local-concavity signal fraction as exactly \(4/7\). Membership of the boundary value itself remains unresolved. We also strengthen the counterexample side: throughout \(q>4/7\), smooth strongly log-concave perturbations can make equal weighting a strict local entropy minimum, with any prescribed finite Gaussian moment matches. A differentiated saddle calculation supplies the curvature constant and a Gaussian local profile on the weight scale \(n^{-1/2}\). These statements answer Luo's near-Gaussian question throughout the strict lower range and his midpoint-curvature question in full. The positive proof upgrades a Gaussian-kernel interpolation estimate to bounded functions and then controls nonlinear derivatives after additional smoothing.

## The question

Let \(X,Y\) be independent copies of a mean-zero, variance-one random variable. Let \(G\) be an independent standard Gaussian and write
\[
F_{q,f}(t)=
h\!\left(\sqrt{qt}\,X+\sqrt{q(1-t)}\,Y+\sqrt{1-q}\,G\right),
\qquad 0\le t\le1,
\]
where \(f\) is the common density and \(h\) is differential entropy with natural logarithms. This equals the entropy obtained by smoothing each input separately with retained-signal fraction \(q\), then combining the independent outputs with weights \(\sqrt t,\sqrt{1-t}\). Both formulations preserve variance one.

The original concavity conjecture is due to Ball, Nayar and Tkocz. Luo's 2026 counterexamples show that symmetry and arbitrarily close approximation to the Gaussian do not guarantee concavity when \(q>4/7\). Schubr's subsequent neighbourhood theorem proves a positive result below \(q=1/2\). The strict interval between those two values is the main gap addressed here.

Let \(\varphi\) be the standard Gaussian density and define
\[
T_\rho g(x)=\mathbb E\,g(\rho x+\sqrt{1-\rho^2}\,G).
\]
On the orthonormal probabilists' Hermite basis \(e_j=H_j/\sqrt{j!}\), this operator satisfies \(T_\rho e_j=\rho^j e_j\); the number operator is \(Ne_j=je_j\). Define
\[
\mathcal C=
\left\{q\in(0,1]:
\begin{array}{l}
\text{there exists }\delta>0\text{ such that every even density}\\
f=\varphi(1+g),\quad \int g\varphi=\int x^2g\varphi=0,\quad
\|g\|_\infty\le\delta,\\
\text{has a concave curve }F_{q,f}\text{ on }[0,1]
\end{array}\right\}.
\]

**Theorem A (local concavity).** For every \(Q\in(0,4/7)\), there is a computable \(\delta_Q>0\) such that every admissible even real bounded perturbation \(g\) with \(\|g\|_\infty\le\delta_Q\) satisfies
\[
F_{q,\varphi(1+g)}''(t)
\le-\frac1{16}
\left\|N\left(T_{\sqrt{qt}}g+
                T_{\sqrt{q(1-t)}}g\right)\right\|_{L^2(\varphi)}^2
\tag{A}
\]
for all \(0<q\le Q\) and \(t\in[0,1]\), with continuous one-sided derivatives at the weight endpoints. Here \(T_\rho\) is the Ornstein--Uhlenbeck operator and \(N\) multiplies Hermite degree \(j\) by \(j\). For nonzero \(g\), the right-hand side is strictly negative.

**Corollary (exact critical value).**
\[
(0,4/7)\subseteq\mathcal C\subseteq(0,4/7],
\qquad
\boxed{\sup\mathcal C=4/7.}
\tag{B}
\]
Indeed Theorem A gives the first inclusion. Luo's counterexamples give the second. The set is downward closed by Gaussian-semigroup contraction. Thus it is either \((0,4/7)\) or \((0,4/7]\); this article does not distinguish them.

The same critical supremum holds when inputs are restricted to smooth strongly log-concave densities with any prescribed finite number of Gaussian moment matches, and when neighbourhoods use any fixed relative \(C_b^m\) norm. The positive result holds already in the larger bounded-relative-density class, while the counterexample sequence below converges in every fixed relative \(C_b^m\) norm.

**Theorem B (local minima above the threshold).** For every integer \(K\ge1\), Luo's moment-corrected Hermite construction gives a single sequence \(f_n\) of smooth, positive, even, strongly log-concave densities matching Gaussian moments through degree \(2K\), with
\[
\|f_n/\varphi-1\|_{C_b^m}\longrightarrow0
\]
for every fixed \(m\). For every fixed \(q\in(4/7,1]\), eventually
\[
F_{q,f_n}''(1/2)>0.
\tag{C}
\]
The smoothed densities have the same approximation and strong-log-concavity properties. With
\(D_{q,n}=h(\varphi)-F_{q,f_n}(1/2)\),
\[
\frac{F_{q,f_n}(1/2+s/\sqrt n)-F_{q,f_n}(1/2)}
     {D_{q,n}}
\longrightarrow
1-\exp\!\left[-\frac{84\sqrt5\,q}{49q-4}s^2\right]
\tag{D}
\]
in local \(C^2\). The proof gives the leading constant in \(F''(1/2)\), not only its sign.

## What is resolved and what remains

| Signal variance fraction | Local statement near the Gaussian | Status established here |
|---|---|---|
| \(0<q<4/7\) | A positive relative \(L^\infty\) neighbourhood has all-weight concavity | Theorem A; extends the previously inspected \(q<1/2\) neighbourhood result |
| \(q=4/7\) | Neighbourhood concavity or arbitrarily close counterexamples | Unresolved |
| \(4/7<q\le1\) | Counterexamples can approach the Gaussian in every fixed relative \(C_b^m\) norm | Nonconcavity already known; Theorem B adds positive midpoint curvature throughout this range and its local profile |

Luo's Problem 10.2 asks whether near-Gaussian counterexample sequences exist in the lower signal range. Theorem A gives a negative answer for each \(q<4/7\), including all extra regularity and moment conditions in that question. The boundary instance remains open. Theorem B gives an affirmative answer to Problem 10.3.

These are local theorems. They do not prove concavity for every symmetric log-concave density at a given positive noise level. They also do not assert that asymmetric inputs behave in the same way.

In the equivalent unnormalised observation \(X+\sigma G\), the signal fraction is
\[
q=\frac{\operatorname{Var}(X)}
        {\operatorname{Var}(X)+\sigma^2}.
\]
Thus \(q=4/7\) corresponds to
\(\sigma^2=(3/4)\operatorname{Var}(X)\). The equality case remains exactly the same unresolved boundary after rescaling.

## Proof architecture and attribution

Part I proves Theorem A. The new step is a weighted version of Luo's interpolation estimate: a weak Gaussian weight makes the input norm finite for bounded perturbations while a strict scalar kernel inequality survives. Additional smoothing then controls weight derivatives of the quadratic convolution term. Schubr's Hermite-energy argument absorbs the nonlinear logarithmic terms.

Part II proves Theorem B. The Hermite construction, its moment correction and its four-variable generating coefficient are Luo's. The new calculation evaluates this coefficient at its moving saddle, controls two derivatives directly and transfers the resulting sign to entropy. A two-weight kernel bound is essential: differentiating the squared norm alone would not control the entropy error.

Both proof parts include the arguments on which they depend rather than treating the cited preprints as unexplained black boxes. Equation labels use P and M to distinguish the two parts; section numbering restarts within each part. The accompanying reproducibility files test exact finite identities and numerical asymptotics. They are not substitutes for the analytic proofs, and no proof-assistant verification or external peer review is claimed.

\newpage


# Part I. Local concavity below four sevenths

## 1. A Gaussian-kernel comparison

Let
\[
K_j(x,y)=c_j\exp[-\alpha_j(x^2+y^2)+\beta_jxy],
\quad c_j>0,\quad 2\alpha_j>\beta_j>0,\qquad j=A,B.
\]
If
\[
\beta_A+2|\alpha_A-\alpha_B|\le\beta_B,                   \tag{P2}
\]
then, in the order of quadratic forms on \(L^2(dx)\),
\[
K_A\preceq
\frac{c_A}{c_B}\frac{\beta_B}{\beta_A}K_B.                \tag{P3}
\]

This is the elementary Gaussian-kernel lemma used in Luo's Appendix F. We include the kernel comparison for completeness.

In Fock space with norm
\[
\|f\|_{\mathcal F}^2=\pi^{-1}\int_{\mathbb C}|f(z)|^2e^{-|z|^2}\,dA(z),
\]
the reproducing kernel is \(k_w(z)=e^{z\overline w}\). Put
\[
\Delta=\alpha_A-\alpha_B,\quad \rho^2=\beta_A/\beta_B,\quad
(Cf)(z)=e^{-\Delta z^2/\beta_B}f(\rho z).
\]
After the change of variables \(w=\rho z\), condition (P2) implies that each real-direction Gaussian decay coefficient is at least one. Therefore
\[
\|Cf\|_{\mathcal F}^2\le\rho^{-2}\|f\|_{\mathcal F}^2.
\]
For real \(x\), set \(\xi_x=e^{-\alpha_Bx^2}k_{\sqrt{\beta_B}x}\). Then
\[
\langle\xi_y,\xi_x\rangle=K_B(x,y)/c_B,\qquad
\langle C^*\xi_y,C^*\xi_x\rangle=K_A(x,y)/c_A.
\]
The norm bound proves (P3) first for finite Gram forms. Riemann-sum approximation and density extend it to \(L^2(dx)\); the Gaussian kernels are Hilbert--Schmidt because \(2\alpha_j>\beta_j\).

## 2. Weighted interpolation for bounded functions

For \(0<u<1\), define
\[
c_u=\frac1{2\pi\sqrt{1-u^2}},\qquad
\alpha_u=\frac1{2(1-u^2)},\qquad
\beta_u=\frac{u}{1-u^2}.
\]
The positive operator \(A_u\) on \(L^2(dx)\), with kernel
\[
A_u(x,y)=c_u e^{-\alpha_u(x^2+y^2)+\beta_uxy},
\]
satisfies
\[
\langle g,A_ug\rangle_{dx}=\|T_{\sqrt u}g\|_2^2.          \tag{P4}
\]
Indeed this kernel is the joint standard Gaussian density with correlation \(u\).

### Lemma 2.1

Suppose
\[
0<r<s<\sqrt{2r(1+r)}-r.                                \tag{P5}
\]
There are computable \(\varepsilon>0\) and \(J_{r,s}<\infty\) such that every real bounded measurable \(g\) satisfies
\[
\boxed{\quad
\|T_{\sqrt s}g\|_2^2
\le J_{r,s}\|g\|_\infty\,\|T_{\sqrt r}g\|_2.
\quad}                                                   \tag{P6}
\]

### Proof and closed-form constants

Put
\[
P_0=\frac1{1-r},\qquad \beta=\beta_r=\frac r{1-r^2},
\]
\[
B_1=\frac1{1-s}-P_0,\qquad
B_2=P_0-\frac1{1+s}.
\]
Both \(B_1,B_2\) are positive because \(r<s\). Define
\[
E_i=P_0-\frac{B_i^2}{2\beta}\quad(i=1,2),\qquad
\boxed{\quad
\varepsilon=\frac12\min\left\{
\frac1{1+r},\,\frac1{1+s},\,E_1,\,E_2
\right\}.
\quad}                                                   \tag{P7}
\]
We first verify that this choice is strictly positive. If
\(c=\sqrt{2r/(1+r)}\), then
\(\sqrt{2\beta P_0}=P_0c\).
Elementary rearrangement gives
\[
B_1<\sqrt{2\beta P_0}
\quad\Longleftrightarrow\quad
s<\frac{r+c}{1+c}=\sqrt{2r(1+r)}-r,
\]
and
\[
B_2<\sqrt{2\beta P_0}
\quad\Longleftrightarrow\quad
s<\frac{c-r}{1-c}=\sqrt{2r(1+r)}+r.
\]
The second upper bound exceeds the first by \(2r\). Assumption (P5) therefore implies \(E_1,E_2>0\), proving positivity of (P7).

Set \(h(x)=e^{-\varepsilon x^2/2}g(x)\). Then
\[
\|h\|_{L^2(dx)}
\le(\pi/\varepsilon)^{1/4}\|g\|_\infty.
\]
Let \(A_u^\varepsilon\) be the kernel \(A_u\) with
\(\alpha_u^\varepsilon=\alpha_u-\varepsilon/2\).
Conjugation is exact:
\[
\langle h,A_u^\varepsilon h\rangle_{dx}
=\|T_{\sqrt u}g\|_2^2.                                  \tag{P8}
\]
The bounds \(\varepsilon<1/(1+r),1/(1+s)\) ensure
\(2\alpha_u^\varepsilon>\beta_u>0\) for \(u=r,s\).

The positive square root of \(A_r^\varepsilon\) has the particularly simple kernel
\[
B_\varepsilon(x,y)=
c_B e^{-P(x^2+y^2)/2+\beta_Bxy},\qquad
P=P_0-\varepsilon,
\]
\[
\beta_B=\sqrt{2\beta P},\qquad
c_B=\sqrt{c_r}\,(P/\pi)^{1/4}.                           \tag{P9}
\]
To verify this directly, Gaussian integration gives
\[
\int_{\mathbb R}B_\varepsilon(x,z)B_\varepsilon(z,y)\,dz
=c_B^2\sqrt{\pi/P}\,
\exp\left[-\left(\frac P2-\frac{\beta_B^2}{4P}\right)
(x^2+y^2)+\frac{\beta_B^2}{2P}xy\right].
\]
The prefactor is \(c_r\), the cross coefficient is \(\beta_r\), and the diagonal coefficient is
\((P-\beta_r)/2=\alpha_r-\varepsilon/2\). Thus \(B_\varepsilon^2=A_r^\varepsilon\).
Moreover \(P>\beta_B>0\): squaring reduces this to
\(P>2\beta_r\), equivalent to \(\varepsilon<1/(1+r)\).
The kernel is bounded and Hilbert--Schmidt. It is a positive operator because expansion of \(e^{\beta_Bxy}\) gives a sum of positive rank-one quadratic forms; approximation by compactly supported functions justifies the expansion and passage to \(L^2(dx)\). Hence (P9) is the positive square root, rather than merely an algebraic square root.

The kernel-comparison condition (P2) between \(A_s^\varepsilon\) and \(B_\varepsilon\) is equivalent to
\[
2\alpha_s^\varepsilon+\beta_s<P+\beta_B,\qquad
P-\beta_B<2\alpha_s^\varepsilon-\beta_s.
\]
Substitute \(2\alpha_s\pm\beta_s=1/(1\mp s)\). These inequalities become
\[
B_1<\beta_B,\qquad B_2<\beta_B.
\]
By (P7), \(\varepsilon<E_i\), so
\(\beta_B^2=2\beta(P_0-\varepsilon)>B_i^2\). Both inequalities are strict.

Apply (P3), then Cauchy--Schwarz:
\[
\|T_{\sqrt s}g\|_2^2
\le C_\varepsilon\langle h,B_\varepsilon h\rangle
\le C_\varepsilon\|h\|_{L^2(dx)}
                  \|B_\varepsilon h\|_{L^2(dx)},
\]
\[
C_\varepsilon=\frac{c_s}{c_B}\frac{\beta_B}{\beta_s}.
\]
By (P8), \(\|B_\varepsilon h\|_{L^2(dx)}=\|T_{\sqrt r}g\|_2\). Thus (P6) holds with
\[
\boxed{\quad
J_{r,s}=\frac{c_s}{c_B}\frac{\beta_B}{\beta_s}
(\pi/\varepsilon)^{1/4}.
\quad}                                                   \tag{P10}
\]
This proves the weighted interpolation inequality with a closed formula for every constant.

## 3. Bilinear derivatives after smoothing

For \(0<\rho<1\), let
\[
\mathcal B_\rho(t;h,h)(x)=
\mathbb E[h(G_1)h(G_2)\mid
\sqrt{\rho t}\,G_1+\sqrt{\rho(1-t)}\,G_2+
\sqrt{1-\rho}\,G_3=x].
\]
Suppose \(h\in L^2(\varphi)\) is even, with Hermite coefficients zero below degree four. Set
\[
\eta_\rho=\frac{\rho^{-1}-1}{4},\qquad
R_\rho=\frac{1+\rho}{2},\qquad
K_\rho=\frac{\sqrt{R_\rho}}{(1-\sqrt{R_\rho})^2}.
\]
Then, uniformly in \(0\le t\le1\),
\[
\|\mathcal B_\rho(t;h,h)\|_2\le\|h\|_2^2,               \tag{P11}
\]
\[
\|N\mathcal B_\rho(t;h,h)\|_2\le K_\rho\|h\|_2^2,        \tag{P12}
\]
\[
\|\partial_t^2\mathcal B_\rho(t;h,h)\|_2
\le2\eta_\rho^{-2}\|h\|_2^2.                            \tag{P13}
\]
The bilinear map is \(C^2\) as an \(L^2\)-valued function through both endpoints.

### Proof

If \(h=\sum b_ne_n\), its output coefficient of degree \(n\) at complex weight \(z\) is
\[
\rho^{n/2}\sum_{j+k=n}b_jb_k\binom nj^{1/2}
z^{j/2}(1-z)^{k/2}.                                    \tag{P14}
\]
Only even \(j,k\ge4\) occur, so the terms are polynomials in \(z\). For \(|z-t|\le\eta_\rho\),
\[
P=\rho(|z|+|1-z|)
\le \rho(1+2\eta_\rho)=R_\rho<1.
\]
Put \(w=|z|/(|z|+|1-z|)\), and let \(h_{\rm abs}\) have coefficients \(|b_n|\). Taking absolute coefficients in (P14) reduces the norm to the real conditional-expectation expression
\(\mathcal B_P(w;h_{\rm abs},h_{\rm abs})\).
Contraction of conditional expectation and independence give
\[
\|\mathcal B_P(w;h_{\rm abs},h_{\rm abs})\|_2
\le\|h_{\rm abs}\|_2^2=\|h\|_2^2.
\]
This proves the bound on the whole disk. To include \(N\), split off \(P^{n/2}\) from the coefficient of the unit-signal conditional expectation. Its norm is bounded by \(\|h\|_2^2\), and
\[
\sup_{n\ge0}nP^{n/2}
\le\sum_{n\ge1}nR_\rho^{n/2}=K_\rho.
\]
There remains a strictly positive smoothing margin beyond the stated disk. Absolute-coefficient majorisation and the exponentially decaying output multipliers give normal convergence, and likewise for every fixed number-operator derivative. Hence the map is Hilbert-valued holomorphic there. Cauchy's derivative estimate proves (P13).

Finite Hermite expansions are dense in \(L^2\), and the bilinear contraction estimate, applied to differences, proves convergence to the actual conditional expectation. The same bounds prove \(C^2\) convergence in the weight, so the assertions include both real endpoints.

## 4. Apply the two lemmas at a fixed signal

Fix \(q\in(0,4/7)\). Put
\[
r=q/2,\qquad b(r)=\sqrt{2r(1+r)}-r,\qquad
s=\frac{q+b(r)}2.
\]
The inequality \(b(r)>q\) is equivalent to \(q<4/7\). Therefore
\[
q<s<b(r).
\]
Let \(J=J_{r,s}\) be (P10), \(\rho=q/s\), and write
\(\eta=\eta_\rho,\ R_\rho,\ K=K_\rho\).

The exact output density has relative perturbation \(L_g(q,t)+M_g(q,t)\), where \(M_g(q,t)=\mathcal B_q(t;g,g)\). This follows by expanding the product of the two input density ratios in the conditional-expectation representation. For the original perturbation \(g\), set \(h=T_{\sqrt s}g\). The Gaussian/Hermite composition identity gives exactly
\[
M_g(q,t)=\mathcal B_\rho(t;h,h).                         \tag{P15}
\]
Lemma 2.1 gives
\[
\|h\|_2^2\le J\delta\|T_{\sqrt{q/2}}g\|_2
\le\frac{J\delta}{2}\|L_g(q,t)\|_2
\le\frac{J\delta}{8}D,                                  \tag{P16}
\]
where \(\delta=\|g\|_\infty\) and \(D=\|NL_g(q,t)\|_2\). The middle inequality follows from
\(t^{n/2}+(1-t)^{n/2}\ge2^{1-n/2}\) in each Hermite degree.

Equations (P11)--(P13) therefore imply
\[
\|M\|_2\le\frac{J\delta}{8}D,\qquad
\|NM\|_2\le\frac{KJ\delta}{8}D,\qquad
\|M_{tt}\|_2\le\frac{J\delta}{4\eta^2}D.                 \tag{P17}
\]

The order of operations matters: absolute Hermite coefficients are taken only after smoothing the actual \(g\). Their \(L^2\) norm equals \(\|T_{\sqrt s}g\|_2\), which is already controlled by weighted interpolation. We do not claim that the absolute-coefficient transform of a bounded \(g\) is bounded.

## 5. Nonlinear curvature and an explicit radius

We reproduce the credited entropy-energy steps to make the argument self-contained. The output ratio is \(1+\mathcal R\), where \(\mathcal R=L+M\), and
\[
\|\mathcal R\|_\infty\le2\delta+\delta^2\le3\delta.
\]
The linear energy obeys
\[
\langle L_{tt},L\rangle\ge D^2/8,\qquad
\|N^{-1}L_{tt}\|_2\le D,\qquad \|L\|_2\le D/4.           \tag{P18}
\]
These follow from the degree-\(n=2k\) multiplier
\(q^k[t^k+(1-t)^k]\) and
\[
n^2\ell_n/8\le\ell_n''\le n^2\ell_n.
\]

For a smooth bounded real \(v\), Gaussian integration by parts and
\(\|v_{xx}\|_2\le\|Nv\|_2\) give
\[
\|(v_x)^2\|_2\le3\|v\|_\infty\|Nv\|_2.                  \tag{P19}
\]
Indeed, if \(A=\int v_x^4\varphi\), then
\[
\int Nv\,v(v_x)^2\varphi
=A+2\int v(v_x)^2v_{xx}\varphi,
\]
and Cauchy--Schwarz gives
\(A\le3\|v\|_\infty A^{1/2}\|Nv\|_2\).

Put \(\psi(u)=\log(1+u)-u\). The number-operator chain rule and (P19) imply, for \(\delta\le1/12\),
\[
\|N\psi(\mathcal R)\|_2
\le30\delta\|N\mathcal R\|_2
\le30\delta(1+KJ\delta/8)D.                             \tag{P20}
\]

Define
\[
A_q=\frac{KJ}{8}+60+\frac{J}{8\eta^2},
\]
and choose
\[
\boxed{\quad
\delta_q=\min\left\{\frac1{12},\frac1J,\frac8{KJ},
\frac1{16A_q}\right\}>0.
\quad}                                                   \tag{P21}
\]
For \(\delta\le\delta_q\), (P20) is at most \(60\delta D\). Also
\[
\frac{\|\mathcal R\|_2}{1-3\delta}
\le D\,\frac{1/4+J\delta/8}{1-3\delta}
\le D/2.
\]
The entropy formula is
\[
F_{tt}=-\langle L_{tt}+M_{tt},\log(1+\mathcal R)\rangle
-\int\frac{\mathcal R_t^2}{1+\mathcal R}\varphi.
\]
Drop the last nonpositive term and expand
\(\log(1+\mathcal R)=L+M+\psi(\mathcal R)\).
The error terms satisfy
\[
|\langle L_{tt},M\rangle|\le(KJ/8)\delta D^2,\qquad
|\langle L_{tt},\psi(\mathcal R)\rangle|\le60\delta D^2,
\]
\[
|\langle M_{tt},\log(1+\mathcal R)\rangle|
\le[J/(8\eta^2)]\delta D^2.
\]
Thus
\[
F_{tt}\le-\frac18D^2+A_q\delta D^2\le-\frac1{16}D^2.
\]
This proves the fixed-\(q\) statement.

The Hilbert-valued regularity established in Section 3 and the corresponding elementary linear regularity imply \(\mathcal R\in C^2([0,1];L^2(\varphi))\). The actual conditional-expectation density ratio, rather than a Hermite truncation, supplies the uniform bounds \(1-3\delta\le1+\mathcal R\le1+3\delta\).

Here is a direct entropy differentiation justification that does not assume an \(L^3\) bound on difference quotients. Let \(\Phi(u)=(1+u)\log(1+u)\). Its first three derivatives are bounded on a fixed interval containing all the values of \(\mathcal R\). The scalar Taylor remainder bounded by \(C|\Delta\mathcal R|^2\), together with \(L^2\) differentiability, gives
\[
\frac{d}{dt}\int\Phi(\mathcal R)\varphi
=\langle\Phi'(\mathcal R),\mathcal R_t\rangle .
\]
For the second derivative, put \(W_h=[\mathcal R(t+h)-\mathcal R(t)]/h\), which tends to \(\mathcal R_t(t)\) in \(L^2\). The mean-value multiplier
\[
A_h=\int_0^1\Phi''\big(\mathcal R(t)+u[\mathcal R(t+h)-\mathcal R(t)]\big)\,du
\]
is uniformly bounded and converges in measure to \(\Phi''(\mathcal R(t))\). Therefore
\(\langle A_hW_h,\mathcal R_t(t)\rangle\) tends to
\(\int\Phi''(\mathcal R(t))\mathcal R_t(t)^2\varphi\):
use \(L^2\) convergence of \(W_h\), then bounded convergence in measure against the fixed integrable density \(\mathcal R_t(t)^2\). The term differentiating \(\mathcal R_t\) converges by \(L^2\) continuity to \(\langle\Phi'(\mathcal R),\mathcal R_{tt}\rangle\). This proves the displayed entropy curvature formula. Continuity of the resulting second derivative follows similarly because the squares of \(L^2\)-convergent first derivatives converge in \(L^1\), and are uniformly integrable. The same argument is one-sided at either endpoint.

Positive Gaussian noise gives smooth spatial representatives; writing \(\mathcal R=T_{\sqrt q}v_t\) with uniformly bounded \(v_t\) makes the spatial derivative bounds explicit. Gaussian Sobolev approximation justifies self-adjointness and (P19)--(P20). No boundedness of finite Hermite truncations is needed.

## 6. One radius for all smaller signal fractions

The fixed-signal theorem automatically implies the assertion uniform in \(0<q\le Q\), using the same \(\delta_Q\) selected at \(Q\).

Indeed set
\[
g_q=T_{\sqrt{q/Q}}g.
\]
The Gaussian semigroup preserves evenness, the mass and variance constraints, and the \(L^\infty\) bound. Composition of independent Gaussian smoothing gives
\[
F_g(q,t)=F_{g_q}(Q,t),\qquad
L_g(q,t)=L_{g_q}(Q,t).
\]
Apply the fixed-\(Q\) theorem to \(g_q\). This proves (A) with one radius on the entire interval \(0<q\le Q\).

## 7. Numerical sanity checks for the parameter recipe

The following values are 100-decimal diagnostics of the closed-form scalar recipe (P7), rounded here for display. They are not interval certificates and are not used in the proof. The checker also verifies the direct Gaussian square root against the independent Mehler parametrisation and direct kernel composition:

| \(q\) | chosen \(s\) | chosen \(\varepsilon\) | strict kernel gap | \(J\) | sufficient radius |
|---:|---:|---:|---:|---:|---:|
| 0.50 | 0.5202847075 | 0.1375827353 | 0.0473458334 | 1.4320544662 | \(1.7134\cdot10^{-5}\) |
| 0.53 | 0.5419050561 | 0.0870079716 | 0.0296201831 | 1.5902847695 | \(4.8520\cdot10^{-6}\) |
| 0.56 | 0.5633202098 | 0.0262722883 | 0.0088148014 | 2.1257516652 | \(2.5717\cdot10^{-7}\) |
| 0.57 | 0.5704164638 | 0.0033875871 | 0.0011301518 | 3.5370963892 | \(2.3569\cdot10^{-9}\) |

The shrinking radius near \(4/7\) is expected: both the interpolation margin and the extra smoothing used for derivatives vanish at that endpoint. The radius is sufficient and deliberately unoptimised.


\newpage

# Part II. Positive midpoint curvature above four sevenths

## 1. Definitions and result

Write
\[
\varphi(x)=(2\pi)^{-1/2}e^{-x^2/2},\qquad
h(p)=-\int_{\mathbb R}p(x)\log p(x)\,dx.
\]
All logarithms are natural. If \(X\) has density \(f\), define \(\mathcal S_qf\) as the density of
\[
\sqrt q\,X+\sqrt{1-q}\,G,
\]
where \(G\) is an independent standard Gaussian. The variable \(q\) is the **signal variance fraction**. For two independent variables of density \(\mathcal S_qf\), put
\[
F_{q,f}(t)=h\bigl(\sqrt t\,X_q+\sqrt{1-t}\,Y_q\bigr),
\quad 0\le t\le1.
\]
The variance remains one, and exchange of the copies gives \(F_{q,f}(t)=F_{q,f}(1-t)\).

The original concavity question is due to Ball, Nayar and Tkocz [1]. Luo [2] already constructed symmetric counterexamples for every \(q>4/7\). His Problem 10.3 asks whether the midpoint can have strictly positive second derivative throughout the interval not covered by his Appendix E. The result below concerns that local shape, and also gives its asymptotic scale. It is not a first refutation of the original conjecture.

**Theorem B (midpoint curvature and local profile).** Fix an integer \(K\ge1\). There is a sequence \(f_n\), indexed by sufficiently large even integers, of strictly positive, smooth, even, variance-one densities with the following properties.

1. Every \(f_n\) matches the standard Gaussian moments through degree \(2K\).
2. Both \(f_n\) and \(\mathcal S_qf_n\) are strongly log-concave for sufficiently large \(n\), uniformly for \(q\in[0,1]\).
3. For each fixed nonnegative integer \(m\), their density ratios converge to one in \(C_b^m(\mathbb R)\), uniformly in \(q\in[0,1]\).
4. For every fixed \(q>4/7\),
\[
F_{q,f_n}''(1/2)>0
\]
for all sufficiently large \(n\). More precisely, with the notation below,
\[
F_{q,f_n}''(1/2)
\sim4n\gamma(q)\epsilon_n^4\|M_{q,n,1/2}\|_2^2,
\tag{M1}
\]
where
\[
\epsilon_n=e^{-\sqrt n},\qquad
\gamma(q)=\frac{\kappa+b}{\kappa-b}
=\frac{21\sqrt5\,q}{49q-4}.
\tag{M2}
\]
All norms in this article, except those explicitly identified as Lebesgue norms, use \(\varphi(x)\,dx\).
5. If
\[
D_{q,n}=h(\varphi)-F_{q,f_n}(1/2),
\]
then \(D_{q,n}>0\) eventually and
\[
\frac{F_{q,f_n}(1/2+s/\sqrt n)-F_{q,f_n}(1/2)}
     {D_{q,n}}
\longrightarrow 1-e^{-4\gamma(q)s^2}
\tag{M3}
\]
in \(C^2\) on each compact interval of the variable \(s\).

The eventual statements and convergence are uniform for \(q\) in compact subintervals of \((4/7,1]\). The same sequence of densities works for all these fixed \(q\). No claim at \(q=4/7\) follows.

Equation (M3) specifies the local geometry: the entropy dip has width of order \(n^{-1/2}\). On this scale its curvature changes sign at \(|s|=1/\sqrt{8\gamma(q)}\) in the limiting profile. This is an asymptotic statement about the constructed sequence, not a universal description of all symmetric inputs.

## 2. The inherited construction and its estimates

We record the construction from [2] explicitly so that the new argument has a fixed referent. The probabilists' Hermite polynomials satisfy
\[
e^{xz-z^2/2}=\sum_{j\ge0}H_j(x)\frac{z^j}{j!}.
\]
Set
\[
a=\frac7{3\sqrt5},\qquad
b=\frac{a-1}{a+1}=\frac{47-21\sqrt5}{2},\qquad
h_n(x)=\sqrt2e^{-ax^2/2}
          \frac{H_n(\sqrt{2a}\,x)}{\sqrt{n!}}.
\tag{M4}
\]
The parameter \(b\) is positive and less than \(1/40\).

For \(0\le j\le K\), define
\[
\psi_j(x)=\sqrt2e^{-x^2/2}
          \frac{H_{2j}(\sqrt2x)}{\sqrt{(2j)!}}.
\]
Choose \(c_{n,j}\) so that
\[
\widetilde h_n=h_n-\sum_{j=0}^Kc_{n,j}\psi_j,\qquad
\int x^{2\ell}\widetilde h_n(x)\varphi(x)\,dx=0
\quad(0\le\ell\le K),
\tag{M5}
\]
and put
\[
f_n=\varphi(1+\epsilon_n\widetilde h_n).
\tag{M6}
\]

For completeness, the correction exists because
\[
\int x^{2\ell}\psi_j(x)\varphi(x)\,dx
=\frac{2^{-\ell}}{\sqrt{(2j)!}}
  \int y^{2\ell}H_{2j}(y)\varphi(y)\,dy
\]
is zero when \(\ell<j\) and equals \(2^{-j}\sqrt{(2j)!}\) when \(\ell=j\). Thus the coefficient matrix in (M5) is triangular and invertible.

Mehler's identity bounds the supremum of \(h_n\) by
\[
B_n=\sqrt{2e}(n+1)^{1/4}.
\]
The derivative recurrence
\[
h_{n,a}'=\sqrt{a/2}
 \left(\sqrt n\,h_{n-1,a}-\sqrt{n+1}\,h_{n+1,a}\right)
\tag{M7}
\]
shows that every fixed spatial derivative grows at most polynomially in \(n\). Integrating the Hermite generating function against \(\varphi\) gives
\[
\int\varphi(x)e^{-ax^2/2}
 e^{\sqrt{2a}xz-z^2/2}\,dx
=(1+a)^{-1/2}e^{bz^2/2}.
\]
Insertion of \(x^{2\ell}\) multiplies the right side by a fixed polynomial of degree \(2\ell\). Coefficient extraction, followed by the fixed triangular inverse, gives, for each fixed \(m\),
\[
\|\widetilde h_n-h_n\|_{C_b^m}
\le C_{K,m}n^K b^{n/2}.
\tag{M8}
\]
These estimates prove the density, moment and approximation assertions in Theorem B. In detail, \(\epsilon_n\|\widetilde h_n\|_{C_b^m}\to0\) for every fixed \(m\), so positivity follows, and
\[
(\log f_n)''
\le -1+\frac{\epsilon_n\|\widetilde h_n''\|_\infty}
              {1-\epsilon_n\|\widetilde h_n\|_\infty}
<-1/2
\]
eventually. Gaussian smoothing contracts the supremum of each derivative of a relative perturbation, and preserves the Gaussian moment matches. It gives the same estimates uniformly for \(q\in[0,1]\).

Let
\[
T_\rho g(x)=\mathbb E\,g(\rho x+\sqrt{1-\rho^2}\,G)
\]
be the Ornstein--Uhlenbeck operator. The output density of the weighted sum is exactly
\[
p_{q,n,t}
=\varphi(1+\epsilon_n\widetilde L_t+\epsilon_n^2\widetilde M_t),
\tag{M9}
\]
where
\[
L_t=T_{\sqrt{qt}}h_n+T_{\sqrt{q(1-t)}}h_n
\]
and
\[
M_t(x)=\mathbb E\left[
h_n(G_1)h_n(G_2)\mid
\sqrt q\,(\sqrt t\,G_1+\sqrt{1-t}\,G_2)
+\sqrt{1-q}\,G_3=x\right].
\tag{M10}
\]
Tildes mean that \(h_n\) is replaced by \(\widetilde h_n\).

The linear generating function is
\[
\sum_{j\ge0}z^j\|T_{\sqrt r}h_j\|_2^2
=\frac{2}
{\sqrt{(1+a)^2-a^2r^2}}
\frac1{\sqrt{(1-\tau(r)z)(1-\nu(r)z)}},
\tag{M11}
\]
where
\[
\tau(r)=\frac{1-a+ar}{1+a-ar},\qquad
\nu(r)=\frac{a-1+ar}{1+a+ar}.
\]
This identity follows by applying Mehler's identity to the two correlated Gaussian arguments in the squared norm. For \(r>2/7\), \(0<\nu(r)<\tau(r)\le1\). The binomial coefficients of \((1-z)^{-1/2}\) are positive and their convolution sums to one; hence
\[
\|T_{\sqrt r}h_n\|_2^2\le C_r\tau(r)^n.
\tag{M12}
\]

Define
\[
\lambda_\pm=\frac{a(1\pm q)-1}{a(1\pm q)+1},\qquad
\alpha=\frac{\lambda_++\lambda_-}{2}-b,\qquad
\beta=\frac{\lambda_+-\lambda_-}{2},\qquad
\kappa=-\lambda_-.
\tag{M13}
\]
For \(q\ge4/7\), put
\[
d=-b-\alpha/2,\qquad e=-\alpha/2,\qquad c=\beta/2.
\tag{M14}
\]
All three are positive. For example,
\[
-\alpha=\frac{2a^2q^2}
 {(1+a)((1+a)^2-a^2q^2)}
\]
increases with \(q\), and
\[
d(4/7)=\frac{203\sqrt5-453}{38}>0.
\]
The last numerator is positive because \(5\cdot203^2-453^2=836>0\).
We have
\[
e-d=b,\qquad d+e+2c=\kappa.
\tag{M15}
\]

Direct substitution in (M11)--(M13) gives the strict exponential comparison
\[
\tau(q/2)<\kappa^2\qquad(q>4/7),
\tag{M16}
\]
with equality at \(q=4/7\). One way to verify the sign is to clear positive denominators; the remaining factor is
\[
-\frac{a}{45}(7q-4)(7q^2-24q+2),
\]
over the positive denominator \((1+a-aq)^2(1+a-aq/2)\). The numerator is positive on \((4/7,1]\). Also \(\sqrt b<\kappa\) on \([4/7,1]\). These are the two margins used below.

## 3. A differentiated coefficient lemma

The following elementary saddle calculation is the main additional analytic step. The method is standard; the derivative control and its application to this coefficient are stated in full.

**Lemma 2 (positive quadratic coefficient).** Let \(B(u)\) be a smooth real symmetric \(4\times4\) matrix, with every entry strictly positive near \(u=0\). Set
\[
C_n(u)=[z_1^nz_2^nz_3^nz_4^n]\,
          e^{z^{\mathsf T}B(u)z/2}.
\]
Let \(r_i(u)>0\) solve
\[
r_i(B(u)r)_i=1\quad(1\le i\le4),
\]
and put
\[
\Phi(u)=2-\sum_i\log r_i(u),\qquad
H(u)=I+\operatorname{diag}(r)B(u)\operatorname{diag}(r).
\]
Then, locally uniformly in \(u\),
\[
C_n(u)\sim
\frac{2}{(2\pi n)^2\sqrt{\det H(u)}}\,
n^{-2n}e^{n\Phi(u)},
\tag{M17}
\]
and
\[
(\log C_n)'=n\Phi'+O(1),\qquad
(\log C_n)''=n\Phi''+O(1).
\tag{M18}
\]
Uniformity also holds for compact parameter families with a common positive lower bound on the entries.

*Proof.* Minimise
\[
f(u,y)=\frac12\sum_{i,j}B_{ij}(u)e^{y_i+y_j}-\sum_i y_i.
\]
The positive diagonal terms make its Hessian positive definite. It is coercive: if the largest coordinate tends to \(+\infty\), its diagonal exponential dominates the linear term; if the largest coordinate remains bounded while \(|y|\to\infty\), the negative sum tends to \(+\infty\). The unique minimiser gives \(r_i=e^{y_i}\). The implicit-function theorem gives smooth dependence on \(u\), and the Hessian at the minimiser is \(H\).

On the coefficient contour \(z_i=\sqrt n\,r_i e^{i\theta_i}\), Cauchy's formula becomes
\[
C_n=\frac{n^{-2n}e^{n\Phi}}{(2\pi)^4}I_n,\qquad
I_n=\int_{\mathbb T^4}e^{n\Psi(u,\theta)}\,d\theta,
\]
where
\[
\Psi=\frac12\sum_{i,j}B_{ij}r_ir_j
   \bigl(e^{i(\theta_i+\theta_j)}-1\bigr)-i\sum_i\theta_i.
\]
Its real part is nonpositive. Equality requires every
\(\theta_i+\theta_j\) to vanish modulo \(2\pi\), including \(i=j\). The only maxima are \(0\) and \(\pi\mathbf1\). Their contributions have equal phase because the extracted degree is \(4n\).

At either maximum, \(\Psi_\theta=0\) and
\(\Psi_{\theta\theta}=-H\). Gaussian rescaling gives
\[
I_n\sim\frac{2(2\pi)^2}{n^2\sqrt{\det H}}.
\tag{M19}
\]
The two maxima are uniformly nondegenerate; away from small neighbourhoods of them the real part is bounded above by a strictly negative constant.

It remains to justify differentiation. Since
\(\Psi(u,0)=0\) and \(\Psi_\theta(u,0)=0\) identically in \(u\),
\[
|\Psi_u|+|\Psi_{uu}|\le C|\theta|^2
\]
near each maximum. The differentiated integrands are
\(n\Psi_u e^{n\Psi}\) and
\((n\Psi_{uu}+n^2\Psi_u^2)e^{n\Psi}\).
Their absolute integrals are \(O(n^{-2})\), by the same rescaling and a Gaussian integrable majorant. The rest of the torus contributes exponentially less. Thus
\[
I_n',I_n''=O(n^{-2}),\qquad I_n\asymp n^{-2}.
\]
It follows that the first two derivatives of \(\log I_n\) are \(O(1)\), proving (M18). No uncontrolled asymptotic remainder has been differentiated. \(\square\)

## 4. The coefficient for the quadratic convolution

Two conditionally independent copies of \((G_1,G_2)\) in (M10) have covariance
\[
\begin{pmatrix}I_2&qv_tv_t^{\mathsf T}\\
qv_tv_t^{\mathsf T}&I_2\end{pmatrix},
\qquad v_t=(\sqrt t,\sqrt{1-t})^{\mathsf T}.
\]
A Gaussian integral of the four Hermite generating functions, followed by the substitution
\((z_1,z_2,z_3,z_4)\mapsto(iz_1,iz_2,-iz_3,-iz_4)\), gives
\[
\|M_t\|_2^2=\frac{4(n!)^2}{\Delta}\,C_n(t-1/2),
\quad
\Delta=(1+a)\sqrt{(1+a)^2-a^2q^2}.
\tag{M20}
\]
The coefficient matrix is
\[
B(u)=
\begin{pmatrix}
d+2eu&ew&c(1+2u)&cw\\
ew&d-2eu&cw&c(1-2u)\\
c(1+2u)&cw&d+2eu&ew\\
cw&c(1-2u)&ew&d-2eu
\end{pmatrix},
\quad w=\sqrt{1-4u^2}.
\tag{M21}
\]
The coefficient is unchanged by the phase substitution, since its phase is
\(i^{2n}(-i)^{2n}=1\). Equations (M20)--(M21) and their multigraph interpretation are already present in [2]. We use them to evaluate the derivative.

**Lemma 3 (explicit saddle and curvature).** Put \(U=\kappa+b\). For the matrix (M21),
\[
\frac{C_n''(0)}{C_n(0)}
=-8n\frac{\kappa+b}{\kappa-b}+O(1).
\tag{M22}
\]
Moreover,
\[
\|M_{1/2}\|_2^2
\sim \frac{4}{\pi\Delta\sqrt{\det H_0}}\,
          \frac{\kappa^{2n}}n,
\quad
H_0=I+B(0)/\kappa,
\tag{M23}
\]
and
\[
\det H_0=
2\left(1-\frac{\lambda_+}{\kappa}\right)
 \left(1-\frac b\kappa\right)^2.
\tag{M24}
\]

*Proof.* At \(u=0\) every row of \(B\) sums to \(\kappa\), so
\(r=\kappa^{-1/2}\mathbf1\). For \(v=(1,-1,1,-1)^{\mathsf T}\),
\[
B(0)v=-bv,\quad B'(0)\mathbf1=Uv,\quad
\mathbf1^{\mathsf T}B''(0)\mathbf1=-8U.
\]
The Hessian formula for the minimum of \(f\) gives
\[
\Phi''(0)=f_{uu}-f_{uy}H_0^{-1}f_{yu}
=-\frac{4U}{\kappa}
 -\frac{4U^2}{\kappa(\kappa-b)}
=-8\frac{\kappa+b}{\kappa-b}.
\]
The coefficient is even in \(u\), so its first derivative at zero vanishes. Lemma 2 proves (M22). The eigenvalues of \(B(0)\) are
\(\kappa,-\lambda_+,-b,-b\), giving (M24). Strict entrywise positivity ensures \(\det H_0>0\). Combining (M17), (M20) and Stirling's formula proves (M23). \(\square\)

There is also an explicit expression for the saddle away from equal weights. Symmetry forces \(r=(x,y,x,y)\). Set
\[
A=-b+Ut,\quad D=-b+U(1-t),\quad E=U\sqrt{t(1-t)}.
\]
The saddle equations are \(Ax^2+Exy=Dy^2+Exy=1\). Consequently
\[
xy=\Lambda(t)^{-1},\quad
\Lambda(t)=\sqrt{U^2t(1-t)-b\kappa}
             +U\sqrt{t(1-t)},
\quad
\Phi(t)=2+2\log\Lambda(t).
\tag{M25}
\]
These formulas hold in a fixed neighbourhood of \(1/2\), where the square roots and entries are positive. They give
\[
\Lambda(1/2)=\kappa,\qquad
(\log\Lambda)''(1/2)=-4\gamma(q).
\]
Equations (M17)--(M18), after setting \(t=1/2+s/\sqrt n\), therefore imply
\[
\frac{\|M_{1/2+s/\sqrt n}\|_2^2}
     {\|M_{1/2}\|_2^2}
\longrightarrow e^{-4\gamma(q)s^2}
\tag{M26}
\]
in local \(C^2\). To see the derivative assertion directly, the first derivative of the residual logarithmic factor is \(O(1)/\sqrt n\) after this change of variable, and its second derivative is \(O(1)/n\). The exponential part is an ordinary Taylor expansion of the smooth function \(\log\Lambda\).

## 5. Derivatives of the density on the correct scale

A scalar norm asymptotic alone does not control the entropy remainder. We need bounds for \(M_t'\) and \(M_t''\) in the same Gaussian \(L^2\) space.

**Lemma 4 (two-weight kernel bound).** On a sufficiently small fixed neighbourhood \(J\) of \(1/2\), and for \(j=0,1,2\),
\[
\|\partial_t^jM_t\|_2\le C(1+n)^j\|M_t\|_2.
\tag{M27}
\]

*Proof.* Define \(G_n(t,s)=\langle M_t,M_s\rangle\). Its four-variable Gaussian covariance has off-diagonal block \(qv_tv_s^{\mathsf T}\); its eigenvalues are \(1+q,1-q,1,1\). Thus the determinant prefactor in its generating function remains \(\Delta\), and its phase-changed coefficient matrix has diagonal blocks
\[
-bI_2-\alpha v_tv_t^{\mathsf T},
\quad
-bI_2-\alpha v_sv_s^{\mathsf T},
\]
and off-diagonal block \(\beta v_tv_s^{\mathsf T}\).
Every entry is strictly positive on a neighbourhood of \((1/2,1/2)\).

The coefficient of total degree \(4n\) in the exponential is the corresponding coefficient in \(Q^{2n}/(2n)!\), where \(Q\) is its quadratic form. Its expansion is a sum of positive monomials in the matrix entries. On the chosen neighbourhood, logarithmic derivatives of each entry through order four are bounded. Differentiating any monomial a total of \(k\le4\) times bounds its absolute derivative by \(C_k(1+n)^k\) times the original monomial. Summing gives
\[
|\partial_t^j\partial_s^jG_n(t,s)|
\le C_j(1+n)^{2j}G_n(t,s).
\]
At \(s=t\), the left side equals \(\|\partial_t^jM_t\|_2^2\), proving (M27). \(\square\)

Fix \(q>4/7\). By (M16), choose \(r_*>q/2\), sufficiently close to \(q/2\), so that
\[
\rho:=\sqrt{\tau(r_*)}<\kappa.
\]
Reduce \(J\) if needed so that \(q\max(t,1-t)<r_*\) for \(t\in J\). Hermite diagonalisation of \(T_\rho\) and the strict margin between these parameters absorb every fixed power of the Hermite degree. Equations (M11)--(M12) then imply
\[
\sup_{t\in J}\|\partial_t^j L_t\|_2\le C_j\rho^n,
\quad 0\le j\le2.
\tag{M28}
\]
This argument is valid although \(h_n\) has an infinite expansion in the standard Gaussian Hermite basis: a factor \(k^m r^k\) is bounded by a constant times \(r_*^k\) whenever \(r<r_*\).

The moment correction also permits weight differentiation. Conditional on the sum in (M10), one may represent the Gaussian pair as
\[
(G_1,G_2)=
\sqrt q\,v_t x+\sqrt{1-q}\,v_t Z+v_t^\perp W,
\]
where \(Z,W\) are independent standard Gaussians. On \(J\), the first two derivatives of \(v_t\) are bounded. Differentiate this representation under the expectation and use (M7)--(M8), then integrate the resulting fixed Gaussian moments. For some fixed exponent \(A\),
\[
\sup_{t\in J}
\left(
\|\partial_t^j(\widetilde L_t-L_t)\|_2+
\|\partial_t^j(\widetilde M_t-M_t)\|_2
\right)
\le C n^A b^{n/2},
\quad j\le2.
\tag{M29}
\]
This representation remains valid at \(q=1\); no inverse of a degenerate covariance is used.

## 6. Entropy transfer and proof of Theorem B

Write the output relative perturbation as
\[
R_t=\epsilon_n\widetilde L_t+\epsilon_n^2\widetilde M_t.
\]
Conditional expectation and the supremum bounds in Section 2 give, uniformly in \(t\),
\[
\|R_t\|_\infty\le
2\epsilon_n\|\widetilde h_n\|_\infty
+\epsilon_n^2\|\widetilde h_n\|_\infty^2
=:\eta_n,
\tag{M30}
\]
where \(n^m\eta_n\to0\) for every fixed \(m\).

Fix a compact interval of the scaled variable \(s\), and set
\(t=1/2+s/\sqrt n\), \(A_n=\|M_{1/2}\|_2\).
Equations (M23), (M26)--(M29), together with
\(\rho<\kappa\) and \(\sqrt b<\kappa\), imply
\[
R_t=\epsilon_n^2M_t+E_n(s),\qquad
\|\partial_s^jE_n(s)\|_2=o(\epsilon_n^2A_n)
\quad(j=0,1,2).
\tag{M31}
\]
The little-oh errors remain negligible after multiplication by any fixed power of \(n\). Indeed their ratios are bounded by a fixed polynomial times
\(e^{\sqrt n}(\rho/\kappa)^n\) or
\(e^{\sqrt n}(\sqrt b/\kappa)^n\).

Lemma 4 and (M26) also give
\[
\|R_t\|_2=O(\epsilon_n^2A_n),\quad
\|\partial_sR_t\|_2=O(n^{1/2}\epsilon_n^2A_n),\quad
\|\partial_s^2R_t\|_2=O(n\epsilon_n^2A_n).
\tag{M32}
\]
These polynomial bounds are sufficient; sharp derivative norms are unnecessary.

For a mean-zero relative perturbation \(R\), define
\[
\mathcal D(R)=\int\varphi\{(1+R)\log(1+R)-R\}.
\]
If \(\|R\|_\infty\le\eta<1/2\), Taylor's theorem gives
\[
\left|\mathcal D(R)-\tfrac12\|R\|_2^2\right|
\le C\eta\|R\|_2^2.
\]
Differentiating along a twice differentiable parameter curve yields errors bounded respectively by
\[
C\eta\|R\|_2\|R'\|_2
\quad\hbox{and}\quad
C\eta\bigl(\|R'\|_2^2+\|R\|_2\|R''\|_2\bigr).
\tag{M33}
\]
These follow directly from
\[
\mathcal D'(R)=\langle\log(1+R),R'\rangle,\qquad
\mathcal D''(R)=
\langle\log(1+R),R''\rangle+
\int\frac{(R')^2}{1+R}\varphi.
\]
For the present finite-degree Gaussian-polynomial densities, differentiation is justified by Gaussian envelopes and the uniform positive lower bound on \(1+R\).

Since the output variance is one,
\(h(\varphi)-F_{q,f_n}(t)=\mathcal D(R_t)\).
Equations (M30)--(M33) show that
\[
\frac{h(\varphi)-F_{q,f_n}(1/2+s/\sqrt n)}
     {\tfrac12\epsilon_n^4 A_n^2}
\longrightarrow e^{-4\gamma(q)s^2}
\tag{M34}
\]
in local \(C^2\). The second derivative error in (M33), after this normalisation, is at most \(C n\eta_n+o(1)\), which tends to zero. This is the point at which bounds on the full derivatives of \(M_t\), rather than only on derivatives of its squared norm, are necessary.

At \(s=0\), (M34) proves
\(D_{q,n}\sim\tfrac12\epsilon_n^4 A_n^2>0\).
Dividing by this value gives (M3). Taking the second \(s\)-derivative at zero, and using \(d^2/ds^2=n^{-1}d^2/dt^2\), gives (M1). Its right side is positive, so exchange symmetry makes the midpoint a strict local minimum. This completes the proof.

The fully explicit leading term is
\[
F_{q,f_n}''(1/2)
\sim
\frac{16\gamma(q)}
{\pi\Delta\sqrt{\det H_0}}\,
e^{-4\sqrt n}\kappa^{2n}.
\tag{M35}
\]
The entropy dip therefore becomes extremely small as \(n\) grows. The theorem concerns its sign and local shape, not a large practical entropy loss.

## 7. Verification and limits

The proof is analytic and has not been checked by a proof assistant. A separate implementation evaluates the coefficient and its first two weight derivatives by a Gaussian-moment recurrence. Exact rational tests also enumerate all admissible four-vertex multigraph configurations independently.

Thirty numerical cases use \(q=4/7,0.575,0.6,0.75,1\) and \(n=2,4,8,12,16,24\). They support the coefficient asymptotic; they do not establish the unrestricted theorem. At \(n=24\), for example, \(C_n''/(nC_n)\) is \(-8.80988058680\) for \(q=0.6\), against the limiting value \(-8.87384457213\).

The \(q=4/7\) computations concern only the generating coefficient. They do not determine the entropy-curvature sign at that signal fraction: the strict linear/quadratic comparison in (M16) is absent there.

A useful exact check is \(b=0\). Then diagonal rescaling gives
\[
C_n(u)=(1-4u^2)^nC_n(0),
\]
so \(C_n''(0)/C_n(0)=-8n\) exactly, agreeing with (M22).

The construction and the exponential comparison in (M16) come from Luo. The added results are the differentiated coefficient asymptotic, explicit prefactor, derivative-scale estimates and entropy profile. No earlier matching resolution was found in the bounded search accompanying this manuscript. That search does not certify historical priority.

The theorem does not settle the uniform concavity question for all symmetric log-concave inputs. It does not settle the boundary \(q=4/7\). Positive local results of Schubr [3] and Abu-Ajamieh [4] concern different input or signal ranges and must be distinguished from the counterexample construction here.



# References

1. Ball, K., Nayar, P., & Tkocz, T. (2016). A reverse entropy power inequality for log-concave random vectors. *Studia Mathematica, 235*(1), 17--30. https://doi.org/10.4064/sm8418-6-2016

2. Luo, C. (2026). *Near-Gaussian counterexamples to the Ball--Nayar--Tkocz entropy concavity conjecture* (Version 2) [Preprint]. arXiv. https://arxiv.org/abs/2609.11418v2

3. Schubr, S. (2026, October 30). *Mixing-weight entropy concavity near the Gaussian below one-half signal* (Public edition 3) [Preprint]. Luminary Research Archive. https://research.luminaryoddities.com/local-and-eventual-entropy-concavity-under-gaussian-smoothing/v/3

4. Abu-Ajamieh, F. (2026). *Noise-induced restoration of entropy concavity near Gaussian equilibrium* (Version 1) [Preprint]. arXiv. https://arxiv.org/abs/2610.22386v1

5. Abbe, E., & Zheng, L. (2012). A coordinate system for Gaussian networks. *IEEE Transactions on Information Theory, 58*(2), 721--733. https://arxiv.org/abs/1002.2271


# Research and verification status

This manuscript was developed in an AI-assisted mathematical research session. Separate agents derived the saddle constant, wrote an independent analytic lemma, searched the literature and implemented coefficient checks. The accompanying review is an internal independent-agent review, not external peer review. No paper has been submitted or published by this session.


---
title: "The critical signal fraction for local entropy concavity"
subtitle: "Gaussian noise restores mixing-weight concavity near the Gaussian exactly below signal fraction 4/7"
author: "Anonymous"
date: "3 November 2026 (revised version; first internal version 2 November 2026)"
lang: en-GB
fontsize: 10pt
geometry:
  - margin=23mm
colorlinks: true
linkcolor: MidnightBlue
urlcolor: MidnightBlue
header-includes:
  - \usepackage{xurl}
---

# Abstract

Take two independent copies of a random signal, add Gaussian noise to each, and mix the two noisy copies with weights \(\sqrt t\) and \(\sqrt{1-t}\). How does the entropy of the mixture depend on the mixing weight \(t\)? Ball, Nayar and Tkocz conjectured that, for log-concave signals, it is always a concave function of \(t\). Luo has since disproved this, and showed that symmetric near-Gaussian counterexamples survive added noise whenever the signal accounts for more than \(4/7\) of the total variance. We show that \(4/7\) is the true local threshold for symmetric signals. For every signal fraction \(q<4/7\), there is an explicit neighbourhood of the Gaussian, defined by a uniform bound on the relative density \(f/\varphi-1\), in which every symmetric variance-one input gives a strictly concave entropy curve (Theorem A). Combined with Luo's counterexamples, this identifies \(4/7\) as the supremum of the signal fractions with this property; in noise terms, it corresponds to a signal-to-noise ratio of \(4/3\). Whether \(q=4/7\) itself has the property remains open. Above the threshold we sharpen the counterexample side: for every \(q\in(4/7,1]\), one sequence of smooth, strongly log-concave inputs makes equal weighting a strict local entropy minimum, with an explicit curvature constant and a Gaussian-shaped local profile (Theorem B). We also explain why the threshold is \(4/7\): it is the fixed point \(q=b(q/2)\) of an Ornstein--Uhlenbeck smoothing inequality whose range we prove to be sharp, and, within the Hermite test functions used to prove sharpness, the extremal member is exactly Luo's counterexample family (Proposition S). Finally, at \(q=4/7\) Luo's family has strictly negative entropy curvature near equal weighting (Proposition E), so it gives no midpoint counterexample at the boundary. The positive-curvature sign on part of the range, Luo's Problem 10.3, was also proposed in an earlier public, unrefereed argument; our additions there are the exact constant, the profile and a complete written proof.

# 1. Introduction

## 1.1 The question in plain terms

Mixing independent signals usually makes them look more Gaussian [2], and Gaussian distributions have the largest entropy for a given variance. A natural guess is therefore that, as one moves the mixing weight from "all of one copy" to "an even blend", entropy rises smoothly and without wiggles. The Ball--Nayar--Tkocz conjecture makes this precise. For independent \(X,Y\) with a common log-concave density, it asserts that
\[
t\longmapsto h\bigl(\sqrt t\,X+\sqrt{1-t}\,Y\bigr)
\]
is concave on \([0,1]\) [1, Conjecture 2]. Here \(h\) is differential entropy. The weights keep the variance fixed. Lieb's form of the entropy power inequality [3] already gives the chord inequality \(h(\sqrt tX+\sqrt{1-t}Y)\ge th(X)+(1-t)h(Y)\); concavity would strengthen it. The conjecture is known for Gaussian mixtures [5], building on entropy comparisons for Gaussian mixtures in [4].

Luo [6] disproved the conjecture. He also asked a sharper question: how much independent Gaussian noise must be added to each copy before concavity returns, at least near the Gaussian? Write \(q\) for the fraction of the variance that is signal, so that each copy is \(\sqrt q\,X+\sqrt{1-q}\,G\), and \(1-q\) is the noise fraction. Luo's counterexamples cover every \(q>4/7\) for symmetric inputs, and every \(q>0\) for asymmetric ones. Luo's original construction, his Problems 10.1--10.3 and the further results we rely on are in [6]; all four 2026 sources [6--9] are unrefereed. Schubr [7] then proved that small symmetric perturbations of the Gaussian cannot destroy concavity when \(q<1/2\). This left the symmetric range \(1/2\le q\le4/7\) open. Luo described \(4/7\) as the boundary of an auxiliary estimate in his method, not as an established critical value.

This paper closes the open interval, except for its right endpoint, and explains where \(4/7\) comes from.

**In one sentence.** Near the Gaussian, Gaussian noise restores entropy concavity for symmetric inputs whenever the noise variance exceeds three quarters of the signal variance (\(q<4/7\)); at every higher signal fraction (\(q>4/7\)), small smooth log-concave bumps can create a genuine local entropy dip at equal weighting.

## 1.2 What is new, and what is inherited

The following list separates our contribution from what we take from earlier work.

| Item | Status | Source of the key ingredient |
|---|---|---|
| Theorem A: local concavity for every \(q<4/7\) | New. It extends Schubr's \(q<1/2\) theorem [7] | Our weighted, bounded-input version (Lemma I.2) of Luo's Gaussian-kernel interpolation [6, Theorem F.1, Corollary F.2]. It replaces the spectral Cauchy--Schwarz step that limits Schubr's argument to \(1/2\). Schubr's entropy-energy argument then absorbs the logarithmic terms |
| Corollary: \(\sup\mathcal C=4/7\) | New. It turns Luo's "method boundary" into the supremum | Theorem A for the lower bound; Luo's counterexamples [6, Theorem 1.2] for the upper bound |
| Theorem B: strict midpoint local minima for every \(q\in(4/7,1]\) | The sign was known for \(q>2-\sqrt2\approx0.586\) [6, Appendix E]. On \((4/7,2-\sqrt2]\), which is Luo's Problem 10.3, the sign was first proposed in an unrefereed public argument dated 1 November 2026 [8]. Here: a complete written proof | Luo's construction and coefficient representation [6]. New here: the exact curvature constant \(\gamma(q)\), one input sequence for all \(q\in(4/7,1]\) (including \(q=1\), which [8] excludes) and the \(C^2\) local profile |
| Proposition S: the interpolation range is sharp | New | Luo's test functions [6, Appendix F] |
| Proposition E: at \(q=4/7\), Luo's family is strictly concave near \(t=1/2\) | New | Parts I and II |
| Remark R: \(\delta_Q\asymp(4/7-Q)^{9/4}\) for the proof's radius | New | Explicit constants of Part I |

Why does Luo's interpolation inequality not already give Theorem A? It has three gaps. First, it is stated for perturbations with finite Lebesgue \(L^2(dx)\) norm, and a small uniform bound \(\|g\|_\infty\) does not make that norm small, or even finite. Second, it compares sizes at a single weight; it does not control the derivatives in \(t\) that concavity requires. Third, the logarithm in the entropy produces nonlinear terms that must be absorbed. Lemma I.2 handles the first gap with a weak Gaussian weight. Extra Gaussian smoothing handles the second. Schubr's Hermite-energy method handles the third.

## 1.3 Why the threshold is four sevenths

Write the density of the noisy mixture as \(\varphi\,(1+L+M)\), where \(\varphi\) is the standard Gaussian density and \(f=\varphi(1+g)\) is the input density. The term \(L\) is linear in the perturbation \(g\). Its Hermite component of degree \(k\) is damped by \((qt)^{k/2}+(q(1-t))^{k/2}\), a convex function of \(t\). This term always pushes the entropy curve towards concavity. The term \(M\) is quadratic in \(g\). It records the interaction of the two copies and can push the other way.

The question is therefore whether the quadratic interaction can ever outweigh the linear term. At equal weighting, the linear term sees the input smoothed to signal level \(r=q/2\). The quadratic term behaves like the input smoothed to a level \(s\) just above \(q\). For bounded inputs, Lemma I.2 and Proposition S show that
\[
\|T_{\sqrt s}g\|_2^2\le J\,\|g\|_\infty\,\|T_{\sqrt r}g\|_2
\quad\text{holds for } s<b(r)=\sqrt{2r(1+r)}-r,
\quad\text{and fails for } s>b(r).
\]
Here \(T_\rho\) is Gaussian (Ornstein--Uhlenbeck) smoothing. With \(r=q/2\) and \(s\) slightly above \(q\), the comparison is available exactly when \(q<b(q/2)\), that is, when \(q<4/7\). The value \(4/7\) is the fixed point of \(q=b(q/2)\). At that point, within the Hermite test functions \(h_{n,a}\) used to prove sharpness, the extremal member is exactly Luo's counterexample family. The positive and negative sides of the threshold are therefore two views of one Gaussian smoothing inequality, not a coincidence of two methods.

This mechanism concerns high Hermite degrees. A single fixed low-degree perturbation can tolerate more signal. For example, Abu-Ajamieh [9] proves concavity for the fixed direction \(\varphi(1+a\,\mathrm{He}_4)\), \(0<a\le10^{-2}\), up to \(q=5/8\). That does not conflict with \(\sup\mathcal C=4/7\), which concerns uniform neighbourhoods containing all directions. Moreover, \(\mathrm{He}_4\) is unbounded, so that family lies in no relative \(L^\infty\) ball.

![The regimes established here (top) and the limiting shape of the entropy dip above the threshold at \(q=3/5\) (bottom). The lower curve is the analytic profile (D) of Theorem B, not a measured finite-density effect.](figures/entropy_transition.pdf){width=82%}

## 1.4 How to read this paper

Section 2 fixes notation and states every result. Part I (Sections I.1--I.7) proves Theorem A. Part II (Sections II.1--II.7) proves Theorem B. Part III (Sections III.1--III.3) proves Propositions S and E and derives Remark R. The Discussion (Sections D.1--D.5) covers the role of each assumption, why the earlier method stopped at \(1/2\), a reading in blind source separation, and open problems. Equations carry the prefix of their part: P (Part I), M (Part II), S, E and R (Part III).

A reader who wants only the ideas can read Sections 1--2, the roadmap at the start of each Part, and the Discussion.

# 2. Setup and main results

## 2.1 Notation

Let \(\varphi(x)=(2\pi)^{-1/2}e^{-x^2/2}\) and let \(h(p)=-\int p\log p\,dx\) be differential entropy, with natural logarithms. Let \(X,Y\) be independent with a common even density \(f\) of mean zero and variance one, and let \(G\) be an independent standard Gaussian. For a **signal fraction** \(q\in(0,1]\), put
\[
F_{q,f}(t)=h\!\left(\sqrt{qt}\,X+\sqrt{q(1-t)}\,Y+\sqrt{1-q}\,G\right),\qquad 0\le t\le1.
\]
This equals the entropy obtained by first smoothing each input separately, \(X_q=\sqrt qX+\sqrt{1-q}G_1\), and then forming \(\sqrt tX_q+\sqrt{1-t}Y_q\). Every variable has variance one, and \(F_{q,f}(t)=F_{q,f}(1-t)\). For an unnormalised observation \(X+\sigma G\), the signal fraction is \(q=\operatorname{Var}(X)/(\operatorname{Var}(X)+\sigma^2)\) and the signal-to-noise ratio is \(q/(1-q)\). Thus \(q=4/7\) means \(\sigma^2=\tfrac34\operatorname{Var}(X)\), or signal-to-noise ratio \(4/3\).

The Ornstein--Uhlenbeck (Mehler) operator is
\[
T_\rho g(x)=\mathbb E\,g\bigl(\rho x+\sqrt{1-\rho^2}\,G\bigr),\qquad 0\le\rho\le1,
\]
so \(T_0g=\int g\varphi\) and \(T_1g=g\). In semigroup time, \(T_\rho=P_{-\log\rho}\). On the orthonormal probabilists' Hermite basis \(e_j=H_j/\sqrt{j!}\) of \(L^2(\varphi)\), \(T_\rho e_j=\rho^je_j\). The number operator is \(Ne_j=je_j\). Norms without a subscript are \(L^2(\varphi)\) norms; Lebesgue norms are marked \(L^2(dx)\).

A perturbation \(g\) is **admissible** if it is real, even and bounded, with \(\int g\varphi=\int x^2g\varphi=0\). Then \(f=\varphi(1+g)\) is an even density of variance one whenever \(\|g\|_\infty\le1\). For an admissible \(g\), the output density of the mixture is \(\varphi(1+L_g+M_g)\), with
\[
L_g(q,t)=T_{\sqrt{qt}}\,g+T_{\sqrt{q(1-t)}}\,g
\tag{L}
\]
and
\[
M_g(q,t)(x)=\mathbb E\bigl[g(G_1)g(G_2)\,\big|\,\sqrt{qt}\,G_1+\sqrt{q(1-t)}\,G_2+\sqrt{1-q}\,G_3=x\bigr].
\tag{Q}
\]
Here \(G_1,G_2,G_3\) are independent standard Gaussians.

## 2.2 The local-concavity set

In words, \(q\) belongs to \(\mathcal C\) if every symmetric variance-one input whose density is uniformly within a small relative distance of the Gaussian has a concave entropy curve at signal fraction \(q\). Formally,
\[
\mathcal C=
\left\{q\in(0,1]:
\begin{array}{l}
\text{there is }\delta>0\text{ such that }F_{q,\varphi(1+g)}\text{ is concave on }[0,1]\\
\text{for every admissible }g\text{ with }\|g\|_\infty\le\delta
\end{array}\right\}.
\]
The set \(\mathcal C\) is downward closed (Section I.6): Gaussian smoothing maps admissible perturbations to admissible perturbations without increasing \(\|g\|_\infty\), and adding more noise to the input is the same as lowering \(q\).

## 2.3 Results

**Theorem A (local concavity below four sevenths).** For every \(Q\in(0,4/7)\), there is an explicitly computable \(\delta_Q>0\) such that every admissible \(g\) with \(\|g\|_\infty\le\delta_Q\) satisfies
\[
F_{q,\varphi(1+g)}''(t)
\le-\frac1{16}
\left\|N\left(T_{\sqrt{qt}}g+T_{\sqrt{q(1-t)}}g\right)\right\|_{2}^2
\tag{A}
\]
for all \(0<q\le Q\) and \(t\in[0,1]\), with continuous one-sided derivatives at \(t=0,1\). For \(g\ne0\), the right-hand side is strictly negative. The constant \(1/16\) is not optimised, and neither is \(\delta_Q\). The radius is tiny: \(\delta_{0.5}\approx1.7\times10^{-5}\), and \(\delta_Q\sim C_\delta(4/7-Q)^{9/4}\) as \(Q\uparrow4/7\) (Remark R). Theorem A is therefore qualitative. Its positive class does not require log-concavity or smoothness of the input.

**Corollary (critical supremum).**
\[
(0,4/7)\subseteq\mathcal C\subseteq(0,4/7],
\qquad
\boxed{\sup\mathcal C=4/7.}
\tag{B}
\]
Theorem A gives the first inclusion. Luo's symmetric counterexamples [6, Theorem 1.2] give the second: for every \(q>4/7\) they converge to \(\varphi\) in relative \(C_b^0\) norm, are even, and can match the Gaussian mass and variance, so they eventually enter every relative \(L^\infty\) ball in the definition of \(\mathcal C\). Theorem B below gives a second proof of the upper inclusion. Thus \(\mathcal C\) is either \((0,4/7)\) or \((0,4/7]\). This article does not decide which.

The same supremum \(4/7\) holds if the neighbourhoods are defined by any fixed relative \(C_b^m\) norm, or if the inputs are restricted to smooth strongly log-concave densities matching any prescribed finite number of Gaussian moments. The positive result holds in the larger class, and the counterexamples lie in the smaller one. Evenness is essential: Luo's asymmetric counterexamples [6, Theorem 1.1] fail concavity near the endpoints for every \(q>0\).

**Theorem B (local minima above four sevenths).** Fix an integer \(K\ge1\). Luo's moment-corrected Hermite construction [6] gives a single sequence \(f_n\), indexed by large even \(n\), of smooth, strictly positive, even, strongly log-concave densities that match the Gaussian moments through degree \(2K\) and satisfy \(\|f_n/\varphi-1\|_{C_b^m}\to0\) for every fixed \(m\). The smoothed densities have the same properties, uniformly in \(q\in[0,1]\). For every fixed \(q\in(4/7,1]\), eventually
\[
F_{q,f_n}''(1/2)>0,
\tag{C}
\]
so equal weighting is a strict local entropy minimum. More precisely, with \(\epsilon_n=e^{-\sqrt n}\), \(\gamma(q)=21\sqrt5\,q/(49q-4)\) and \(D_{q,n}=h(\varphi)-F_{q,f_n}(1/2)>0\),
\[
\frac{F_{q,f_n}(1/2+s/\sqrt n)-F_{q,f_n}(1/2)}
     {D_{q,n}}
\longrightarrow
1-\exp\!\left[-4\gamma(q)\,s^2\right]
=1-\exp\!\left[-\frac{84\sqrt5\,q}{49q-4}s^2\right]
\tag{D}
\]
in \(C^2\) on compact \(s\)-intervals, locally uniformly in \(q\in(4/7,1]\). The leading term of \(F_{q,f_n}''(1/2)\) is given explicitly in (M1) and (M35). The dip is minute. For example, the leading-order formula gives \(D_{q,n}\approx2\times10^{-13}\) nats at \(q=0.6\), \(n=8\), and it decreases rapidly with \(n\).

**Proposition S (the interpolation range is sharp).** For \(0<r<1\) and \(b(r)<s\le1\), the ratio \(\|T_{\sqrt s}g\|_2^2/(\|g\|_\infty\|T_{\sqrt r}g\|_2)\) is unbounded over bounded \(g\ne0\), even within admissible inputs. Within the Hermite test functions \(h_{n,a}\), the extremal parameter on the boundary \(s=b(r)\) is \(a=(1-r^2)^{-1/2}\). At \((r,s)=(2/7,4/7)\) this is Luo's parameter \(a=7/(3\sqrt5)\). The endpoint \(s=b(r)\) itself is not decided.

**Proposition E (Luo's family at the boundary).** At \(q=4/7\), Luo's moment-corrected sequence with any amplitudes \(\epsilon_n\to0\) satisfying \(\epsilon_nn^{1/4}\to0\) (including \(\epsilon_n=e^{-\sqrt n}\)) has \(F_{4/7,f_n}''<0\) on a fixed neighbourhood of \(t=1/2\) for all large \(n\). This does not decide whether \(4/7\in\mathcal C\).

## 2.4 Summary of the regimes

| Signal fraction \(q\) | Near-Gaussian behaviour in uniform relative-density neighbourhoods | Source |
|---|---|---|
| \(0<q<1/2\) | Concave in a neighbourhood | Schubr [7]; also Theorem A |
| \(1/2\le q<4/7\) | Concave in a neighbourhood | Theorem A (new) |
| \(q=4/7\) | Open. Luo's family is strictly concave near the midpoint | Proposition E |
| \(4/7<q\le2-\sqrt2\) | Not concave; strict midpoint local minima | Luo [6] for non-concavity; Theorem B for the midpoint minima (also proposed in [8]) |
| \(2-\sqrt2<q\le1\) | Not concave; strict midpoint local minima | Luo [6, Appendix E]; Theorem B adds \(\gamma(q)\) and the profile (D) |

These are local theorems. They do not settle concavity for all symmetric log-concave inputs at any positive noise level (Luo's Problem 10.1). In the language of Luo's Problem 10.2 [6], Theorem A answers the question negatively for every \(q<4/7\), under all the extra regularity and moment conditions posed there; the instance \(q=4/7\) remains open. Theorem B gives a complete proof of an affirmative answer to Problem 10.3, which was proposed earlier in [8].



# Part I. Proof of Theorem A (local concavity below four sevenths)

**Roadmap.** The second derivative of the entropy curve is governed by the output perturbation \(\mathcal R=L+M\) of (L) and (Q). The linear part \(L\) gives a strictly negative contribution of size \(D^2=\|NL\|_2^2\) (inequality (P18)). The proof shows that everything else is at most a small multiple of \(D^2\). It has five steps.

1. Section I.1 recalls a comparison principle for Gaussian kernels.
2. Section I.2 uses it to prove Lemma I.2: for bounded inputs, a smoothed norm at level \(s\) is controlled by the norm at level \(r\) whenever \(s<b(r)\). This is what makes the quadratic term \(M\) small compared with \(L\).
3. Section I.3 shows that, after a little extra smoothing, \(M\) can be differentiated twice in the weight \(t\), with explicit bounds.
4. Section I.4 combines the two lemmas at a fixed \(q<4/7\). This is the only place where \(q<4/7\) is used.
5. Sections I.5 and I.6 absorb the logarithmic nonlinearity, following Schubr [7, Section 2], fix an explicit radius, and make it uniform in \(q\le Q\). Section I.7 evaluates the radius numerically.

## I.1 A Gaussian-kernel comparison

Let
\[
K_j(x,y)=c_j\exp[-\alpha_j(x^2+y^2)+\beta_jxy],
\quad c_j>0,\quad 2\alpha_j>\beta_j>0,\qquad j=A,B.
\tag{P1}
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

This is the Gaussian-kernel lemma of Luo [6, Appendix F, (F.4)--(F.5)], which rests on weighted composition operators on Fock space [10, 11]. We include a proof sketch for completeness.

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

## I.2 Weighted interpolation for bounded functions

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

### Lemma I.2 (weighted interpolation)

Luo [6, Theorem F.1 and Corollary F.2] proves the following inequality with \(\|g\|_{L^2(dx)}\) in place of \(\|g\|_\infty\), including the endpoint \(s=b(r)\). The version below admits bounded inputs, whose Lebesgue norm may be infinite. It needs the strict inequality in (P5). Proposition S in Part III shows that the range (P5) cannot be enlarged.

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

The kernel-comparison condition (P2) between \(A_s^\varepsilon\) and \(B_\varepsilon\) is the pair of non-strict inequalities corresponding to the following display. We prove the strict versions, which suffice:
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

## I.3 Bilinear derivatives after smoothing

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
The same estimate holds on the larger disk \(|z-t|\le\tfrac32\eta_\rho\), where \(P\le\rho(1+3\eta_\rho)=(3+\rho)/4<1\); this is the smoothing margin. Absolute-coefficient majorisation and the exponentially decaying output multipliers give normal convergence, and likewise for every fixed number-operator derivative. Hence the map is Hilbert-valued holomorphic there. Cauchy's derivative estimate proves (P13).

Finite Hermite expansions are dense in \(L^2\), and the bilinear contraction estimate, applied to differences, proves convergence to the actual conditional expectation. The same bounds prove \(C^2\) convergence in the weight, so the assertions include both real endpoints.

## I.4 Apply the two lemmas at a fixed signal

Fix \(q\in(0,4/7)\). Put
\[
r=q/2,\qquad b(r)=\sqrt{2r(1+r)}-r,\qquad
s=\frac{q+b(r)}2.
\]
The inequality \(b(r)>q\) is equivalent to \(q<4/7\): squaring \(\sqrt{q(1+q/2)}>3q/2\) gives \(q(7q-4)<0\). Therefore
\[
r<q<s<b(r)<1,
\]
where \(b(r)<1\) because \(2r(1+r)<(1+r)^2\).
Let \(J=J_{r,s}\) be (P10), \(\rho=q/s\), and write
\(\eta=\eta_\rho,\ R_\rho,\ K=K_\rho\).

The exact output density has relative perturbation \(L_g(q,t)+M_g(q,t)\), as in (L) and (Q), and \(M_g(q,t)=\mathcal B_q(t;g,g)\). This follows by expanding the product of the two input density ratios in the conditional-expectation representation. For the original perturbation \(g\), set \(h=T_{\sqrt s}g\). The Gaussian/Hermite composition identity (in the Hermite coordinates of [12]) gives exactly
\[
M_g(q,t)=\mathcal B_\rho(t;h,h).                         \tag{P15}
\]
Lemma I.2 gives
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

## I.5 Nonlinear curvature and an explicit radius

The steps in this section follow Schubr's entropy-energy argument [7, Section 2]. We reproduce them so that the proof is self-contained. The output ratio is \(1+\mathcal R\), where \(\mathcal R=L+M\), and
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

Put \(\psi(u)=\log(1+u)-u\). For \(N=-\partial_{xx}+x\partial_x\), the chain rule gives
\[
N\psi(\mathcal R)=\psi'(\mathcal R)\,N\mathcal R-\psi''(\mathcal R)\,(\mathcal R_x)^2 .
\]
For \(\delta\le1/12\), \(|\mathcal R|\le3\delta\le1/4\), so \(|\psi'(\mathcal R)|=|\mathcal R|/(1+\mathcal R)\le4\delta\) and \(|\psi''(\mathcal R)|=(1+\mathcal R)^{-2}\le16/9\). With (P19) and \(\|\mathcal R\|_\infty\le3\delta\), this gives \(\|N\psi(\mathcal R)\|_2\le(4+16)\delta\|N\mathcal R\|_2\). We use the weaker constant \(30\):
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
The pairing with \(N^{-1}\) below is legitimate because \(L_{tt}\) has no Hermite component of degree zero. Drop the last nonpositive term and expand
\(\log(1+\mathcal R)=L+M+\psi(\mathcal R)\).
The error terms satisfy
\[
|\langle L_{tt},M\rangle|\le(KJ/8)\delta D^2,\qquad
|\langle L_{tt},\psi(\mathcal R)\rangle|=|\langle N^{-1}L_{tt},N\psi(\mathcal R)\rangle|\le60\delta D^2,
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

The Hilbert-valued regularity established in Section I.3 and the corresponding elementary linear regularity imply \(\mathcal R\in C^2([0,1];L^2(\varphi))\). The actual conditional-expectation density ratio, rather than a Hermite truncation, supplies the uniform bounds \(1-3\delta\le1+\mathcal R\le1+3\delta\).

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

## I.6 One radius for all smaller signal fractions

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

## I.7 Numerical values of the radius

The following values are 100-decimal diagnostics of the closed-form scalar recipe (P7), rounded here for display. They are not interval certificates and are not used in the proof. The checker also verifies the direct Gaussian square root against the independent Mehler parametrisation and direct kernel composition:

| \(q\) | chosen \(s\) | chosen \(\varepsilon\) | strict kernel gap | \(J\) | sufficient radius |
|---:|---:|---:|---:|---:|---:|
| 0.50 | 0.5202847075 | 0.1375827353 | 0.0473458334 | 1.4320544662 | \(1.7134\cdot10^{-5}\) |
| 0.53 | 0.5419050561 | 0.0870079716 | 0.0296201831 | 1.5902847695 | \(4.8520\cdot10^{-6}\) |
| 0.56 | 0.5633202098 | 0.0262722883 | 0.0088148014 | 2.1257516652 | \(2.5717\cdot10^{-7}\) |
| 0.57 | 0.5704164638 | 0.0033875871 | 0.0011301518 | 3.5370963892 | \(2.3569\cdot10^{-9}\) |

The shrinking radius near \(4/7\) is expected: both the interpolation margin and the extra smoothing used for derivatives vanish at that endpoint. The radius is sufficient and deliberately unoptimised. Remark R in Section III.3 shows that \(\delta_Q\sim C_\delta(4/7-Q)^{9/4}\), with exponent \(2\) from the derivative smoothing and \(1/4\) from the weighted interpolation. The checker `checks/local_positive_checks.py` reproduces every entry of this table to the displayed digits.



# Part II. Proof of Theorem B (midpoint minima above four sevenths)

**Roadmap.** We use Luo's explicit Hermite construction. Its quadratic interaction term \(M\) is a four-variable generating-function coefficient. Section II.2 records the construction and two exponential rates whose comparison (M16) is strict exactly when \(q>4/7\). Section II.3 proves a saddle-point lemma that controls two derivatives of such coefficients. Section II.4 evaluates the curvature constant \(\gamma(q)\). Section II.5 bounds the derivatives of the density itself, which is what the entropy needs. Section II.6 transfers the result to entropy, and Section II.7 records checks and limits.

## II.1 Detailed statement

We restate Theorem B in the form used in the proof, with the notation of Section 2. If \(X\) has density \(f\), define \(\mathcal S_qf\) as the density of
\[
\sqrt q\,X+\sqrt{1-q}\,G,
\]
where \(G\) is an independent standard Gaussian. The variable \(q\) is the signal fraction. For two independent variables of density \(\mathcal S_qf\), put
\[
F_{q,f}(t)=h\bigl(\sqrt t\,X_q+\sqrt{1-t}\,Y_q\bigr),
\quad 0\le t\le1.
\]
The variance remains one, and exchange of the copies gives \(F_{q,f}(t)=F_{q,f}(1-t)\).

Luo [6] already constructed symmetric counterexamples for every \(q>4/7\), and proved positive midpoint curvature for \(q>2-\sqrt2\) [6, Appendix E]. His Problem 10.3 asks for positive midpoint curvature on the remaining interval \((4/7,2-\sqrt2]\). An earlier proposed argument for that interval, using the same construction and a differentiated saddle-point mechanism, was posted publicly on 1 November 2026 [8]; it is unrefereed. It bounds the curvature from below but does not give the exact constant, and it excludes \(q=1\). The statement below gives the exact constant \(\gamma(q)\), the local profile and a single sequence for all \(q\in(4/7,1]\). It is not a first refutation of the original conjecture.

**Theorem B (midpoint curvature and local profile).** Fix an integer \(K\ge1\). There is a sequence \(f_n\), indexed by sufficiently large even integers, of strictly positive, smooth, even, variance-one densities with the following properties.

1. Every \(f_n\) matches the standard Gaussian moments through degree \(2K\).
2. Both \(f_n\) and \(\mathcal S_qf_n\) are strongly log-concave for sufficiently large \(n\), with a threshold for \(n\) that is uniform in \(q\in[0,1]\).
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

Items 4 and 5 hold uniformly for \(q\) in compact subintervals of \((4/7,1]\); this is a different uniformity from items 2 and 3. The same sequence of densities works for all these fixed \(q\). No claim at \(q=4/7\) follows.

Equation (M3) specifies the local geometry: the entropy dip has width of order \(n^{-1/2}\). On this scale its curvature changes sign at \(|s|=1/\sqrt{8\gamma(q)}\) in the limiting profile. This is an asymptotic statement about the constructed sequence, not a universal description of all symmetric inputs.

## II.2 The inherited construction and its estimates

We record the construction from [6] explicitly so that the new argument has a fixed referent. The probabilists' Hermite polynomials satisfy
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
This identity follows by applying Mehler's identity to the two correlated Gaussian arguments in the squared norm (a complete proof for every \(a>0\) is Lemma S.1 in Part III). For \(r>2/7\), \(0<\nu(r)<\tau(r)\le1\). At \(r=2/7\) the two rates coincide; for \(r<2/7\) the order is reversed. Only \(r>2/7\) is used in Part II. The binomial coefficients of \((1-z)^{-1/2}\) are positive and their convolution sums to one; hence
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

## II.3 A differentiated coefficient lemma

The following elementary saddle calculation is the main additional analytic step. The method is standard; the derivative control and its application to this coefficient are stated in full.

**Lemma II.3 (positive quadratic coefficient).** Let \(B(u)\) be a smooth real symmetric \(4\times4\) matrix, with every entry strictly positive near \(u=0\). Set
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

## II.4 The coefficient for the quadratic convolution

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
\(i^{2n}(-i)^{2n}=1\). Equations (M20)--(M21) and their multigraph interpretation are already present in [6]. The shipped checker `checks/identity_checks.py` verifies them independently against a direct Hermite expansion. We use them to evaluate the derivative.

**Lemma II.4 (explicit saddle and curvature).** Put \(U=\kappa+b\). For the matrix (M21),
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
The coefficient is even in \(u\), so its first derivative at zero vanishes. Lemma II.3 proves (M22). The eigenvalues of \(B(0)\) are
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

## II.5 Derivatives of the density on the correct scale

A scalar norm asymptotic alone does not control the entropy remainder. We need bounds for \(M_t'\) and \(M_t''\) in the same Gaussian \(L^2\) space.

**Lemma II.5 (two-weight kernel bound).** On a sufficiently small fixed neighbourhood \(J\) of \(1/2\), and for \(j=0,1,2\),
\[
\|\partial_t^jM_t\|_2\le C(1+n)^j\|M_t\|_2.
\tag{M27}
\]

*Proof.* For fixed \(n\), \(t\mapsto M_t\) is \(C^\infty\) on \(J\) as an \(L^2(\varphi)\)-valued map: in the conditional representation displayed before (M29), the integrand is a Gaussian factor times a polynomial with coefficients smooth in \(t\), so one may differentiate under the expectation. Hence \(\partial_t^j\partial_s^j\langle M_t,M_s\rangle|_{s=t}=\|\partial_t^jM_t\|_2^2\). Define \(G_n(t,s)=\langle M_t,M_s\rangle\). Its four-variable Gaussian covariance has off-diagonal block \(qv_tv_s^{\mathsf T}\); its eigenvalues are \(1+q,1-q,1,1\). Thus the determinant prefactor in its generating function remains \(\Delta\), and its phase-changed coefficient matrix has diagonal blocks
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
where \(Z,W\) are independent standard Gaussians. On \(J\), the first two derivatives of \(v_t\) are bounded. Differentiating \(\widetilde h_n\) or \(h_n\) evaluated at the represented coordinates produces spatial derivatives multiplied by the bounded vectors \(\partial_tv_t,\partial_t^2v_t\) and by polynomials in \(x,Z,W\). By (M7), each spatial derivative of \(h_{n,a}\) is a combination of neighbouring \(h_{n\pm i,a}\) with coefficients \(O(n^{i/2})\). Apply this to the difference \(\widetilde h_n-h_n\), whose size is (M8), and integrate the resulting fixed Gaussian moments. For some fixed exponent \(A\),
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

## II.6 Entropy transfer and proof of Theorem B

Write the output relative perturbation as
\[
R_t=\epsilon_n\widetilde L_t+\epsilon_n^2\widetilde M_t.
\]
Conditional expectation and the supremum bounds in Section II.2 give, uniformly in \(t\),
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

Lemma II.5 and (M26) also give
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
For the present densities, whose relative perturbations are Gaussian factors times polynomials (not finite Hermite expansions), differentiation is justified by these Gaussian envelopes and the uniform positive lower bound on \(1+R\).

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
The entropy dip therefore becomes extremely small as \(n\) grows. For example, the leading-order value of \(D_{q,n}\) is about \(2\times10^{-13}\) nats at \(q=0.6\), \(n=8\), and about \(2\times10^{-7}\) nats at \(q=1\), \(n=8\). The theorem concerns the sign and local shape of the dip, not a large practical entropy loss.

## II.7 Verification and limits

The proof is analytic and has not been checked by a proof assistant. A separate implementation evaluates the coefficient and its first two weight derivatives by a Gaussian-moment recurrence. Exact rational tests also enumerate all admissible four-vertex multigraph configurations independently.

Thirty numerical cases in `checks/coefficient_checks.py` use \(q=4/7,0.575,0.6,0.75,1\) and \(n=2,4,8,12,16,24\). For \(q>4/7\) the checker requires \(|C_n''/(nC_n)+8\gamma(q)|\le2/n\) and a matching tolerance for the prefactor in (M23). These cases support the coefficient asymptotic; they do not establish the unrestricted theorem. At \(n=24\), for example, \(C_n''/(nC_n)\) is \(-8.80988058680\) for \(q=0.6\), against the limiting value \(-8.87384457213\).

The \(q=4/7\) computations concern only the generating coefficient. They do not determine the entropy-curvature sign at that signal fraction: the strict linear/quadratic comparison in (M16) is absent there.

A useful exact check is \(b=0\). Then diagonal rescaling gives
\[
C_n(u)=(1-4u^2)^nC_n(0),
\]
so \(C_n''(0)/C_n(0)=-8n\) exactly, agreeing with (M22).

The construction and the exponential comparison in (M16) come from Luo [6]. The sign of the midpoint curvature on \((4/7,2-\sqrt2]\) was also the subject of the earlier public proposal [8]. The additions here are the exact differentiated coefficient asymptotic (M22), the explicit prefactor (M35), the derivative-scale estimates of Section II.5, the entropy profile (M3), and a single argument that also covers \(q=1\) (excluded in [8]).

The theorem does not settle the uniform concavity question for all symmetric log-concave inputs. It does not settle the boundary \(q=4/7\). Positive local results of Schubr [7] and Abu-Ajamieh [9] concern different input classes or signal ranges and must be distinguished from the counterexample construction here.


# Part III. Why the method stops at four sevenths

**Roadmap.** Section III.1 shows that the interpolation range of Lemma I.2 cannot be enlarged, and that at the critical point the extremal member of the Hermite test family is Luo's counterexample family (Proposition S). In words: the bounded test functions \(h_{n,a}\) can break the inequality exactly when \(s>b(r)\), and on the borderline only one parameter \(a\) comes close. Section III.2 shows that this same family, at \(q=4/7\) exactly, is strictly concave near equal weighting (Proposition E). Section III.3 records how fast the radius of Theorem A shrinks near \(4/7\) (Remark R). Equation labels use S, E and R. Notation is that of Section 2: \(T_\rho\) is the Ornstein--Uhlenbeck operator, \(e_k=H_k/\sqrt{k!}\) is the orthonormal probabilists' Hermite basis of \(L^2(\varphi)\), and all norms without subscript are \(L^2(\varphi)\) norms. Throughout,
\[
h_{n,a}(x)=\sqrt2\,e^{-ax^2/2}\frac{H_n(\sqrt{2a}\,x)}{\sqrt{n!}},\qquad a>0,
\]
\[
\tau_a(u)=\frac{1-a+au}{1+a-au},\qquad
\nu_a(u)=\frac{a-1+au}{1+a+au},\qquad
\varrho_a(u)=\max\{|\tau_a(u)|,|\nu_a(u)|\},
\]
\[
\Pi_a(u)=(1+a)^2-a^2u^2,\qquad
c_j=\binom{2j}{j}4^{-j},\qquad
P_{n,a}(u)=\sum_{j+l=n}c_jc_l\,\tau_a(u)^j\nu_a(u)^l .
\]
For \(0<r<1\) put
\[
b(r)=\sqrt{2r(1+r)}-r,\qquad
a_r=(1-r^2)^{-1/2},\qquad
\lambda_r=\frac{r}{1+\sqrt{1-r^2}}=\tanh\!\big(\tfrac12\operatorname{artanh} r\big).
\tag{S0}
\]

## III.1 Sharpness of Lemma I.2 for bounded inputs

### III.1.1 Three elementary facts about \(h_{n,a}\)

**Lemma S.1 (squared smoothing norms).** For every \(a>0\), \(n\ge0\) and \(u\in[0,1]\),
\[
\|T_{\sqrt u}h_{n,a}\|_2^2=\frac{2}{\sqrt{\Pi_a(u)}}\,P_{n,a}(u).
\tag{S1}
\]
Consequently \(\|T_{\sqrt u}h_{n,a}\|_2^2\le2\Pi_a(u)^{-1/2}\varrho_a(u)^n\), and if \(\tau_a(u),\nu_a(u)\ge0\) then
\[
\|T_{\sqrt u}h_{n,a}\|_2^2\ge\frac{2c_n}{\sqrt{\Pi_a(u)}}\max\{\tau_a(u),\nu_a(u)\}^n
\ge\frac{\varrho_a(u)^n}{\sqrt{n\,\Pi_a(u)}}\quad(n\ge1).
\tag{S2}
\]

*Proof.* Since \(\langle T_{\sqrt u}f,T_{\sqrt u}g\rangle=\langle f,T_ug\rangle\), the left side of (S1) is \(\mathbb E[h_{n,a}(X)h_{n,a}(Y)]\), where \((X,Y)\) is centred Gaussian with covariance \(\Sigma=\begin{pmatrix}1&u\\u&1\end{pmatrix}\) (for \(u=1\), \(X=Y\)). The Hermite generating function gives
\(\sum_n h_{n,a}(x)z^n/\sqrt{n!}=\sqrt2\,e^{-ax^2/2+\sqrt{2a}xz-z^2/2}\). By (S6) below each \(|h_{j,a}|\) is bounded by a polynomial in \(j\), so the double series converges absolutely and uniformly in \((x,y)\) for complex \(z,w\), and
\[
\sum_{j,k}\mathbb E[h_{j,a}(X)h_{k,a}(Y)]\frac{z^jw^k}{\sqrt{j!k!}}
=2\,\mathbb E\exp\Big(-\tfrac a2(X^2+Y^2)+\sqrt{2a}(zX+wY)\Big)e^{-(z^2+w^2)/2}.
\]
For positive semidefinite \(\Sigma\) and \(\ell\in\mathbb C^2\),
\(\mathbb E\exp(-\frac a2|V|^2+\ell^{\mathsf T}V)=\det(I+a\Sigma)^{-1/2}\exp(\frac12\ell^{\mathsf T}\Sigma(I+a\Sigma)^{-1}\ell)\).
Here \(\det(I+a\Sigma)=\Pi_a(u)\), and \(\Sigma(I+a\Sigma)^{-1}\) has eigenvalues \((1\pm u)/(1+a(1\pm u))\) on \((1,\pm1)/\sqrt2\). With \(y=a(1\pm u)\) and \(\frac{y}{1+y}-\frac12=\frac{y-1}{2(1+y)}\), the exponent becomes
\[
\tfrac14\big[\nu_a(u)(z+w)^2-\tau_a(u)(z-w)^2\big].
\]
Hence \(\|T_{\sqrt u}h_{n,a}\|^2=2\Pi_a(u)^{-1/2}\,n!\,[z^nw^n]F\), with \(F\) the exponential of this quadratic form. Let \(W=X'+iY'\) be a standard complex Gaussian (\(X',Y'\) independent \(N(0,1/2)\)), so that \(\mathbb E[W^j\overline W^k]=k!\,\delta_{jk}\). For \(|\zeta|\) small, absolute convergence (the majorant \(\exp(C|\zeta||W|^2)\) is integrable) permits termwise expectation:
\[
\sum_n\zeta^n\,n!\,[z^nw^n]F=\mathbb E\,F(\sqrt\zeta W,\sqrt\zeta\,\overline W)
=\mathbb E\,e^{\zeta(\nu X'^2+\tau Y'^2)}
=(1-\zeta\nu)^{-1/2}(1-\zeta\tau)^{-1/2},
\]
using \((W+\overline W)^2=4X'^2\) and \((W-\overline W)^2=-4Y'^2\). Expanding \((1-x)^{-1/2}=\sum c_jx^j\) and comparing coefficients gives (S1). This is (M11) of Part II, now for every \(a>0\) and \(u\in[0,1]\).

Since \(c_j>0\) and \(\sum_{j+l=n}c_jc_l=1\) (the square of \((1-x)^{-1/2}\) is \((1-x)^{-1}\)), \(|P_{n,a}(u)|\le\varrho_a(u)^n\). If both rates are nonnegative, all terms are nonnegative and \(P_{n,a}\) is at least the single term \(c_n\max\{\tau,\nu\}^n\). Finally \(c_n^2=\frac14\prod_{i=2}^n\big(\frac{2i-1}{2i}\big)^2\ge\frac14\prod_{i=2}^n\frac{i-1}{i}=\frac1{4n}\), because \((2i-1)^2\ge 2i(2i-2)\). \(\square\)

We shall also use \(c_m\le(2m+1)^{-1/2}\), which follows from \(((2i-1)/(2i))^2\le(2i-1)/(2i+1)\) and telescoping.

**Lemma S.2 (hyperbolic form of the rate).** For \(a>0\) and \(0\le u<1\),
\[
\varrho_a(u)=\tanh\Theta_a(u),\qquad
\Theta_a(u)=\tfrac12\Big[\operatorname{artanh}u+\big|\log a-\log a_u\big|\Big].
\tag{S3}
\]
Moreover, for \(u>0\), \(\varrho_a(u)=\tau_a(u)>0\) if \(a\le a_u\), and \(\varrho_a(u)=\nu_a(u)>0\) if \(a\ge a_u\). At \(u=0\) the rates are nonnegative, and they vanish only at \(a=1\). In particular \(\min_a\varrho_a(u)=\varrho_{a_u}(u)=\tau_{a_u}(u)=\nu_{a_u}(u)=\lambda_u\).

*Proof.* Put \(m(x)=(1-x)/(1+x)\), so that \(|m(x)|=\tanh(\frac12|\log x|)\), \(\tau_a(u)=m(a(1-u))\) and \(\nu_a(u)=-m(a(1+u))\). With \(A=\log a\), the two numbers \(|A+\log(1-u)|\) and \(|A+\log(1+u)|\) are the distances from \(A\) to \(-\log(1-u)\) and \(-\log(1+u)\). These points have midpoint \(-\frac12\log(1-u^2)=\log a_u\) and half-distance \(\frac12\log\frac{1+u}{1-u}=\operatorname{artanh}u\). The larger distance is therefore \(\operatorname{artanh}u+|A-\log a_u|\), and \(\tanh\) is increasing. The sign statements follow because \(a(1-u)<1\) when \(a\le a_u\) and \(a(1+u)>1\) when \(a\ge a_u\) (note \(a_u(1-u)<1<a_u(1+u)\)). At \(a=a_u\), \(\Theta=\frac12\operatorname{artanh}u\), and \(\tanh(\frac12\operatorname{artanh}u)=u/(1+\sqrt{1-u^2})\). \(\square\)

**Lemma S.3 (supremum bound).** For every \(a>0\) and \(n\ge0\),
\[
\|h_{n,a}\|_\infty\le B_n=\sqrt{2e}\,(n+1)^{1/4},
\qquad
\|h_{n,a}\|_\infty\ge\|h_{n,a}\|_2\ge\Big(\frac{2c_n}{\sqrt{1+2a}}\Big)^{1/2}\quad(a\ge\tfrac12).
\tag{S6}
\]

*Proof.* The substitution \(y=\sqrt{2a}\,x\) shows that \(\|h_{n,a}\|_\infty=\sqrt2\sup_ye^{-y^2/4}|H_n(y)|/\sqrt{n!}\), independent of \(a\). Mehler's formula on the diagonal, \(\sum_m\rho^mH_m(y)^2/m!=(1-\rho^2)^{-1/2}\exp(\rho y^2/(1+\rho))\) for \(0<\rho<1\), keeps only the \(m=n\) term to give
\[
e^{-y^2/2}\frac{H_n(y)^2}{n!}\le\rho^{-n}(1-\rho^2)^{-1/2}e^{-(1-\rho)y^2/(2(1+\rho))}\le\rho^{-n}(1-\rho^2)^{-1/2}.
\]
For \(n\ge1\) take \(\rho=n/(n+1)\): then \(\rho^{-n}\le e\) and \(1-\rho^2\ge1/(n+1)\). Hence \(h_{n,a}^2\le2e(n+1)^{1/2}\). For \(n=0\), \(|h_{0,a}|\le\sqrt2\). The lower bound uses \(\varphi\,dx\) as a probability measure and (S2) at \(u=1\), where \(\tau_a(1)=1\), \(\nu_a(1)=(2a-1)/(2a+1)\ge0\) and \(\Pi_a(1)=1+2a\). \(\square\)

*Remark.* The \(n^{1/4}\) in \(B_n\) is an artefact of this elementary proof. A classical uniform bound for Hermite functions (Cramér's inequality) removes the factor \((n+1)^{1/4}\). We do not use it in any proof below. Every polynomial factor is harmless for the exponential divergence in Proposition S. It matters only for the amplitude condition in Proposition E; see Remark E.6.

### III.1.2 The exact envelope of the rate condition

**Proposition S.4 (envelope).** Let \(0<r<s\le1\).

1. If \(s\le b(r)\), then \(\varrho_a(s)^2\le\varrho_a(r)\) for every \(a>0\). Equality holds for some \(a\) if and only if \(s=b(r)\), and then only at \(a=a_r\).
2. If \(s>b(r)\), then \(\tau_{a_r}(s)^2>\lambda_r=\varrho_{a_r}(r)\).
3. In closed form, with \(w=\sqrt{\lambda_r}\),
\[
\boxed{\;
s>b(r)\iff\frac{s-r}{2-r-s}>\frac{w}{1+w+w^2},
\qquad
b(r)=\frac{2w(1-w+w^2)}{1+w^4}.
\;}
\tag{S4}
\]

Thus the set of \(s\) for which some member of the family \(\{h_{n,a}\}\) has an exponentially favourable rate is exactly \((b(r),1]\). The extremal parameter on the envelope \(s=b(r)\) is \(a=a_r=(1-r^2)^{-1/2}\), which is Luo's choice in [6, Appendix F].

*Proof.* We use one scalar fact. For \(\Theta>0\) and \(\delta>0\), write \(p=\tanh\Theta=w^2\) with \(0<w<1\), and \(T=\tanh\delta\). Then
\[
\tanh(\Theta+\delta)-w=\frac{w^2+T}{1+w^2T}-w
=\frac{(1-w)\,[\,T(1+w+w^2)-w\,]}{1+w^2T}.
\]
Therefore \(\tanh^2(\Theta+\delta)>,=,<\tanh\Theta\) according as \(T>,=,<\phi(w)\), where \(\phi(w)=w/(1+w+w^2)\). The function \(\phi\) is strictly increasing on \((0,1)\), since \(\phi'(w)=(1-w^2)/(1+w+w^2)^2\).

First let \(s<1\), and set \(A=\log a\), \(A_u=\log a_u\). Since \(a_u\) is increasing in \(u\), \(A_r<A_s\). By (S3), \(\Theta_a(u)\) is piecewise linear in \(A\) with slopes \(\pm\frac12\) and a kink at \(A_u\). We consider three regions.

*Region I, \(A\le A_r\).* Here \(\Theta_a(u)=\frac12(-\log(1-u)-A)\) for \(u=r,s\), because \(\operatorname{artanh}u+\log a_u=-\log(1-u)\). Hence \(\Theta_a(s)=\Theta_a(r)+\delta_I\) with
\(\tanh\delta_I=(s-r)/(2-r-s)\), independent of \(a\). Also \(\Theta_a(r)\ge\Theta_{a_r}(r)\), with equality only at \(a=a_r\). Hence \(w(a)=\sqrt{\varrho_a(r)}\ge\sqrt{\lambda_r}\), and so \(\phi(w(a))\ge\phi(\sqrt{\lambda_r})\), with equality only at \(a_r\).

*Region II, \(A_r\le A\le A_s\).* Here \(\Theta_a(r)\) increases and \(\Theta_a(s)\) decreases with \(A\). Hence \(\varrho_a(s)^2/\varrho_a(r)\) is strictly decreasing on this interval.

*Region III, \(A\ge A_s\).* Here \(\Theta_a(u)=\frac12(\log(1+u)+A)\), so \(\Theta_a(s)=\Theta_a(r)+\delta_{III}\) with \(\tanh\delta_{III}=(s-r)/(2+r+s)\), again independent of \(a\). The function \(w(a)\) is increasing in \(A\).

Now suppose \(s\le b(r)\). Since \(s\mapsto(s-r)/(2-r-s)\) is increasing, the defining equation of the envelope (proved below) gives \(\tanh\delta_I\le\phi(\sqrt{\lambda_r})\), with equality iff \(s=b(r)\). In Region I the scalar fact gives \(\varrho_a(s)^2\le\varrho_a(r)\), with equality only if \(a=a_r\) and \(s=b(r)\). On Region II the ratio is strictly below its value at \(a_r\), which is at most one. In particular \(\varrho_{a_s}(s)^2<\varrho_{a_s}(r)\), that is \(\tanh\delta_{III}<\phi(w(a_s))\). Since \(w\) increases on Region III, \(\tanh\delta_{III}<\phi(w(a))\) there as well, so the ratio is strictly less than one. This proves item 1. If \(s>b(r)\), the same equivalence at \(a=a_r\) gives item 2.

It remains to solve \((s-r)/(2-r-s)=\phi(w)\), \(w=\sqrt{\lambda_r}\). Its solution is \(s=[r(1-\phi)+2\phi]/(1+\phi)\). Using \(1-\phi=(1+w^2)/(1+w+w^2)\), \(1+\phi=(1+w)^2/(1+w+w^2)\), and \(r=2w^2/(1+w^4)\) (the double-angle formula for \(\tanh\)), one gets
\[
s=\frac{r(1+w^2)+2w}{(1+w)^2}
=\frac{2w(1+w)(1+w^3)}{(1+w)^2(1+w^4)}
=\frac{2w(1-w+w^2)}{1+w^4}.
\]
On the other hand \(2r(1+r)=4w^2(1+w^2)^2/(1+w^4)^2\), so
\(b(r)=[2w(1+w^2)-2w^2]/(1+w^4)\), which is the same number. This proves (S4).

For \(s=1\), \(\tau_a(1)=1\), so \(\varrho_a(1)=1>\varrho_a(r)\) for every \(a\). Also \(b(r)<1\) for \(r<1\), since \(2r+2r^2<(1+r)^2\). Thus \(s=1\) belongs to case 2. \(\square\)

*Remark (Luo's direct form).* Item 2 can also be written as \(\tau_{a_r}(s)>\sqrt{\lambda_r}\), which is the form in which Luo states the necessity condition in [6, Appendix F]. As \(a\to0\) or \(a\to\infty\), \(\varrho_a(s)^2/\varrho_a(r)\to1\). Thus the supremum over \(a\) of the ratio is always at least one. Item 1 asserts that, for \(s\le b(r)\), the value one is not exceeded and is attained only at the single point stated.

### III.1.3 Divergence of the interpolation quotient

**Proposition S (sharpness of Lemma I.2).** Let \(0<r<1\) and \(b(r)<s\le1\). Then
\[
\sup_{g\ \text{real, bounded},\ g\ne0}
\frac{\|T_{\sqrt s}g\|_2^2}{\|g\|_\infty\,\|T_{\sqrt r}g\|_2}=+\infty .
\tag{S5}
\]
More precisely, with \(g_n=h_{n,a_r}\),
\[
\frac{\|T_{\sqrt s}g_n\|_2^2}{\|g_n\|_\infty\|T_{\sqrt r}g_n\|_2}
\ \ge\ c_{r,s}\,n^{-3/4}\left(\frac{\tau_{a_r}(s)}{\sqrt{\lambda_r}}\right)^{n}
\longrightarrow\infty,
\tag{S7}
\]
where \(c_{r,s}>0\) and \(\tau_{a_r}(s)/\sqrt{\lambda_r}>1\) by Proposition S.4. The same divergence holds along even functions \(\widetilde g_n\), \(n\) even, satisfying \(\int x^{2\ell}\widetilde g_n\varphi=0\) for \(0\le\ell\le K\), for any fixed \(K\ge1\). In particular it holds inside the admissible class of Theorem A.

*Proof.* Write \(a=a_r\). Since \(a\ge1\), both \(\nu_a(s)\) and \(\tau_a(s)\) are positive: \(a(1+s)>1\), and \(\tau_a(s)>\tau_a(r)=\lambda_r>0\) because \(\tau_a\) increases in \(u\). By (S2), \(\|T_{\sqrt s}g_n\|^2\ge\Pi_a(s)^{-1/2}n^{-1/2}\tau_a(s)^n\). By (S1) at \(u=r\), where \(\tau_a=\nu_a=\lambda_r\),
\[
\|T_{\sqrt r}g_n\|_2^2=\frac{2}{\sqrt{\Pi_a(r)}}\lambda_r^{\,n}\sum_{j+l=n}c_jc_l=\frac{2}{\sqrt{\Pi_a(r)}}\lambda_r^{\,n}
\tag{S8}
\]
exactly. Together with \(\|g_n\|_\infty\le B_n\), these give (S7).

For the admissible version, let \(\psi_j(x)=\sqrt2e^{-x^2/2}H_{2j}(\sqrt2x)/\sqrt{(2j)!}\) and define \(\widetilde g_n=g_n-\sum_{j\le K}c_{n,j}\psi_j\) by the triangular system (M5). Integration of the generating function gives
\[
\int\varphi(x)x^{2\ell}e^{-ax^2/2+\sqrt{2a}xz-z^2/2}\,dx=(1+a)^{-1/2}e^{b_az^2/2}Q_\ell(z),
\qquad b_a=\frac{a-1}{a+1}.
\]
Here \(Q_\ell\) is a polynomial of degree \(2\ell\) whose coefficients depend only on \(a\) and \(\ell\): they are moments of a Gaussian of mean \(\sqrt{2a}z/(1+a)\) and variance \(1/(1+a)\). Since \(\sqrt{n!}\,[z^{n-i}]e^{b_az^2/2}\) is bounded by \(n^{i/2}|b_a|^{(n-i)/2}\sqrt{\binom{2p}{p}2^{-2p}}\), with \(2p=n-i\), one obtains
\(|\int x^{2\ell}g_n\varphi|\le C n^{\ell}|b_a|^{n/2}\). The fixed triangular inverse gives
\[
\|\widetilde g_n-g_n\|_\infty+\|\widetilde g_n-g_n\|_2\le C_K\,n^K|b_a|^{n/2}.
\]
At \(a=a_r\), \(b_a=(1-\sqrt{1-r^2})/(1+\sqrt{1-r^2})=\lambda_r^2\). Contraction of \(T_\rho\) on \(L^2(\varphi)\) gives
\[
\|T_{\sqrt s}\widetilde g_n\|\ge\|T_{\sqrt s}g_n\|-Cn^K\lambda_r^n,\qquad
\|T_{\sqrt r}\widetilde g_n\|\le\|T_{\sqrt r}g_n\|+Cn^K\lambda_r^n .
\]
Since \(\lambda_r^2<\lambda_r<\tau_a(s)\), the corrections are exponentially smaller than \(\tau_a(s)^{n/2}\) and \(\lambda_r^{n/2}\) respectively. Also \(\|\widetilde g_n\|_\infty\le B_n+o(1)\). Hence (S7) holds for \(\widetilde g_n\) with a smaller constant. \(\square\)

**Endpoint \(s=b(r)\).** We cannot decide whether (P6) holds at \(s=b(r)\). For each fixed \(a\) treated below, the sequence \(h_{n,a}\) does not refute it:

* Along \(a=a_r\), (S1) gives \(P_{n,a}(s)\le\tau^n[(n+1)^{-1/2}(1-\nu/\tau)^{-1/2}+n(\nu/\tau)^{n/2}]\), with \(\tau=\tau_a(s)>\nu=\nu_a(s)>0\) (\(\tau_a'>\nu_a'\) on \((r,s]\)). For this, split \(j\ge n/2\) and \(j<n/2\) and use \(c_m\le(2m+1)^{-1/2}\). With \(\tau^2=\lambda_r\), (S6) and (S8), the quotient in (S5) is at most \(C n^{-1/2}/n^{-1/4}=Cn^{-1/4}\to0\).
* For every \(a\in(1/(1+r),1/(1-r))\mathbin{\backslash}\{a_r\}\), all four rates \(\tau_a(r),\nu_a(r),\tau_a(1),\nu_a(1)\) are positive. Then (S2), (S6) and the strict inequality of Proposition S.4 bound the quotient by \(Cn^{1/2}(\varrho_a(s)^2/\varrho_a(r))^{n/2}\to0\).

For \(a\) outside that window, the denominator \(P_{n,a}(r)\) has terms of both signs, and we have not proved a matching lower bound. Luo's endpoint inequality ([6, Theorem F.1], with \(\|g\|_{L^2(dx)}\) in place of \(\|g\|_\infty\)) does not imply the bounded-input version, because \(\|g\|_{L^2(dx)}\) is not controlled by \(\|g\|_\infty\).

### III.1.4 Consequences for the method of Part I

**Corollary S.5 (optimality of the strict range; saturation at \(4/7\)).**

1. The strict range (P5) in Lemma I.2 is optimal up to the endpoint: (P6) holds for \(r<s<b(r)\) (Lemma I.2), fails for every \(s\in(b(r),1]\) (Proposition S), and is undecided at \(s=b(r)\). Failure persists within even inputs with any finite number of vanishing Gaussian moments.
2. Part I applies Lemma I.2 with \(r=q/2\), the parameter of the midpoint linear energy in (P16), and needs extra smoothing \(s>q\) so that \(\rho=q/s<1\) in (P15)--(P17). (A smaller \(r\) would only make the required inequality stronger, because \(\|T_{\sqrt r}g\|\) increases with \(r\) and \(b\) is increasing.) Since \(b\) is strictly increasing (\(b'(r)=(1+2r)/\sqrt{2r(1+r)}-1>0\)), the admissible interval \(q<s<b(q/2)\) is nonempty iff \(q<b(q/2)\), that is iff \(q<4/7\). For \(q>4/7\) every \(s>q\) lies in the failure range \(s>b(q/2)\) of Proposition S. At \(q=4/7\) every \(s>q=b(2/7)\) fails, and \(s=q\) is the undecided endpoint, which also gives \(\rho=q/s=1\) and no smoothing margin. The value \(4/7\) is the unique fixed point \(q=b(q/2)\): squaring \(\sqrt{q(1+q/2)}=3q/2\) gives \(1+q/2=9q/4\).
3. At \((r,s)=(2/7,4/7)\), the extremal parameter is \(a_{2/7}=(1-4/49)^{-1/2}=7/(3\sqrt5)\). This is exactly the parameter of Luo's counterexamples (M4). With \(\kappa=(3-\sqrt5)/2\),
\[
\tau_a(2/7)=\nu_a(2/7)=\lambda_{2/7}=\kappa^2=\frac{7-3\sqrt5}{2},\qquad
\tau_a(4/7)=\kappa,
\qquad b_a=\frac{a-1}{a+1}=\kappa^4.
\tag{S9}
\]
Indeed \(\lambda_{2/7}=2/(7+3\sqrt5)=(7-3\sqrt5)/2\), and \(a\cdot\frac37=1/\sqrt5\) gives \(\tau_a(4/7)=(\sqrt5-1)/(\sqrt5+1)=\kappa\). This number equals \(-\lambda_-\) of (M13) at \(q=4/7\).
4. For every \(a>0\), the quantity \(\kappa_a:=-\lambda_-(q)\) of (M13) equals \(\tau_a(q)\). Proposition S.4 with \((r,s)=(q/2,q)\) therefore shows that, for every \(q\le4/7\) and every \(a>0\),
\[
\tau_a(q)^2\le\varrho_a(q)^2\le\varrho_a(q/2),
\tag{S10}
\]
with equality throughout only for \(q=4/7\) and \(a=7/(3\sqrt5)\). The left side is the quantity \(\kappa_a^2\) that appears as the exponential rate of the quadratic term in (M23), wherever the saddle analysis of Part II applies to the parameter \(a\). The right side is the exact exponential rate of the linear midpoint energy \(\|L_{1/2}\|^2=4\|T_{\sqrt{q/2}}h_{n,a}\|^2\) along even \(n\). The upper bound is (S1). For the lower bound, the generating function of \(P_{n,a}\) gives
\[
P_{n,a}(u)=\frac1\pi\int_0^\pi\bigl(\tau_a(u)\cos^2\theta+\nu_a(u)\sin^2\theta\bigr)^n\,d\theta ,
\]
because \(\frac1\pi\int_0^\pi(1-x(\tau\cos^2\theta+\nu\sin^2\theta))^{-1}d\theta=((1-\tau x)(1-\nu x))^{-1/2}\). For even \(n\) the integrand is nonnegative, and integrating over a window of width \(n^{-1/2}\) around \(\theta=0\) or \(\theta=\pi/2\) (whichever endpoint attains \(\varrho_a(u)\)) gives \(P_{n,a}(u)\ge c_{a,u}n^{-1/2}\varrho_a(u)^n\). Thus the strict gap (M16) that drives Theorem B is unavailable at \(q\le4/7\) for every choice of \(a\). At \(q=4/7\) the two rates coincide, and they coincide only at Luo's parameter.

So the Part I method saturates at \(q=4/7\), and it does so on the same Hermite family \(h_{n,7/(3\sqrt5)}\) that produces Luo's counterexamples above \(4/7\). Item 4 is a statement about exponential rates only. Whether some other construction beats the linear term at \(q=4/7\) is not addressed.

## III.2 Luo's family at the boundary \(q=4/7\)

### III.2.1 Setting and statement

Fix \(q=4/7\), \(r=q/2=2/7\), \(a=7/(3\sqrt5)\), \(\kappa=(3-\sqrt5)/2\), an integer \(K\ge1\), and even \(n\). Let \(h_n=h_{n,a}\), let \(\widetilde h_n\) be its moment correction (M5), and for amplitudes \(\epsilon_n\in(0,1)\) put
\[
f_n=\varphi(1+\epsilon_n\widetilde h_n),\qquad
R_t=\epsilon_n\widetilde L_t+\epsilon_n^2\widetilde M_t,
\tag{E1}
\]
with \(\widetilde L_t,\widetilde M_t\) as in (M9)--(M10) at \(q=4/7\). If \(\epsilon_nB_n\to0\), then \(f_n\) is eventually a strictly positive even density matching Gaussian moments through degree \(2K\), and \(F_{4/7,f_n}(t)=h(\varphi)-\mathcal D(R_t)\). Denote by \(\hat h_k\) the coefficients of \(h_n\) in the standard basis \(e_k\) (only even \(k\) occur), and put
\[
S_n(u)=\|T_{\sqrt u}h_n\|_2^2=\sum_k\hat h_k^2u^k=\frac{2}{\sqrt{\Pi_a(u)}}P_{n,a}(u).
\]

**Proposition E (no midpoint convexity at the boundary).** Let \(\epsilon_n\to0\) with \(\epsilon_nn^{1/4}\to0\); for example \(\epsilon_n=e^{-\sqrt n}\) as in Theorem B. There are a closed interval \(J\) containing \(1/2\) in its interior, a constant \(c>0\) and \(n_0\) such that, for every even \(n\ge n_0\), \(t\mapsto F_{4/7,f_n}(t)\) is \(C^2\) on \(J\) and
\[
F_{4/7,f_n}''(t)\le-c\,\epsilon_n^2n^2\|L_t\|_2^2\le-c\sqrt{\tfrac{2}{1+a}}\;\epsilon_n^2n^2\kappa^{2n}<0
\qquad(t\in J).
\tag{E2}
\]
At the midpoint the exact leading term is
\[
F_{4/7,f_n}''(1/2)=-\frac{136}{15\sqrt{2(1+a)}}\;\epsilon_n^2\,n^2\kappa^{2n}\,\bigl(1+o(1)\bigr).
\tag{E3}
\]
Thus, under these amplitude conditions, Luo's construction at \(q=4/7\) has strictly negative curvature on a fixed neighbourhood of equal weighting. This is the opposite sign to (C), and the equal-weight point is a strict local maximum.

**What this does not show.** Proposition E does not decide whether \(4/7\in\mathcal C\). It shows only that Luo's construction, with this \(a\) and any vanishing amplitude satisfying \(\epsilon_nn^{1/4}\to0\), cannot produce a counterexample at the boundary with the failure near \(t=1/2\). It says nothing about \(t\notin J\), about other perturbation families, or about amplitudes not tending to zero.

### III.2.2 The linear term

Write the Hermite degree as \(k=2m\). Then \(L_t=\sum_k\hat h_k\,\ell_k(t)e_k\), with
\(\ell_k(t)=q^m(t^m+(1-t)^m)\). The derivative of \(\ell_k\) is exact: \(\ell_k''=q^mm(m-1)(t^{m-2}+(1-t)^{m-2})\).

**Lemma E.1 (degree multipliers).** For \(t\in[0,1]\) and \(m\ge0\),
\[
m(m-1)\,\ell_k\le\ell_k''\le4m(m-1)\,\ell_k\le k^2\ell_k,\qquad |\ell_k'|\le k\,\ell_k,
\tag{E4}
\]
with \(\ell_k''(1/2)=4m(m-1)\ell_k(1/2)=k(k-2)\ell_k(1/2)\) and \(\ell_k'(1/2)=0\).

*Proof.* Let \(x=t\), \(y=1-t\), \(x+y=1\). The lower bound is \(x^{m-2}\ge x^m\). For the upper bound, expand \((x+y)^2(x^{m-2}+y^{m-2})\). The four mixed terms are bounded using \((x^{m-1}-y^{m-1})(x-y)\ge0\) and \((x^{m-2}-y^{m-2})(x^2-y^2)\ge0\), giving \(x^{m-2}+y^{m-2}\le4(x^m+y^m)\) for \(m\ge2\). The case \(m\le1\) is trivial since \(\ell_k''=0\). Similarly \(x^{m-1}+y^{m-1}=(x+y)(x^{m-1}+y^{m-1})\le2(x^m+y^m)\) gives the first-derivative bound. The midpoint identities are direct. \(\square\)

**Lemma E.2 (Hermite-degree distribution of \(h_n\)).** Let \(U=[0.2,0.4]\). There are constants \(\mu_U>0\) and \(C\) such that, for \(u\in U\) and all \(n\):

1. the mean degree \(\mathsf m_n(u)=\sum_kk\hat h_k^2u^k/S_n(u)=uS_n'(u)/S_n(u)\) satisfies \(\mathsf m_n(u)\ge\mu_Un\);
2. \(S_n(ue^{1/n})\le C\,S_n(u)\) whenever \(ue^{1/n}\in[0.2,0.45]\);
3. \(\sum_k k^{2p}\hat h_k^2u^k\le C_p\,n^{2p}S_n(u)\) for each fixed \(p\), when \(ue^{1/n}\le0.45\).

At \(u=r\) the mean is exact:
\[
\mathsf m_n(r)=\frac{n(\alpha+\beta)}2+\frac{a^2r^2}{\Pi_a(r)},\quad
\alpha=\frac{r\tau_a'(r)}{\kappa^2}=\frac3{\sqrt5},\quad
\beta=\frac{r\nu_a'(r)}{\kappa^2}=\frac{\sqrt5}3,\quad
\frac{\alpha+\beta}2=a .
\tag{E5}
\]

*Proof.* On \([0.2,0.45]\), \(\tau_a\) and \(\nu_a\) are positive: \(\tau_a(u)>0\) iff \(u>1-1/a\approx0.042\), and \(\nu_a(u)>0\) for all \(u\ge0\) since \(a>1\). The derivatives \(\tau_a'=2a/(1+a-au)^2\) and \(\nu_a'=2a/(1+a+au)^2\) are positive. Put \(A=u\tau_a'/\tau_a\), \(B=u\nu_a'/\nu_a\). Then
\[
\frac{uP_{n,a}'}{P_{n,a}}=\sum_jw_j\,(jA+(n-j)B),\qquad
w_j\propto c_jc_{n-j}\tau_a^j\nu_a^{n-j}>0,\quad\sum w_j=1,
\]
so it lies between \(n\min(A,B)\) and \(n\max(A,B)\). Also \(u\frac{d}{du}\log\Pi_a^{-1/2}=a^2u^2/\Pi_a\in(0,1)\). Item 1 follows with \(\mu_U=\min_U\min(A,B)\). Item 2 follows from
\(\log S_n(ue^{1/n})-\log S_n(u)=\int_u^{ue^{1/n}}\mathsf m_n(v)\,dv/v\le n^{-1}\sup\mathsf m_n\le\sup(A,B)+1\).
Item 3 uses \(k^{2p}\le(2pn/e)^{2p}e^{k/n}\), which gives \(\sum k^{2p}\hat h_k^2u^k\le(2pn/e)^{2p}S_n(ue^{1/n})\), followed by item 2. At \(u=r\), \(\tau_a=\nu_a=\kappa^2\) by (S9), so \(w_j=c_jc_{n-j}\), and \(\sum_jjc_jc_{n-j}=n/2\) by symmetry. The values of \(\alpha,\beta\) follow by substitution. \(\square\)

**Lemma E.3 (linear energy near the midpoint).** Let \(J\) be a closed interval around \(1/2\) small enough that \(qt,\ q(1-t),\ q\sqrt{t(1-t)}\in U\) and \(e^{1/n_0}\cdot q\max(t,1-t)\le0.45\) for \(t\in J\). Then for \(t\in J\) and \(n\ge n_0\ge4/\mu_U\),
\[
\|L_t\|^2=S_n(qt)+S_n(q(1-t))+2S_n\big(q\sqrt{t(1-t)}\big)\ \ge\ S_n(r)=\sqrt{\tfrac2{1+a}}\;\kappa^{2n},
\tag{E6}
\]
\[
\langle L_t,L_t''\rangle\ge\frac{\mu_U^2}{8}\,n^2\|L_t\|^2,\qquad
\|\partial_t^jL_t\|\le C\,n^j\|L_t\|\quad(j=1,2).
\tag{E7}
\]
At the midpoint, exactly,
\[
\langle L_{1/2},L_{1/2}''\rangle
=4\sum_kk(k-2)\hat h_k^2r^k
=\frac{8\kappa^{2n}}{\sqrt{2(1+a)}}\Big(\tfrac{17}{15}n^2+\big(\tfrac{17}{15}-\tfrac{7\sqrt5}5\big)n+\tfrac{47}{30}-\tfrac{7\sqrt5}{10}\Big).
\tag{E8}
\]

*Proof.* Since \(k\) is even, \(\ell_k(t)^2=(qt)^k+(q(1-t))^k+2(q\sqrt{t(1-t)})^k\). This gives the identity in (E6). Each \(S_n\) is increasing, so \(\|L_t\|^2\ge S_n(q\max(t,1-t))\ge S_n(r)\). Moreover \(S_n(r)=2\Pi_a(r)^{-1/2}\kappa^{2n}\) by (S8), and \(\Pi_a(r)=1+2a+a^2(1-r^2)=2(1+a)\) because \(a^2(1-r^2)=1\).

By (E4), \(\langle L_t,L_t''\rangle=\sum\hat h_k^2\ell_k\ell_k''\ge\frac14\sum_k k(k-2)\hat h_k^2\ell_k^2\). Expand \(\ell_k^2\) as above. For each of the three values \(u\), Jensen's inequality and Lemma E.2 give \(\sum k(k-2)\hat h_k^2u^k\ge(\mathsf m_n^2-2\mathsf m_n)S_n(u)\ge\frac12\mu_U^2n^2S_n(u)\), since \(\mathsf m_n\ge\mu_Un\ge4\). Summing gives the first part of (E7). For the derivative bounds, (E4) gives \(\|L_t''\|\le\|N^2L_t\|\) and \(\|L_t'\|\le\|NL_t\|\). Lemma E.2(3), applied to the three terms of \(\ell_k^2\), gives \(\|N^pL_t\|\le Cn^p\|L_t\|\).

For (E8), use \(\ell_k\ell_k''(1/2)=k(k-2)\ell_k(1/2)^2=4k(k-2)r^k\). The sum is then \(4[(\vartheta^2-2\vartheta)S_n](r)\), where \(\vartheta=u\,d/du\). Write \(S_n=2\Pi_a^{-1/2}P_{n,a}\), with \(\vartheta\log\Pi_a^{-1/2}=a^2u^2/\Pi_a=:\pi_1\) and \(\vartheta^2\Pi_a^{-1/2}/\Pi_a^{-1/2}=2a^2u^2/\Pi_a+3a^4u^4/\Pi_a^2=:\pi_2\). At \(u=r\), using \(\tau_a=\nu_a\) and the weights \(c_jc_{n-j}\),
\[
\frac{\vartheta P}{P}=\frac{n(\alpha+\beta)}{2},\qquad
\frac{\vartheta^2P}{P}=\mathbb E\,(j\alpha+(n-j)\beta)^2+\frac n2(\vartheta A+\vartheta B).
\]
The convolution identities
\(\sum_jj(j-1)c_jc_{n-j}=[x^n]\tfrac34x^2(1-x)^{-3}=\tfrac38n(n-1)\) and \(\sum_jjc_jc_{n-j}=n/2\) give \(\mathbb Ej^2=(3n^2+n)/8\) and \(\mathbb E\,j(n-j)=(n^2-n)/8\). Direct differentiation gives \(\vartheta A(r)=\vartheta B(r)=0\) for this \(a\), and \(\pi_1(r)=\frac{7\sqrt5}{30}-\frac12\). Collecting terms in
\((\vartheta^2-2\vartheta)S/S=\vartheta^2P/P+2\pi_1\vartheta P/P+\pi_2-2\pi_1-2\vartheta P/P\) gives (E8). Its leading coefficient is \((3\alpha^2+3\beta^2+2\alpha\beta)/8=17/15\). \(\square\)

### III.2.3 The quadratic term and the moment correction

**Lemma E.4 (inherited estimates, valid at \(q=4/7\)).** Shrinking \(J\) if necessary, there is \(C\) such that, for \(t\in J\), \(j\le2\) and even \(n\),
\[
\|M_t\|\le C\,n^{-1/2}\kappa^{n},\qquad
\|\partial_t^jM_t\|\le C(1+n)^j\|M_t\|,\qquad\partial_tM_t\big|_{t=1/2}=0 .
\tag{E9}
\]

*Proof.* Lemma II.3, (M17), (M20), (M25) and Lemma II.5 (M27) of Part II require only strict positivity of the entries of (M21), and of the two-weight matrix in the proof of Lemma II.5, near \(u=0\). Those entries are built from \(d,e,c\) of (M14), which are positive at \(q=4/7\) because \(d(4/7)>0\). By (M20) and (M17), with Stirling's formula, \(\|M_t\|^2\le Cn^{-1}\Lambda(t)^{2n}\) uniformly on \(J\). By (M25), \(\Lambda(1/2)=\kappa\), \((\log\Lambda)'(1/2)=0\) and \((\log\Lambda)''(1/2)=-4\gamma(4/7)=-2\sqrt5<0\). Hence \(\Lambda\le\kappa\) on a small \(J\). (The radicand \(U^2t(1-t)-b\kappa\) equals \((\kappa-b)^2/4>0\) at \(t=1/2\).) The derivative bound is (M27). The vanishing first derivative follows from the exchange symmetry \(M_t=M_{1-t}\) for even \(h_n\): the Hermite coefficients (P14) are symmetric in \(t\leftrightarrow1-t\). \(\square\)

**Lemma E.5 (bilinear smoothing bound; moment correction).** For even \(f,g\in L^2(\varphi)\) and \(0<\rho<1\), the map \(t\mapsto\mathcal B_\rho(t;f,g)\) belongs to \(C^2([0,1];L^2(\varphi))\), and
\[
\sup_t\|\partial_t^j\mathcal B_\rho(t;f,g)\|\le j!\,\eta_\rho^{-j}\|f\|\,\|g\|\qquad(j=0,1,2).
\tag{E10}
\]
Consequently, with \(\delta_n=h_n-\widetilde h_n\), for \(t\in J\) and \(j\le2\),
\[
\|\partial_t^j(\widetilde L_t-L_t)\|+\|\partial_t^j(\widetilde M_t-M_t)\|\le C\,n^K\kappa^{2n}.
\tag{E11}
\]

*Proof.* The proof of (P11) and (P13) in Section I.3 uses the hypothesis "coefficients zero below degree four" only to make the output coefficients (P14) polynomials in the complex weight \(z\). For even inputs, only even \(j,k\) occur, so evenness suffices. The absolute-coefficient majorisation and Cauchy's estimate on the disk of radius \(\eta_\rho\) apply verbatim to the bilinear form (P14) with coefficients \(b_j,b'_k\). This gives \(\|\partial_t^j\mathcal B_\rho(t;h,h)\|\le j!\eta_\rho^{-j}\|h\|^2\) for even \(h\). The same majorisation with the two coefficient sequences \(|b_j|\) and \(|b'_k|\) reduces the bilinear form to \(\mathbb E[f_{\rm abs}(G_1)g_{\rm abs}(G_2)\mid\cdot]\), whose norm is at most \(\|f\|\,\|g\|\) by conditional-expectation contraction; Cauchy's estimate then gives (E10). The form is not symmetric: \(\mathcal B_\rho(t;g,f)=\mathcal B_\rho(1-t;f,g)\).

For (E11), use (M8) at \(a=7/(3\sqrt5)\): \(\|\delta_n\|_\infty+\|\delta_n\|_2\le Cn^Kb^{n/2}\) with \(b=\kappa^4\). The linear difference is \(-\sum_jc_{n,j}(T_{\sqrt{qt}}+T_{\sqrt{q(1-t)}})\psi_j\). On \(J\subset(0,1)\) its \(t\)-derivatives are fixed \(L^2\) functions multiplied by \(c_{n,j}=O(n^K\kappa^{2n})\). For the quadratic difference, \(\widetilde M_t-M_t=-\mathcal B_q(t;h_n,\delta_n)-\mathcal B_q(t;\delta_n,h_n)+\mathcal B_q(t;\delta_n,\delta_n)\). Apply (E10) with \(\rho=q=4/7<1\), together with \(\|h_n\|_2^2=S_n(1)\le2(1+2a)^{-1/2}\), which follows from (S1) and \(|P_{n,a}(1)|\le1\). \(\square\)

### III.2.4 Proof of Proposition E

Write \(\epsilon=\epsilon_n\), and for fixed \(t\in J\) write \(\mathcal L=\|L_t\|\). By (E6), \(\kappa^n\le C\mathcal L\). Hence the moment-correction errors (E11) are at most \(Cn^K\kappa^n\mathcal L=o(n^{-p}\mathcal L)\) for every \(p\).

*Regularity and the entropy formula.* By Lemma E.1, Lemma E.2 and (E10), \(R\in C^2(J;L^2(\varphi))\). Conditional expectation and (S6) give
\[
\|R_t\|_\infty\le\eta_n:=2\epsilon\|\widetilde h_n\|_\infty+\epsilon^2\|\widetilde h_n\|_\infty^2\le3\epsilon B_n
\]
eventually, and \(\eta_n\to0\) because \(\epsilon n^{1/4}\to0\). The direct differentiation argument at the end of Section I.5 then yields
\[
F''=-\langle R_{tt},\log(1+R)\rangle-\int\frac{R_t^2}{1+R}\varphi
\le-\langle R_{tt},R\rangle+\|R_{tt}\|\,\|\psi(R)\|,
\tag{E12}
\]
with \(\psi(u)=\log(1+u)-u\) and \(|\psi(u)|\le u^2\) for \(|u|\le\frac12\). At \(t=1/2\), the first derivative \(\partial_tR\) vanishes exactly: \(\widetilde L_t\) and \(\widetilde M_t\) are symmetric in \(t\leftrightarrow1-t\) because \(\widetilde h_n\) is even. Hence at the midpoint no term is dropped: \(F''(1/2)=-\langle R_{tt},R\rangle-\langle R_{tt},\psi(R)\rangle\), with \(|\langle R_{tt},\psi(R)\rangle|\le\|R_{tt}\|\,\|\psi(R)\|\).

*Sizes.* From (E7), (E9) and (E11), for \(n\ge n_0\),
\[
\|\widetilde L\|\le2\mathcal L,\quad
\|\widetilde L_{tt}\|\le Cn^2\mathcal L,\quad
\|\widetilde M\|\le Cn^{-1/2}\mathcal L,\quad
\|\widetilde M_{tt}\|\le Cn^{3/2}\mathcal L,
\]
\[
\langle\widetilde L,\widetilde L_{tt}\rangle\ge\langle L,L_{tt}\rangle-\|\widetilde L-L\|\|\widetilde L_{tt}\|-\mathcal L\|\partial_t^2(\widetilde L-L)\|
\ge\frac{\mu_U^2}{16}n^2\mathcal L^2 .
\]

*Assembly.* Expanding \(\langle R_{tt},R\rangle\),
\[
-\langle R_{tt},R\rangle
\le-\epsilon^2\frac{\mu_U^2}{16}n^2\mathcal L^2
+\epsilon^3\big(\|\widetilde L\|\|\widetilde M_{tt}\|+\|\widetilde M\|\|\widetilde L_{tt}\|\big)
+\epsilon^4\|\widetilde M\|\|\widetilde M_{tt}\|
\le-\epsilon^2n^2\mathcal L^2\Big[\frac{\mu_U^2}{16}-C\epsilon n^{-1/2}-C\epsilon^2n^{-1}\Big].
\]
The nonlinear remainder is bounded by
\[
\|R_{tt}\|\,\|\psi(R)\|\le\eta_n\|R\|\,\|R_{tt}\|
\le\eta_n(2\epsilon\mathcal L+C\epsilon^2\mathcal L)(C\epsilon n^2\mathcal L+C\epsilon^2n^{3/2}\mathcal L)
\le C\eta_n\,\epsilon^2n^2\mathcal L^2 .
\]
Since \(\epsilon\to0\) and \(\eta_n\to0\), (E12) gives \(F''(t)\le-(\mu_U^2/32)\epsilon^2n^2\mathcal L^2\) for all \(t\in J\) once \(n\) is large. The bound is uniform in \(t\), since all constants are. Combined with (E6), this proves (E2).

At \(t=1/2\), the same bounds show that every term other than \(-\epsilon^2\langle\widetilde L,\widetilde L_{tt}\rangle\) is \(O(\epsilon n^{-1/2}+\eta_n)\,\epsilon^2n^2\kappa^{2n}\). Moreover \(\langle\widetilde L,\widetilde L_{tt}\rangle=\langle L,L_{tt}\rangle(1+O(n^{K}\kappa^{n}))\). Then (E8) gives (E3): \(8\cdot\frac{17}{15}/\sqrt{2(1+a)}=136/(15\sqrt{2(1+a)})\). \(\square\)

### III.2.5 Remarks

**Remark E.6 (why the boundary is different; amplitude condition).** At \(q>4/7\), Theorem B compares \(\epsilon^4\|M\|^2\asymp\epsilon^4\kappa^{2n}/n\) with a linear energy of exponential order \(\tau(q/2)^n<\kappa^{2n}\). The quadratic term wins exponentially, so the exponential comparison survives any \(\epsilon_n=e^{-o(n)}\) (the density and entropy-remainder conditions would then need separate checking). At \(q=4/7\), (S8) makes the one-copy smoothed energy exactly \(S_n(r)=\sqrt{2/(1+a)}\,\kappa^{2n}\), and \(\|L_{1/2}\|^2=4S_n(r)\): it is of the same exponential order and loses no factor \(n^{-1/2}\). It also carries the curvature factor \(n^2\) of (E8). The \(\epsilon^2\) linear curvature then dominates the \(\epsilon^3\) cross terms (relative size \(\epsilon n^{-1/2}\)) and the \(\epsilon^4\) quadratic term (relative size \(\epsilon^2n^{-1}\)) for any \(\epsilon\to0\). The condition \(\epsilon n^{1/4}\to0\) enters only through \(\|R\|_\infty\le3\epsilon B_n\) in the logarithmic remainder. Using Cramér's uniform bound (Remark after Lemma S.3) instead of \(B_n\) would reduce the hypothesis to \(\epsilon_n\to0\). We record this only as a conditional remark, because we do not prove that inequality here.

**Remark E.7 (scope).** (i) Proposition E is local in \(t\). We have not proved negativity on all of \([0,1]\): near \(t\in\{0,1\}\), the entries of (M21) lose positivity (\(w\to0\)), and Lemma II.5 is not available. (ii) No log-concavity is asserted for these amplitudes. (iii) By (S10) and the lower bound in Corollary S.5(4), no other value of \(a\) recovers the exponential gap at \(q=4/7\). This is a statement about rates, not a proof that every \(a\) gives negative curvature. (iv) The result is consistent with both possibilities \(\mathcal C=(0,4/7)\) and \(\mathcal C=(0,4/7]\).

## III.3 Decay of the Part I radius at the boundary

**Remark R (order of \(\delta_Q\)).** Let \(x=4/7-Q\downarrow0\), and let \(\delta_Q\) be the radius (P21) built from (P7) and (P10) with \(s=(Q+b(Q/2))/2\). Then
\[
\delta_Q\sim C_\delta\,x^{9/4},\qquad
C_\delta=\frac{(49/96)^2\,(343/144)^{1/4}}{64\,C_J}\approx5.92873\times10^{-3},
\tag{R1}
\]
where
\[
C_J=\frac{c_{4/7}}{c_{2/7}^{1/2}(7/(5\pi))^{1/4}}\cdot\frac{14/15}{28/33}\cdot\pi^{1/4}\approx0.852984,
\qquad c_u=\frac1{2\pi\sqrt{1-u^2}} .
\]
So the exponent is \(p=9/4\).

*Derivation.* Each step is a first-order expansion of an explicit elementary function; all limits are at \((r,s)\to(2/7,4/7)\).

1. From \(b'(2/7)=5/6\), \(\frac{d}{dQ}[b(Q/2)-Q]=-7/12\). Hence \(b(r)-Q\sim\frac7{12}x\) and \(s-Q=b(r)-s\sim\frac7{24}x\).
2. In (P7), \(1/(1+r)\to7/9\), \(1/(1+s)\to7/11\) and \(E_2\to P_0-B_2^2/(2\beta)>0.46\) stay bounded below. Meanwhile \(E_1(r,b(r))=0\) identically, and \(\partial_sE_1=-B_1/(\beta(1-s)^2)\to-49/3\) (with \(P_0=7/5\), \(\beta=14/45\), \(B_1=14/15\)). Hence \(\varepsilon=\frac12E_1\sim\frac{343}{144}x\).
3. In (P10) every factor other than \((\pi/\varepsilon)^{1/4}\) has a finite positive limit, so \(J\sim C_J\varepsilon^{-1/4}\).
4. With \(\rho=Q/s\), \(1-\rho\sim\frac{49}{96}x\). Also \(\eta_\rho\sim(1-\rho)/4\) and \(K_\rho\sim16/(1-\rho)^2\). Hence \(A_Q\sim KJ/8+J/(8\eta^2)\sim4J/(1-\rho)^2\).
5. In (P21), \(1/(16A_Q)\sim(1-\rho)^2/(64J)\) is the smallest entry. The others are \(1/12\), \(1/J\asymp x^{1/4}\) and \(8/(KJ)\sim(1-\rho)^2/(2J)\).

Combining these gives (R1). The exponent splits as \(2\) from the derivative smoothing (\(K_\rho\) and \(\eta_\rho^{-2}\), both \(\asymp(1-\rho)^{-2}\)) plus \(1/4\) from the weighted interpolation (\((\pi/\varepsilon)^{1/4}\) with \(\varepsilon\asymp x\)).

*Numerical confirmation (diagnostic).* In 60-digit Decimal arithmetic, `checks/extension_checks.py` finds \(\delta_Q/x^{9/4}=C_\delta(1+1.365\,x+\dots)\) for \(x=10^{-3},\dots,10^{-14}\). The local exponent between successive decades is \(2.25000000000005\) at \(x=10^{-14}\). The script also reproduces the four radii of the Part I table. Both quantities concern the sufficient radius of this proof, not the true size of the concavity neighbourhood.




# Discussion

## D.1 Which assumptions matter

- **Variance one** is a normalisation. Concavity in \(t\) is unchanged by rescaling, since \(h(cZ)=h(Z)+\log c\), and the noise level is absorbed into \(q\).
- **Identical inputs** are used essentially. They give the symmetry \(t\leftrightarrow1-t\) that makes \(t=1/2\) a critical point, which Theorem B exploits. Concavity for two different input laws is a different question.
- **Evenness** is necessary for a positive threshold. Luo's asymmetric counterexamples [6, Theorem 1.1] destroy concavity near the endpoints for every \(q>0\).
- **The neighbourhood.** The identity \(\sup\mathcal C=4/7\) holds for uniform relative-density neighbourhoods (\(\|f/\varphi-1\|_\infty\) small), and equally for any fixed relative \(C_b^m\) norm. The upper bound \(4/7\) also holds for weaker notions of closeness, such as relative entropy or Wasserstein distance, because Luo's sequences converge in those senses too. The lower bound does not transfer automatically: a perturbation can be small in relative entropy while its density ratio is far from one at some points. Fixed directions behave differently again; see the discussion of Abu-Ajamieh's \(\mathrm{He}_4\) family in Section 1.3.

## D.2 Why Schubr's argument stops at one half

Schubr's neighbourhood theorem [7] uses, in place of Lemma I.2, a semigroup bound of the form
\[
\|T_{\sqrt s}g\|_2^2=\langle T_{\sqrt r}g,\,T_{s/\sqrt r}g\rangle\le\|g\|_2\,\|T_{\sqrt r}g\|_2,
\]
which requires \(s\le\sqrt r\). With \(r=q/2\) and \(s>q\), this needs \(q<\sqrt{q/2}\), that is \(q<1/2\). The restriction is genuine for any estimate of this \(L^2(\varphi)\) form. The Hermite polynomial \(g=e_k\) gives \(\|T_{\sqrt s}e_k\|_2^2/\|T_{\sqrt r}e_k\|_2=(s/\sqrt r)^k\), which is unbounded when \(s>\sqrt r\). The improvement to \(4/7\) uses boundedness of the perturbation: the polynomials \(e_k\) are unbounded, and the bounded test functions of Proposition S push the threshold out to \(b(r)\ge\sqrt r\), with equality only at \(r\in\{0,1\}\).

## D.3 A reading in blind source separation

Because a rotation of independent standard Gaussians is again a pair of independent standard Gaussians, for \(t=\cos^2\theta\),
\[
\sqrt t\,X_q+\sqrt{1-t}\,Y_q=\sqrt q\,(\cos\theta\,X+\sin\theta\,Y)+\sqrt{1-q}\,G',
\]
with \(G'\) standard Gaussian. So \(\theta\mapsto F_{q,f}(\cos^2\theta)\) is the entropy of a one-dimensional projection of a whitened two-source mixture observed in isotropic Gaussian noise at signal-to-noise ratio \(q/(1-q)\). Minimising the entropy of projections is the classical contrast of independent component analysis [16]. The map \(\theta\mapsto t=\cos^2\theta\) is a monotone bijection from \([0,\pi/2]\) to \([0,1]\), and at \(\theta=\pi/4\) one has \(dt/d\theta=-1\) and \(F'(1/2)=0\), so \(d^2F/d\theta^2=F''(1/2)\) there.

The two theorems then read as follows, for symmetric near-Gaussian sources.

- **Signal-to-noise ratio below \(4/3\).** For non-Gaussian sources within the neighbourhood of Theorem A, the projected entropy has no local minimum strictly inside \((0,\pi/2)\): a strictly concave function of \(t\) has no interior local minimum. The only minima on the quarter circle are the source directions. (Gaussian sources give a constant contrast.)
- **Signal-to-noise ratio above \(4/3\).** By Theorem B, there are smooth, strongly log-concave sources, matching any prescribed number of Gaussian moments, for which the \(45^\circ\) mixing direction is a strict local minimum of the entropy contrast. Such spurious mixing minima were previously known for strongly multimodal sources [17]. Because the sources match Gaussian moments to any fixed order, contrasts built from finitely many cumulants cannot see them.

Three caveats apply. The effect is minute (Section II.6). It concerns the exact population contrast, not an estimator. And in noisy models the minimiser of the contrast need not coincide with the noiseless source direction.

## D.4 Relation to other Gaussian tools

The classical concavity of entropy power along the heat flow [13, 14] and the I-MMSE relation [15] move along the noise direction: they vary \(q\) with the input fixed. Here \(q\) is fixed and the mixing weight varies. The obstruction is a competition between the linear and quadratic Hermite terms of the output density, which those identities do not separate. Theorem A has the form of a local quantitative stability statement near the Gaussian, in the spirit of stability results for the entropy power inequality [18], but for a different functional.

## D.5 Open problems

1. **The boundary.** Is \(4/7\in\mathcal C\)? Proposition E shows that Luo's family cannot decide this near equal weighting: at \(q=4/7\) its linear term still dominates there. Lemma I.2 is undecided at its endpoint \(s=b(r)\), and the smoothing margin \(s-q\) of Part I vanishes there. A proof that \(4/7\in\mathcal C\) would need a derivative estimate directly against \(\|NL\|\), avoiding the extra smoothing. A counterexample would need a family whose quadratic interaction beats the linear term without the exponential gap of (M16), which (S10) shows is unavailable, at the level of exponential rates, for every Hermite family \(h_{n,a}\).
2. **The global problem.** Luo's Problem 10.1, concavity for every symmetric log-concave input at a fixed noise level, is untouched.
3. **Other functionals and more copies.** Theorem B's leading term is a \(\chi^2\)-type quantity, which suggests that the same threshold may govern Rényi entropies near the Gaussian. Mixtures of \(k\ge3\) copies may have a different threshold. Neither question is addressed here.
4. **Sharper radii.** The radius \(\delta_Q\) of Theorem A is far from optimal, and its decay rate \((4/7-Q)^{9/4}\) is a property of this proof. The true size of the concavity neighbourhood near \(4/7\) is unknown.




# References

1. Ball, K., Nayar, P., & Tkocz, T. (2016). A reverse entropy power inequality for log-concave random vectors. *Studia Mathematica*, 235(1), 17--30. <https://doi.org/10.4064/sm8418-6-2016>

2. Artstein, S., Ball, K., Barthe, F., & Naor, A. (2004). Solution of Shannon's problem on the monotonicity of entropy. *Journal of the American Mathematical Society*, 17(4), 975--982. <https://doi.org/10.1090/S0894-0347-04-00459-X>

3. Lieb, E. H. (1978). Proof of an entropy conjecture of Wehrl. *Communications in Mathematical Physics*, 62(1), 35--41. <https://doi.org/10.1007/BF01940328>

4. Eskenazis, A., Nayar, P., & Tkocz, T. (2018). Gaussian mixtures: Entropy and geometric inequalities. *The Annals of Probability*, 46(5), 2908--2945. <https://doi.org/10.1214/17-AOP1242>

5. Eskenazis, A., & Gavalakis, L. (2024). On the entropy and information of Gaussian mixtures. *Mathematika*, 70(2), e12246. <https://doi.org/10.1112/mtk.12246>

6. Luo, C. (2026). *Near-Gaussian counterexamples to the Ball-Nayar-Tkocz entropy concavity conjecture* (Version 2, 14 September 2026) [Preprint]. arXiv:2609.11418. <https://arxiv.org/abs/2609.11418v2>

7. Schubr, S. (2026, October 30). *Mixing-weight entropy concavity near the Gaussian below one-half signal* (Public edition 3) [Preprint, not peer reviewed]. Luminary Research Archive. <https://research.luminaryoddities.org/local-and-eventual-entropy-concavity-under-gaussian-smoothing/v/3>

8. the-omega-institute/trureturing (2026, November 1). *Open problem (preregistration): prove Congyi Luo Problem 10.3 (arXiv:2609.11418v2): positive midpoint entropy curvature for every \(4/7<q\le2-\sqrt2\)* [GitHub issue #25275, opened 1 November 2026; proposed argument, not refereed]. <https://github.com/the-omega-institute/trureturing/issues/25275>

9. Abu-Ajamieh, F. (2026). *Noise-induced restoration of entropy concavity near Gaussian equilibrium* (Version 1) [Preprint]. arXiv:2610.22386. <https://arxiv.org/abs/2610.22386v1>

10. Le, T. (2014). Normal and isometric weighted composition operators on the Fock space. *Bulletin of the London Mathematical Society*, 46(4), 847--856. <https://doi.org/10.1112/blms/bdu046>

11. Carroll, T., & Gilmore, C. (2021). Weighted composition operators on Fock spaces and their dynamics. *Journal of Mathematical Analysis and Applications*, 502(1), 125234. <https://doi.org/10.1016/j.jmaa.2021.125234>

12. Abbe, E., & Zheng, L. (2012). A coordinate system for Gaussian networks. *IEEE Transactions on Information Theory*, 58(2), 721--733. <https://doi.org/10.1109/TIT.2011.2169536>

13. Costa, M. (1985). A new entropy power inequality. *IEEE Transactions on Information Theory*, 31(6), 751--760. <https://doi.org/10.1109/TIT.1985.1057105>

14. Villani, C. (2000). A short proof of the "concavity of entropy power". *IEEE Transactions on Information Theory*, 46(4), 1695--1696. <https://doi.org/10.1109/18.850718>

15. Guo, D., Shamai, S., & Verdu, S. (2005). Mutual information and minimum mean-square error in Gaussian channels. *IEEE Transactions on Information Theory*, 51(4), 1261--1282. <https://doi.org/10.1109/TIT.2005.844072>

16. Comon, P. (1994). Independent component analysis, A new concept? *Signal Processing*, 36(3), 287--314. <https://doi.org/10.1016/0165-1684(94)90029-9>

17. Vrins, F., Pham, D.-T., & Verleysen, M. (2007). Mixing and non-mixing local minima of the entropy contrast for blind source separation. *IEEE Transactions on Information Theory*, 53(3), 1030--1042. <https://doi.org/10.1109/TIT.2006.890716>

18. Courtade, T. A., Fathi, M., & Pananjady, A. (2018). Quantitative stability of the entropy power inequality. *IEEE Transactions on Information Theory*, 64(8), 5691--5703. <https://doi.org/10.1109/TIT.2018.2808161>

# Research and verification status

This is an unrefereed research release, version 0.1.0-candidate, archived at Zenodo with DOI [10.5281/zenodo.23910629](https://doi.org/10.5281/zenodo.23910629) (all versions: [10.5281/zenodo.23910628](https://doi.org/10.5281/zenodo.23910628)); source and checks: <https://github.com/ipitchford/entropy-concavity-four-sevenths>. The manuscript was developed in an AI-assisted mathematical research session; its author is listed as Anonymous. The proofs are analytic and have not been checked by a proof assistant.

Before release, the manuscript received an internal AI-assisted review: five independent reviewer agents, a separate adversarial review by a different model, a targeted audit of the new and rewritten proofs, and a clean replay of every shipped script from a fresh copy. This is internal checking, not external peer review or independent human reproduction.

The shipped checkers (`checks/`) use only the Python standard library. They test exact finite identities, the closed-form constants, the radius table and finite-\(n\) asymptotics, and they fail if any check fails, including under `python -O`. `checks/negative_controls.py` confirms that 16 deliberately broken copies of the first three checkers fail in both modes; `checks/extension_checks.py` carries its own internal negative controls. The checks support the analytic arguments; they do not replace them.

**Version history.** The first internal version (2 November 2026) contained Theorems A and B. This revision adds Propositions S and E and Remark R, the discussion of the threshold and of blind source separation, the attribution to the public proposal [8] for Luo's Problem 10.3, and an expanded bibliography. It also rewrites the reproducibility scripts to fail closed.

The literature comparison was a bounded search. It found the earlier public proposal [8] for Problem 10.3, cited above. It found no earlier statement of Theorem A, the Corollary, or Propositions S and E. A bounded search does not certify priority.



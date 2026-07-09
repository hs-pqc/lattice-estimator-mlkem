# From Empirical Formulas to Search Resolution: Rethinking ζ_optimal in Dual Hybrid Attacks on LWE

**Status:** Research Note (Work in Progress)  
**Author:** Hangshin Cho (hs-pqc)  
**Date:** June 2026 (initial formulas); July 2026 (reframing)

---

## Abstract

We initially set out to find a closed-form empirical formula for the
optimal splitting dimension ζ in dual hybrid attacks on LWE, as a function
of lattice dimension n, modulus q, and secret distribution parameter η.
Systematic experiments across 42 parameter sets (ML-KEM, ML-DSA, NTRU+,
HAETAE) produced two candidate formulas with reasonable fit (max_err ≤ 5).

However, closer analysis of *why* the two cost-model call paths in
lattice-estimator — `LWE.dual_hybrid` (indirect, MATZOV) and
`dual_hybrid` (direct) — disagree on ζ by as much as 20–40 units for the
same parameters revealed that **treating ζ_optimal as a function of
(n, q, η) alone is a misframing**. ζ is determined jointly by (a) which
call path is used, since the two paths encode different assumptions about
whether an FFT distinguisher is realistically implementable, and (b) the
optimizer's search structure — in particular, the indirect call's coarse
`opt_step` grid (~8–10 unit steps), which appears to report ζ=0 for
ML-KEM-512 and ML-KEM-1024 as a search-resolution artifact rather than a
true optimum. A parameter-indexed formula cannot capture either effect,
because it silently assumes the search that produced its training data was
resolved finely enough to trust.

This reframes the useful question from *"what is ζ_optimal(n, q, η)?"* to
*"how does search resolution affect the reliability of the resulting rop
(security-bit) estimate?"* — a question about the estimator's behavior,
not about the lattice parameters. Under this lens, the ~10-bit spread we
observe in ML-KEM-768's security margin between the two call paths is
better understood as an artifact of assumption/resolution mismatch than
as a property of ML-KEM-768 itself. The original formulas are retained
below (Section 5) as a documented first attempt, superseded by this
reframing.

## 1. Introduction

The dual hybrid attack combines lattice reduction (BKZ) with exhaustive search
over a subset of secret coordinates. The splitting dimension ζ controls the
trade-off between BKZ cost and guessing cost:

- **Large ζ:** fewer dimensions for BKZ → cheaper lattice reduction,
  but more secret coordinates to guess
- **Small ζ:** more dimensions for BKZ → expensive lattice reduction,
  but fewer coordinates to guess

The optimal ζ minimizes total attack cost. In practice, lattice-estimator
finds ζ numerically via binary search (local_minimum in lwe_dual.py).

**Original Research Question (superseded):**  
Can ζ_optimal be expressed as a closed-form function of (n, q, η)?

**Current Research Question:**  
How does the resolution of the ζ search (grid step, call path, cost-model
assumptions) affect the reliability of the resulting rop (security-bit)
estimate? Put differently: when lattice-estimator reports a security
level, how much of that number reflects the cryptographic hardness of the
instance, and how much reflects the coarseness of the search that produced
it?

---

## 2. Preliminaries

### 2.1 LWE Parameters

LWE instance: (A, b = As + e) where:
- n: lattice dimension
- q: modulus
- Xs = Xe = CenteredBinomial(η): secret/error distribution
  - Support: {-η, ..., η}
  - Variance: σ² = η/2
  - Standard deviation: σ = √(η/2)

### 2.2 Dual Hybrid Attack Cost

Total cost = BKZ cost + Guessing cost
BKZ cost    ≈ 2^(0.292β)

Guessing cost ≈ (2η+1)^ζ
Optimal ζ satisfies:
d/dζ [BKZ_cost(n-ζ) + Guessing_cost(ζ)] = 0
### 2.3 lattice-estimator

Open-source tool by Albrecht et al.
Uses local_minimum binary search to find optimal (β, ζ) numerically.

**Key finding:** We use `dual_hybrid` direct call (not `LWE.dual_hybrid`)
for accurate results. `LWE.dual_hybrid` uses MATZOV cost model which
gives estimates consistent with the official MATZOV security analysis.

---

## 3. Experimental Setup

- Tool: SageMath + lattice-estimator 0.1.0
- Method: `dual_hybrid` direct call
- Parameters tested: 42 parameter sets across:
  - ML-KEM (512/768/1024)
  - ML-DSA (44/65/87)
  - NTRU+ (576/768/864/1152)
  - Custom parameters (various n, q, η)
  - FHE-type parameters (q=2^32, 2^60)
  - Extreme parameters (η=1, 8, 16; small q=1021)

---

## 4. Key Observations

### 4.1 Primal vs Dual

Accurate cost model shows primal attack is more efficient than dual
across all three standards:

| Standard | Primal (bits) | Dual (bits) |
|---|---|---|
| ML-KEM-768 | 204.9 | 206.4 |
| ML-DSA-65 | 231.3 | 233.9 |
| NTRU+768 | 197.7 | 199.1 |

### 4.2 Cost Model Comparison: ζ values

Two cost models give significantly different ζ values:

| Parameter | LWE.dual_hybrid (MATZOV) | dual_hybrid direct | Difference |
|---|---|---|---|
| ML-KEM-512 | 0 | 20 | +20 |
| ML-KEM-768 | 20 | 32 | +12 |
| ML-KEM-1024 | 0 | 41 | +41 |
| ML-DSA-44 | 10 | 23 | +13 |
| ML-DSA-65 | 20 | 25 | +5 |
| ML-DSA-87 | 30 | 44 | +14 |
| NTRU+576 | 20 | 25 | +5 |
| NTRU+1152 | 30 | 43 | +13 |

MATZOV uses FFT distinguisher which makes ζ=0 optimal for some parameters.
dual_hybrid direct call finds larger ζ values without FFT assumption.
ML-KEM-512 and ML-KEM-1024 show ζ=0 under MATZOV — meaning no dimension
reduction is optimal when FFT is available.

**Caveat (see Section 5.0 and Open Problem 0):** the indirect call's
default `opt_step` (~8–10) means ζ=0 could also simply be the coarsest
grid point evaluated, rather than a confirmed optimum. We have not yet
re-run this comparison at `opt_step=1` to rule this out — this is the
single most important unresolved check in this note, since it would
determine whether the FFT-assumption story above is the full explanation
or only part of it.
### 4.3 ζ scales with n/log2(q)

ζ/n ≈ 0.030 + 0.096/√(n/log2(q))

This suggests ζ ≈ 0.030n + 0.096√(n × log2(q)) for small q.

### 4.4 C(η) follows η^(1/3)

When normalizing ζ × log2(q) / n, the resulting constant C varies
with η following approximately η^(1/3):

| η | C = ζ×log2(q)×η^(1/3)/n |
|---|---|
| 1 | 0.625 |
| 2 | 0.597 |
| 3 | 0.587 |
| 4 | 0.575 |
| 8 | 0.561 |

C converges to ≈ 0.597 for η ≥ 2.

---

## 5. Initial Formulas (Superseded) and Why the Approach Was Misframed

### 5.0 Why a Parameter-Indexed Formula Doesn't Work

The formulas below were fit to ζ values produced by the `dual_hybrid`
direct call at a fixed `opt_step`. Two problems only became visible after
the fact:

1. **Call-path dependence.** `LWE.dual_hybrid` (indirect, MATZOV) and
   `dual_hybrid` (direct) do not disagree because one is "more accurate" —
   they encode different assumptions about FFT distinguisher
   implementability (Section 4.2). A formula fit to one call path's output
   is a formula for *that assumption*, not for ζ_optimal in general. Fitting
   across both paths in one formula is not meaningful, since they are not
   estimating the same quantity.

2. **Search-resolution dependence.** The indirect call's default
   `opt_step` (~8–10) is coarse enough that for ML-KEM-512 and ML-KEM-1024
   it likely never evaluates the true optimum and reports ζ=0 by default
   instead of a small nonzero value. Any formula trained on this output
   inherits the artifact as if it were signal.

Because both effects are about the *search and cost-model machinery*,
not about (n, q, η), no formula indexed only on parameter values can be
correct in principle — it will always be a formula for "whatever the
grid happened to find," which changes if the grid changes. This is the
core reframing of this note: the two formulas below are kept as a
documented first attempt, not as the intended contribution.

### Formula 1: Practical Formula

> **ζ ≈ floor((-0.076 + 0.701 × η^0.2) × n / (log2(q) × log2(η+1))) + 2**

- Derived via least-squares fitting
- Valid: n≥512, η≥1, q=1021~2^60
- Accuracy: max_err=3, avg_err=1.00 (42 parameter sets)

### Formula 2: Mathematically Motivated

> **ζ ≈ floor(0.597 × n / (log2(q) × η^(1/3))) + 2**

- Derived from minimizing variance of C(η) over α
- Optimal α = 1/3 (cube root)
- Valid: n≥512, η≥2, q=1021~2^60
- Accuracy: max_err=5, avg_err=1.29 (31 parameter sets)

---

## 6. Mathematical Intuition

### 6.1 Why η^(1/3)?

CenteredBinomial(η) has standard deviation σ = √(η/2).

η^(1/3) = (η/2)^(1/3) × 2^(1/3) ∝ σ^(2/3)

**Conjecture:** The guessing cost in dual hybrid attack scales as σ^(2/3),
not σ or σ².

This is consistent with the balance condition:
0.292 × β(n-ζ) ≈ ζ × log2(guessing_cost_per_coordinate)
where guessing_cost_per_coordinate ∝ σ^(2/3) = η^(1/3).

### 6.2 Why C ≈ 0.597?

0.597 ≈ 0.292 × 2.044

0.292 is the BKZ exponent (Albrecht et al.).
The factor 2.044 may relate to the ratio between
lattice dimension and effective guessing dimension.

**Open Question:** What is the exact mathematical origin of C ≈ 0.597?

---

## 7. Open Problems

0. **[Primary] Search resolution vs. rop reliability:**
   For a fixed instance, how does varying `opt_step` (grid coarseness)
   change the reported rop? Does the reported security level converge
   as `opt_step` → 1, or does it plateau early? If ML-KEM-512/1024's
   ζ=0 result is a resolution artifact, does refining the grid change
   their reported bit-security, or only ζ itself (i.e., is the cost
   function flat enough near ζ=0 that rop is robust even when ζ is
   not)? This determines whether the ~10-bit ML-KEM-768 spread
   (Section 11) is a real modeling uncertainty or a fixable estimator
   artifact.

1. **Mathematical proof of η^(1/3):**
   Prove that guessing cost in dual hybrid scales as σ^(2/3).
   Direct verification shows LHS/RHS ≈ 2~3 (not 1), suggesting
   the balance condition is more complex than BKZ_cost = Guess_cost.

2. **Origin of C ≈ 0.597:**
   0.597 ≈ 0.292 × 2.044. The factor 2.044 may relate to
   repetition cost, sieving dimension, or memory cost in BKZ.
   Deriving C analytically from the full cost function remains open.

3. **η=1 anomaly:**
   Formula 2 has larger errors for η=1 (max_err=5).
   Empirically, ζ(η=1)/ζ(η=2) ≈ 1.38 consistently across n.
   This ratio is larger than η^(1/3) predicts (2^(1/3) ≈ 1.26).
   Conjecture: η=1 sparsity introduces an additional factor of (4/3)^(1/3).

4. **Cost model independence:**
   All experiments use MATZOV cost model (default in lattice-estimator 0.1.0).
   ADPS16, Kyber, and quantum cost models are not compatible
   with dual_hybrid in this version. Verification under other
   cost models remains open.

5. **SMAUG-T analysis:**
   SMAUG-T uses sparse secret (HWT distribution) with small h.
   lattice-estimator 0.1.0 does not support SparseTernary(n, h)
   for small h values. Separate analysis needed.

6. **Connection to security margins:**
   How does ζ_optimal relate to actual security margins
   in NIST/KPQC PQC standards?

7. **Generalization beyond CBD:**
   All experiments use CenteredBinomial(η).
   Does the formula extend to DiscreteGaussian or other distributions?

## 8. Experimental Data

All experimental results are in the `results/` directory.

### Core results (`results/1_core/`)

| File | Content |
|---|---|
| `all_security_revised.txt` | Primal/dual comparison (all standards) |
| `zeta_accurate.txt` | MATZOV vs direct call ζ comparison |
| `zeta_final.txt` | Formula 1 blind test |
| `zeta_verify3.txt` | Formula 1 validation (17 parameter sets) |
| `zeta_cube_verify.txt` | Formula 2 validation (31 parameter sets) |
| `kpqc_test.txt` | HAETAE validation |
| `hints_threshold.txt` | ML-KEM-512 hint threshold analysis |
| `hints_all.txt` | ML-KEM hint threshold comparison |
| `C_vs_nlogq.png` | C vs n/log2(q) visualization |

### Analysis (`results/2_analysis/`)
Intermediate experiments and fitting results.

### Deprecated (`results/3_deprecated/`)
Early results using LWE.dual_hybrid (MATZOV) — kept for reference.

## 9. References

1. Albrecht et al., "On the Concrete Hardness of LWE", JMC 2015
2. MATZOV, "Report on the Security of LWE", 2022
3. Espitau, Joux, Kharchenko, "On a Dual/Hybrid Approach", INDOCRYPT 2020
4. NIST FIPS 203 (ML-KEM), 2024
5. NIST FIPS 204 (ML-DSA), 2024
6. lattice-estimator, github.com/malb/lattice-estimator
7. Cheon et al., "HAETAE: Shorter Lattice-Based Fiat-Shamir Signatures", CHES 2024
8. Cheon et al., "SMAUG-T: Post-Quantum KEM", KpqC 2025
9. Hhan, Hong, Kim, Lee, Lee, "From Perfect to Approximate Hints: Efficient LWE Secret Recovery Leveraging Low Hamming Weight", S&P 2026 (ePrint 2026/1081) — co-authored by Changmin Lee
10. Kim, Lee, Kim, Lee, "SQIsign with Fixed-Precision Integer Arithmetic", PKC 2026 (ePrint 2025/1649) — Changmin Lee, corresponding author

## 10. Connection to Approximate Hints

### Motivation

If FFT distinguisher (MATZOV) is realistic, ML-KEM-512 has only
11.7-bit security margin above 128 bits. Combined with approximate
hints (Hhan et al., S&P 2026), this margin may shrink further.

### Experiment: hint threshold analysis

We approximate hint effects by reducing effective dimension n → n-h.
This is a rough model; accurate analysis requires DBDD framework.

| Parameter | Base security | Margin | Threshold |
|---|---|---|---|
| ML-KEM-512 | 139.7 bits | 11.7 bits | ~44 hints |
| ML-KEM-768 | 196.4 bits | **4.4 bits** | **~16 hints** |
| ML-KEM-1024 | 262.3 bits | 6.3 bits | ~21 hints |

## 11. FFT Distinguisher and Security Margin Uncertainty

### Key Claim

The practical security of ML-KEM-768 critically depends on whether
the FFT distinguisher (MATZOV) is realistically implementable.

| Parameter | Target | FFT available | FFT unavailable | Uncertainty |
|---|---|---|---|---|
| ML-KEM-512 | 128 bits | 11.7 bits | 17.5 bits | 5.9 bits |
| ML-KEM-768 | 192 bits | **4.4 bits** | 14.4 bits | **10.0 bits** |
| ML-KEM-1024 | 256 bits | 6.3 bits | 21.5 bits | 15.2 bits |

### Interpretation

- If FFT is realistic → ML-KEM-768 margin is only 4.4 bits (most vulnerable)
- If FFT is unrealistic → ML-KEM-768 margin is 14.4 bits (safe)
- ML-KEM-768 has the largest security uncertainty (10 bits) among all three

### Open Question

What is the realistic implementation cost of the FFT distinguisher?
This single question determines whether ML-KEM-768 is the weakest
or the middle parameter set among the three ML-KEM variants.

Combined with approximate hints (Section 10), if FFT is realistic
and 16 hints are obtainable, ML-KEM-768 NIST Level 3 security
may be at risk.

### Key Finding

**ML-KEM-768 has the smallest margin (4.4 bits) among all three
parameter sets.** Only ~16 approximate hints are needed to break
NIST Level 3 security (192 bits).

### Open Question

Can an attacker realistically obtain 16 approximate hints against
ML-KEM-768 in practice? If so, NIST Level 3 security may be
at risk under MATZOV + approximate hints combined attack.

This connects directly to Hhan et al. S&P 2026 which analyzes
the cost reduction from approximate vs perfect hints.
## Appendix: MATZOV ζ Search Resolution — Mechanism and Verification (2026-07)

### Background

`LWE.dual_hybrid` in lattice-estimator resolves to `MATZOV.__call__`
(`estimator/lwe_dual.py`), which searches over ζ (k_enum) and t (k_fft)
using `early_abort_range(..., step=10)` — a **hardcoded** grid step,
independent of the public `opt_step` parameter. `early_abort_range` is a
greedy search: it terminates as soon as the cost increases once. This
appendix documents a full reconstruction of the ζ=20 vs ζ=32 discrepancy
originally reported in Section 4, and identifies the exact mechanism
behind it.

### Experiment 1 — Isolating the MATZOV ζ search

For ML-KEM-768, the default MATZOV call finds ζ=20 (log2(rop)=196.4).
Adjusting the internal grid step for ζ alone (with t held at the default
step of 10) still gives ζ=4 — worse, not better, than the default. This
rules out simple "finer grid → better answer" as the story; the search
is path-dependent, not merely under-resolved.

### Experiment 2 — Full ζ enumeration under MATZOV (ML-KEM-768)

Scanning ζ=0..768 (step=4), with p and t re-optimized at each ζ:
- Local non-convexity for ζ=0..32 (amplitude ≈1.7 bits)
- Monotonic increase beyond ζ≈36
- **True global minimum: ζ=24, log2(rop)=196.133**
- Default MATZOV output (ζ=20) is only 0.23 bits above the true optimum

### Experiment 3 — Reproducing the original ζ=32 (direct call, no FFT)

The original ζ=20 vs 32 comparison in Section 4 was between MATZOV
(FFT distinguisher assumed) and the direct `dual_hybrid` call (FFT
distinguisher **not** assumed, `fft=False`), not between two grid
resolutions. A full enumeration of the direct call reproduces the
recorded optimum exactly:

**Global minimum: ζ=32, log2(rop)=206.357** — matches Section 4's
recorded value (ζ=32, rop=206.4) exactly.

The direct-call cost curve is smooth and strictly unimodal (monotonic
decrease then increase around ζ=32), unlike the MATZOV curve.

### Experiment 4 — opt_step robustness of the direct call

Running the direct call with `opt_step` ∈ {1,2,4,8,16,32} on ML-KEM-768:
all six runs converge to the identical result (ζ=32, log2(rop)=206.357).
**The direct call is fully robust to opt_step; the grid-resolution
problem is specific to the MATZOV cost model, not to dual_hybrid search
in general.**

### Experiment 5 — MATZOV ζ=0 anomaly in ML-KEM-512 / ML-KEM-1024

`LWE.dual_hybrid` (MATZOV) reports ζ=0 for both ML-KEM-512 and
ML-KEM-1024. Full enumeration (ζ=0..60, step=2) shows this is *not* the
true optimum:

| Parameter | MATZOV default | True global minimum (full scan) | Gap |
|---|---|---|---|
| ML-KEM-512  | ζ=0, log2(rop)=139.656  | ζ=16, log2(rop)=139.298 | 0.36 bits |
| ML-KEM-768  | ζ=20, log2(rop)=196.366 | ζ=24, log2(rop)=196.133 | 0.23 bits |
| ML-KEM-1024 | ζ=0, log2(rop)=262.336  | ζ=34, log2(rop)=261.257 | **1.08 bits** |

### Experiment 6 — Mechanism: discrete re-optimization of (p, t)

Fine-grained scans (step=1) of ζ, recording p (FFT modulus) and t
(k_fft) alongside cost, reveal the source of the non-convexity.

ML-KEM-1024 (ζ=0..40, step=1) shows a clear **sawtooth pattern**: t
takes only the values {120, 110, 100, 90, 80, 70, ...} because
`early_abort_range` steps t in units of 10. Within each t-block, cost
strictly decreases as ζ increases; at each block boundary (t drops by
10), cost jumps back up. Example:

```
zeta=8,  t=110, log2(rop)=261.840   <- local minimum within block
zeta=9,  t=110, log2(rop)=262.310   <- cost rises (still same block)
zeta=10, t=100, log2(rop)=263.103   <- block boundary, cost jumps up
...
zeta=34, t=80,  log2(rop)=261.257   <- true global minimum, several blocks later
```

The greedy search terminates at the first local rise (ζ=0→1 for
ML-KEM-1024, within the very first block) and never reaches the lower
minima in subsequent t-blocks.

### Conclusion

The ζ=20-vs-32 discrepancy originally reported is **not** a pure
grid-resolution artifact of `opt_step`. It stems from a difference in
cost model (FFT distinguisher assumed vs not). However, full enumeration
uncovered a **separate, confirmed issue**: MATZOV's hardcoded
`early_abort_range(step=10)` search over ζ and t produces a sawtooth
cost surface (caused by the 10-unit discretization of t), and the
greedy termination rule can stop at a local rise before reaching a lower
minimum in a later block. This is confirmed across all three ML-KEM
parameter sets, with the gap ranging 0.23–1.08 bits (worst case:
ML-KEM-1024). The direct `dual_hybrid` call (no FFT assumption) does not
exhibit this issue and is robust to `opt_step` across 1–32.

**Practical implication:** security-bit estimates produced by
`LWE.dual_hybrid` (MATZOV) for ML-KEM-1024-class parameters may
underestimate attacker cost by up to ~1 bit due to this grid artifact,
independent of the separate FFT-assumption question already discussed
in Section 5 / Section 11.

### Reproducibility

All experiments use lattice-estimator (main branch, as of 2026-07),
SageMath. Scripts:
- `verify_zeta_isolated.sage` — Experiment 1
- `verify_zeta_full_scan.sage` — Experiment 2
- `verify_zeta_direct_scan.sage` — Experiment 3
- `verify_opt_step_direct.sage` — Experiment 4
- `verify_zeta_512_1024.sage` — Experiment 5
- `verify_mechanism_512.sage`, `verify_mechanism_1024.sage` — Experiment 6

### Open follow-up

- Confirm whether NIST/official MATZOV-derived security-bit claims for
  ML-KEM-1024 are affected in practice (the 1.08-bit gap is within
  typical safety margins but should be checked against the specific
  claim being cited).
- Consider filing this as a GitHub issue against `malb/lattice-estimator`
  (see draft below).

---

## Draft GitHub issue (malb/lattice-estimator)

**Title:** MATZOV ζ/t search (`early_abort_range`, step=10) can miss the
global optimum due to sawtooth cost structure

**Body:**

`LWE.dual_hybrid` (which dispatches to `MATZOV.__call__` in
`lwe_dual.py`) searches ζ (k_enum) and t (k_fft) using
`early_abort_range(..., step=10)`, a greedy search that stops at the
first cost increase. Because t is only re-optimized in steps of 10, the
cost-vs-ζ curve has a sawtooth shape (cost decreases within a t-block,
then jumps up at each block boundary). If the first step from ζ=0
happens to land on a local rise, the search terminates immediately and
reports ζ=0, even when much lower-cost points exist at larger ζ in
later blocks.

Reproduced for ML-KEM-512, ML-KEM-768, and ML-KEM-1024 (full ζ
enumeration vs default `LWE.dual_hybrid` output):

| Parameter | Default (MATZOV) ζ | True optimum ζ (full scan) | log2(rop) gap |
|---|---|---|---|
| ML-KEM-512  | 0  | 16 | 0.36 |
| ML-KEM-768  | 20 | 24 | 0.23 |
| ML-KEM-1024 | 0  | 34 | 1.08 |

Minimal repro (ML-KEM-1024): full enumeration of ζ=0..40 (step=1) shows
t taking only values {120,110,...,70}, with cost rising at each block
boundary and a global minimum at ζ=34 that the default greedy search
never reaches because it stops at ζ=0→1.

Happy to share the full scan scripts/data if useful.

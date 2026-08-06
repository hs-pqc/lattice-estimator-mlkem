# From Empirical Formulas to Search Resolution: Rethinking ζ_optimal in Dual Hybrid Attacks on LWE

**Status:** Research Note (Work in Progress)
**Author:** Hangshin Cho (hs-pqc)
**Date:** June 2026 (initial formulas); July 2026 (reframing); July 2026 (search-resolution mechanism resolved, see Appendix)

---

## Abstract

We initially set out to find a closed-form empirical formula for the optimal splitting dimension ζ in dual hybrid attacks on LWE, as a function of lattice dimension n, modulus q, and secret distribution parameter η. Systematic experiments across 42 parameter sets (ML-KEM, ML-DSA, NTRU+, HAETAE) produced two candidate formulas with reasonable fit (max_err ≤ 5).

However, closer analysis of why the two cost-model call paths in lattice-estimator — `LWE.dual_hybrid` (indirect, MATZOV) and `dual_hybrid` (direct) — disagree on ζ by as much as 20–40 units for the same parameters revealed that treating ζ_optimal as a function of (n, q, η) alone is a misframing. ζ is determined jointly by (a) which call path is used, since the two paths encode different assumptions about whether an FFT distinguisher is realistically implementable, and (b) the optimizer's search structure — in particular, the indirect call's hardcoded internal grid step (10 units, independent of the public `opt_step` parameter), which is confirmed (see Appendix) to report ζ=0 for ML-KEM-512 and ML-KEM-1024 as a search-resolution artifact rather than a true optimum. A parameter-indexed formula cannot capture either effect, because it silently assumes the search that produced its training data was resolved finely enough to trust.

This reframes the useful question from "what is ζ_optimal(n, q, η)?" to "how does search resolution affect the reliability of the resulting rop (security-bit) estimate?" — a question about the estimator's behavior, not about the lattice parameters. Under this lens, the ~10-bit spread we observe in ML-KEM-768's security margin between the two call paths is better understood as an artifact of assumption/resolution mismatch than as a property of ML-KEM-768 itself. The original formulas are retained below (Section 5) as a documented first attempt, superseded by this reframing. The search-resolution question itself is fully resolved in the Appendix: the two call paths differ primarily because of the FFT-assumption (not grid resolution), but MATZOV's internal ζ/t search independently has a confirmed, quantified grid-resolution artifact (**0.60–1.19 bits** across ML-KEM-512/768/1024) with an identified mechanism.

---

## 1. Introduction

The dual hybrid attack combines lattice reduction (BKZ) with exhaustive search over a subset of secret coordinates. The splitting dimension ζ controls the trade-off between BKZ cost and guessing cost:

- **Large ζ:** fewer dimensions for BKZ → cheaper lattice reduction, but more secret coordinates to guess
- **Small ζ:** more dimensions for BKZ → expensive lattice reduction, but fewer coordinates to guess

The optimal ζ minimizes total attack cost. In practice, lattice-estimator finds ζ numerically via binary search (`local_minimum` in `lwe_dual.py`) for the direct call, and via a hardcoded greedy grid search (`early_abort_range`) for the MATZOV call path (see Appendix).

**Original Research Question (superseded):** Can ζ_optimal be expressed as a closed-form function of (n, q, η)?

**Current Research Question (resolved, see Appendix):** How does the resolution of the ζ search (grid step, call path, cost-model assumptions) affect the reliability of the resulting rop (security-bit) estimate? Put differently: when lattice-estimator reports a security level, how much of that number reflects the cryptographic hardness of the instance, and how much reflects the coarseness of the search that produced it?

---

## 2. Preliminaries

### 2.1 LWE Parameters

LWE instance: (A, b = As + e) where:

- n: lattice dimension
- q: modulus
- Xs = Xe = CenteredBinomial(η): secret/error distribution
  * Support: {-η, ..., η}
  * Variance: σ² = η/2
  * Standard deviation: σ = √(η/2)

### 2.2 Dual Hybrid Attack Cost

Total cost = BKZ cost + Guessing cost
BKZ cost ≈ 2^(0.292β)
Guessing cost ≈ (2η+1)^ζ
Optimal ζ satisfies: d/dζ [BKZ_cost(n-ζ) + Guessing_cost(ζ)] = 0

### 2.3 lattice-estimator

Open-source tool by Albrecht et al. Uses `local_minimum` binary search to find optimal (β, ζ) numerically.

**Key finding (revised, 2026-07):** MATZOV (indirect, `LWE.dual_hybrid`) and `dual_hybrid(fft=True)` are NOT the same technique despite both being labeled "FFT distinguisher." `dual_hybrid(fft=True)` uses the Walsh–Hadamard transform-based distinguisher from Guo-Johansson [AC:GuoJoh21, 2021], a distinct published algorithm from MATZOV's own FFT-based guessing cost model (2022).

For ML-KEM-768:

- MATZOV (indirect): log2(rop) = 196.37 (cheapest/strongest attack found)
- dual_hybrid(fft=True) (Guo-Johansson): 203.79 (ζ=31, t=120)
- dual_hybrid(fft=False): 206.36

MATZOV finds an attack **8.2 bits cheaper** than the Guo-Johansson FFT option. This corrects our earlier (incorrect) conclusion that indirect/MATZOV estimates were "overly optimistic" — in fact MATZOV identifies the stronger attack here, and any security-margin comparison must specify which of the two distinguishers is being used.

---

## 3. Experimental Setup

- Tool: SageMath + lattice-estimator 0.1.0
- Method: `dual_hybrid` direct call
- Parameters tested: 42 parameter sets across:
  * ML-KEM (512/768/1024)
  * ML-DSA (44/65/87)
  * NTRU+ (576/768/864/1152)
  * Custom parameters (various n, q, η)
  * FHE-type parameters (q=2^32, 2^60)
  * Extreme parameters (η=1, 8, 16; small q=1021)

---

## 4. Key Observations

### 4.1 Primal vs Dual

**Note:** The "Dual (bits)" column below uses `dual_hybrid(fft=True)` (Guo-Johansson distinguisher). Using MATZOV's own cost model instead yields dual estimates up to 8.2 bits lower (see §2.3) — i.e. this table likely understates the strength of the best-known dual attack.

| Standard   | Primal (bits) | Dual, Guo-Johansson FFT (bits) | Dual, MATZOV (bits) |
| ---------- | -------------- | -------------------------------- | ---------------------- |
| ML-KEM-768 | 204.9          | 206.4                            | 196.37 (default) / 195.55 (full ζ×t search) |

### 4.2 Cost Model Comparison: ζ values

Two cost models give significantly different ζ values:

| Parameter   | LWE.dual_hybrid (MATZOV) | dual_hybrid direct | Difference |
| ----------- | ------------------------- | -------------------- | ---------- |
| ML-KEM-512  | 0                          | 20                    | +20        |
| ML-KEM-768  | 20                         | 32                    | +12        |
| ML-KEM-1024 | 0                          | 41                    | +41        |
| ML-DSA-44   | 10                         | 23                    | +13        |
| ML-DSA-65   | 20                         | 25                    | +5         |
| ML-DSA-87   | 30                         | 44                    | +14        |
| NTRU+576    | 20                         | 25                    | +5         |
| NTRU+1152   | 30                         | 43                    | +13        |

MATZOV uses FFT distinguisher which makes ζ=0 optimal for some parameters. `dual_hybrid` direct call finds larger ζ values without FFT assumption. ML-KEM-512 and ML-KEM-1024 show ζ=0 under MATZOV — meaning no dimension reduction is optimal when FFT is available.

**Resolved (see Appendix, 2026-07):** the indirect call's ζ=0 output for ML-KEM-512/1024 is confirmed to be partly a search-resolution artifact. Full ζ×t enumeration shows the true MATZOV-model optimum is ζ=14 (ML-KEM-512, gap 0.60 bits), ζ=23 (ML-KEM-768, gap 0.81 bits), and ζ=32 (ML-KEM-1024, gap 1.19 bits).

### 4.3 ζ scales with n/log2(q)

ζ/n ≈ 0.030 + 0.096/√(n/log2(q))

This suggests ζ ≈ 0.030n + 0.096√(n × log2(q)) for small q.

### 4.4 C(η) follows η^(1/3)

When normalizing ζ × log2(q) / n, the resulting constant C varies with η following approximately η^(1/3):

| η | C = ζ×log2(q)×η^(1/3)/n |
| --- | ----------------------- |
| 1 | 0.625                   |
| 2 | 0.597                   |
| 3 | 0.587                   |
| 4 | 0.575                   |
| 8 | 0.561                   |

C converges to ≈ 0.597 for η ≥ 2.

---

## 5. Initial Formulas (Superseded) and Why the Approach Was Misframed

### 5.0 Why a Parameter-Indexed Formula Doesn't Work

The formulas below were fit to ζ values produced by the `dual_hybrid` direct call at a fixed `opt_step`. Two problems only became visible after the fact:

**Call-path dependence.** `LWE.dual_hybrid` (indirect, MATZOV) and `dual_hybrid` (direct) do not disagree because one is "more accurate" — they encode different assumptions about FFT distinguisher implementability (Section 4.2). A formula fit to one call path's output is a formula for that assumption, not for ζ_optimal in general. Fitting across both paths in one formula is not meaningful, since they are not estimating the same quantity.

**Search-resolution dependence.** The indirect call's hardcoded internal grid step (10 units) is confirmed (Appendix) to be coarse enough that for ML-KEM-512 and ML-KEM-1024 it does not evaluate the true optimum and reports ζ=0 by default instead of ζ=14 / ζ=32 respectively. Any formula trained on this output would inherit the artifact as if it were signal.

Because both effects are about the search and cost-model machinery, not about (n, q, η), no formula indexed only on parameter values can be correct in principle — it will always be a formula for "whatever the grid happened to find," which changes if the grid changes. This is the core reframing of this note: the two formulas below are kept as a documented first attempt, not as the intended contribution.

**Formula 1: Practical Formula**

> ζ ≈ floor((-0.076 + 0.701 × η^0.2) × n / (log2(q) × log2(η+1))) + 2

- Derived via least-squares fitting
- Valid: n≥512, η≥1, q=1021~2^60
- Accuracy: max_err=3, avg_err=1.00 (42 parameter sets)

**Formula 2: Mathematically Motivated**

> ζ ≈ floor(0.597 × n / (log2(q) × η^(1/3))) + 2

- Derived from minimizing variance of C(η) over α
- Optimal α = 1/3 (cube root)
- Valid: n≥512, η≥2, q=1021~2^60
- Accuracy: max_err=5, avg_err=1.29 (31 parameter sets)

---

## 6. Mathematical Intuition

### 6.1 Why η^(1/3)?

CenteredBinomial(η) has standard deviation σ = √(η/2).

η^(1/3) = (η/2)^(1/3) × 2^(1/3) ∝ σ^(2/3)

**Conjecture:** The guessing cost in dual hybrid attack scales as σ^(2/3), not σ or σ².

This is consistent with the balance condition:
0.292 × β(n-ζ) ≈ ζ × log2(guessing_cost_per_coordinate)
where guessing_cost_per_coordinate ∝ σ^(2/3) = η^(1/3).

### 6.2 Why C ≈ 0.597?

0.597 ≈ 0.292 × 2.044

0.292 is the BKZ exponent (Albrecht et al.). The factor 2.044 may relate to the ratio between lattice dimension and effective guessing dimension.

**Open Question:** What is the exact mathematical origin of C ≈ 0.597?

---

## 7. Open Problems

1. **[Primary — RESOLVED, see Appendix (2026-07)] Search resolution vs. rop reliability:** For a fixed instance, how does varying `opt_step` (grid coarseness) change the reported rop? Does the reported security level converge as `opt_step` → 1, or does it plateau early? If ML-KEM-512/1024's ζ=0 result is a resolution artifact, does refining the grid change their reported bit-security, or only ζ itself?
   → **Answer:** For the direct call, rop is fully converged and `opt_step`-independent (tested `opt_step` ∈ {1,2,4,8,16,32}, identical result each time). For the MATZOV call path, the grid artifact is real. Mechanism identified: MATZOV's internal t (FFT dimension) search is hardcoded to step in units of 10, producing a sawtooth cost surface. Full ζ×t resolution (Experiment 7) gives confirmed gaps of 0.60 bits (ML-KEM-512), 0.81 bits (ML-KEM-768), and 1.19 bits (ML-KEM-1024).

2. **Mathematical proof of η^(1/3):** Prove that guessing cost in dual hybrid scales as σ^(2/3). Direct verification shows LHS/RHS ≈ 2~3 (not 1), suggesting the balance condition is more complex than BKZ_cost = Guess_cost. (Still open.)

3. **Origin of C ≈ 0.597:** 0.597 ≈ 0.292 × 2.044. The factor 2.044 may relate to repetition cost, sieving dimension, or memory cost in BKZ. Deriving C analytically from the full cost function remains open. (Still open.)

4. **η=1 anomaly:** Formula 2 has larger errors for η=1 (max_err=5). Empirically, ζ(η=1)/ζ(η=2) ≈ 1.38 consistently across n. This ratio is larger than η^(1/3) predicts (2^(1/3) ≈ 1.26). Conjecture: η=1 sparsity introduces an additional factor of (4/3)^(1/3). (Still open.)

5. **Cost model independence:** All experiments use MATZOV cost model (default in lattice-estimator 0.1.0). ADPS16, Kyber, and quantum cost models are not compatible with `dual_hybrid` in this version. Verification under other cost models remains open. (Still open.)

6. **SMAUG-T analysis:** SMAUG-T uses sparse secret (HWT distribution) with small h. lattice-estimator 0.1.0 does not support SparseTernary(n, h) for small h values. Separate analysis needed. (Still open — see also Section 10 connection to approximate hints.)

7. **Connection to security margins:** How does ζ_optimal relate to actual security margins in NIST/KPQC PQC standards? (Still open.)

8. **Generalization beyond CBD:** All experiments use CenteredBinomial(η). Does the formula extend to DiscreteGaussian or other distributions? (Still open.)

---

## 8. Experimental Data

All experimental results are in the `results/` directory.

**Core results (results/1_core/)**

| File | Content |
| --- | --- |
| all_security_revised.txt | Primal/dual comparison (all standards) |
| zeta_accurate.txt | MATZOV vs direct call ζ comparison |
| zeta_final.txt | Formula 1 blind test |
| zeta_verify3.txt | Formula 1 validation (17 parameter sets) |
| zeta_cube_verify.txt | Formula 2 validation (31 parameter sets) |
| kpqc_test.txt | HAETAE validation |
| hints_threshold.txt | ML-KEM-512 hint threshold analysis |
| hints_all.txt | ML-KEM hint threshold comparison |
| C_vs_nlogq.png | C vs n/log2(q) visualization |

**Search-resolution verification (results/1_core/, added 2026-07)**

| File | Content |
| --- | --- |
| verify_zeta_isolated.sage | ζ-only grid isolation (ML-KEM-768) |
| verify_zeta_full_scan.sage | Full ζ enumeration under MATZOV (ML-KEM-768) |
| verify_zeta_direct_scan.sage | Full ζ enumeration, direct call (reproduces ζ=32) |
| verify_opt_step_direct.sage | opt_step robustness test, direct call |
| verify_zeta_512_1024.sage | Full ζ enumeration under MATZOV (ML-KEM-512/1024) |
| verify_mechanism_512.sage | p/t trace, ML-KEM-512 |
| verify_mechanism_1024.sage | p/t trace, ML-KEM-1024 (sawtooth mechanism) |

**Analysis (results/2_analysis/)**

Intermediate experiments and fitting results.

**Deprecated (results/3_deprecated/)**

Early results using `LWE.dual_hybrid` (MATZOV) — kept for reference.

---

## 9. References

1. Albrecht et al., "On the Concrete Hardness of LWE", JMC 2015
2. MATZOV, "Report on the Security of LWE", 2022
3. Espitau, Joux, Kharchenko, "On a Dual/Hybrid Approach", INDOCRYPT 2020
4. NIST FIPS 203 (ML-KEM), 2024
5. NIST FIPS 204 (ML-DSA), 2024
6. lattice-estimator, github.com/malb/lattice-estimator
7. Guo, Johansson, "A New Sieving-Style Information-Set Decoding Algorithm", ASIACRYPT 2021
8. Hhan, Hong, Kim, Lee, Lee, "From Perfect to Approximate Hints: Efficient LWE Secret Recovery Leveraging Low Hamming Weight", S&P 2026 (ePrint 2026/1081) — co-authored by Changmin Lee
9. Kim, Lee, Kim, Lee, "SQIsign with Fixed-Precision Integer Arithmetic", ePrint 2025/1649 — Changmin Lee, co-author (first paper in this line)
10. Kim, Lee, Yoo, "Compact Quaternion Algorithms for SQIsign", ePrint 2026/1031 — Changmin Lee, co-author (follow-up, 74% precision reduction over [9])

---

## 10. Connection to Approximate Hints

**Note on methodology:** The hint model below (reducing effective dimension n → n−h) is a rough approximation motivated by the existence of hints, not a reproduction of Hhan et al.'s actual algorithm — their Algorithm 3 is a correlation-based method that does not use lattice reduction at all. A more accurate analysis would require the DBDD framework.

### Motivation

If FFT distinguisher (MATZOV) is realistic, ML-KEM-512 has only 11.7-bit security margin above 128 bits. Combined with approximate hints (Hhan et al., S&P 2026), this margin may shrink further.

### Experiment: hint threshold analysis

We approximate hint effects by reducing effective dimension n → n-h. This is a rough model; accurate analysis requires DBDD framework.

| Parameter | Base security | Margin | Threshold |
| --- | --- | --- | --- |
| ML-KEM-512 | 139.7 bits | 11.7 bits | ~44 hints |
| ML-KEM-768 | 196.4 bits | 4.4 bits | ~16 hints |
| ML-KEM-1024 | 262.3 bits | 6.3 bits | ~21 hints |

---

## 11. FFT Distinguisher and Security Margin Uncertainty

### Key Claim

The practical security of ML-KEM-768 critically depends on whether the FFT distinguisher (MATZOV) is realistically implementable.

| Parameter | Target | FFT available | FFT unavailable | Uncertainty |
| --- | --- | --- | --- | --- |
| ML-KEM-512 | 128 bits | 11.7 bits | 17.5 bits | 5.9 bits |
| ML-KEM-768 | 192 bits | 4.4 bits | 14.4 bits | 10.0 bits |
| ML-KEM-1024 | 256 bits | 6.3 bits | 21.5 bits | 15.2 bits |

### Interpretation

- If FFT is realistic → ML-KEM-768 margin is only 4.4 bits (most vulnerable)
- If FFT is unrealistic → ML-KEM-768 margin is 14.4 bits (safe)
- ML-KEM-768 has the largest security uncertainty (10 bits) among all three
- Of this uncertainty, ≤0.81 bits (ML-KEM-768) is attributable to the MATZOV grid-resolution artifact (Appendix); the remainder is a genuine open question about FFT distinguisher realizability — the estimator artifact contributes only a small fraction of that spread.

### Open Question

What is the realistic implementation cost of the FFT distinguisher? This single question determines whether ML-KEM-768 is the weakest or the middle parameter set among the three ML-KEM variants.

Combined with approximate hints (Section 10), if FFT is realistic and 16 hints are obtainable, ML-KEM-768 NIST Level 3 security may be at risk.

### Key Finding

ML-KEM-768 has the smallest margin (4.4 bits) among all three parameter sets. Only ~16 approximate hints are needed to break NIST Level 3 security (192 bits).

### Open Question

Can an attacker realistically obtain 16 approximate hints against ML-KEM-768 in practice? If so, NIST Level 3 security may be at risk under MATZOV + approximate hints combined attack.

This connects directly to Hhan et al. S&P 2026 which analyzes the cost reduction from approximate vs perfect hints.

---

## Appendix: MATZOV ζ Search Resolution — Mechanism and Verification (2026-07)

### Background

`LWE.dual_hybrid` in lattice-estimator resolves to `MATZOV.__call__` (`estimator/lwe_dual.py`), which searches over ζ (k_enum) and t (k_fft) using `early_abort_range(..., step=10)` — a hardcoded grid step, independent of the public `opt_step` parameter. `early_abort_range` is a greedy search: it terminates as soon as the cost increases once. This appendix documents a full reconstruction of the ζ=20 vs ζ=32 discrepancy originally reported in Section 4, and identifies the exact mechanism behind it.

**Experiment 1 — Isolating the MATZOV ζ search**

For ML-KEM-768, the default MATZOV call finds ζ=20 (log2(rop)=196.4). Adjusting the internal grid step for ζ alone (with t held at the default step of 10) still gives ζ=4 — worse, not better, than the default. This rules out simple "finer grid → better answer" as the story; the search is path-dependent, not merely under-resolved.

**Experiment 2 — Full ζ enumeration under MATZOV (ML-KEM-768)**

Scanning ζ=0..768 (step=4), with p and t re-optimized at each ζ:
- Local non-convexity for ζ=0..32 (amplitude ≈1.7 bits)
- Monotonic increase beyond ζ≈36
- True global minimum: ζ=24, log2(rop)=196.133
- Default MATZOV output (ζ=20) is only 0.23 bits above this local scan's optimum (superseded by Experiment 7's fully-resolved figure below)

**Experiment 3 — Reproducing the original ζ=32 (direct call, no FFT)**

The original ζ=20 vs 32 comparison in Section 4 was between MATZOV (FFT distinguisher assumed) and the direct `dual_hybrid` call (FFT distinguisher not assumed, `fft=False`), not between two grid resolutions. A full enumeration of the direct call reproduces the recorded optimum exactly:

Global minimum: ζ=32, log2(rop)=206.357 — matches Section 4's recorded value (ζ=32, rop=206.4) exactly.

The direct-call cost curve is smooth and strictly unimodal (monotonic decrease then increase around ζ=32), unlike the MATZOV curve.

**Experiment 4 — opt_step robustness of the direct call**

Running the direct call with `opt_step` ∈ {1,2,4,8,16,32} on ML-KEM-768: all six runs converge to the identical result (ζ=32, log2(rop)=206.357). The direct call is fully robust to `opt_step`; the grid-resolution problem is specific to the MATZOV cost model, not to `dual_hybrid` search in general.

**Experiment 5 — MATZOV ζ=0 anomaly in ML-KEM-512 / ML-KEM-1024**

`LWE.dual_hybrid` (MATZOV) reports ζ=0 for ML-KEM-512 and ML-KEM-1024, and ζ=20 for ML-KEM-768. An initial full enumeration of ζ alone (step=2, with t still searched via the default `early_abort_range(step=10)`) found gaps of 0.36 / 0.23 / 1.08 bits. However, this scan inherited the same coarse t-resolution the issue is about — see Experiment 7 below for the corrected, fully-resolved figures.

**Experiment 6 — Mechanism: discrete re-optimization of (p, t)**

Fine-grained scans (step=1) of ζ, recording p (FFT modulus) and t (k_fft) alongside cost, reveal the source of the non-convexity.

ML-KEM-1024 (ζ=0..40, step=1) shows a clear sawtooth pattern: t takes only the values {120, 110, 100, 90, 80, 70, ...} because `early_abort_range` steps t in units of 10. Within each t-block, cost strictly decreases as ζ increases; at each block boundary (t drops by 10), cost jumps back up. Example:

```
zeta=8,  t=110, log2(rop)=261.840   <- local minimum within block
zeta=9,  t=110, log2(rop)=262.310   <- cost rises (still same block)
zeta=10, t=100, log2(rop)=263.103   <- block boundary, cost jumps up
...
zeta=34, t=80,  log2(rop)=261.257   <- lower minimum, several blocks later (superseded by Experiment 7's global figure below)
```

The greedy search terminates at the first local rise (ζ=0→1 for ML-KEM-1024, within the very first block) and never reaches the lower minima in subsequent t-blocks.

**Experiment 7 — Full resolution in both ζ and t (2026-07, added)**

The scan in Experiment 5 fixed ζ exhaustively but still searched t via `early_abort_range(step=10)` internally, so it inherited the same sawtooth artifact it was trying to measure. Re-running with both ζ and t scanned at step=1 (SageMath, lattice-estimator main branch, 2026-07) gives:

| Parameter | MATZOV default | True global optimum (ζ, t step=1) | Gap |
| --- | --- | --- | --- |
| ML-KEM-512 | ζ=0, log2(rop)=139.656 | ζ=14, log2(rop)=139.057, t=34, β=385 | 0.599 bits |
| ML-KEM-768 | ζ=20, log2(rop)=196.366 | ζ=23, log2(rop)=195.554, t=59, β=586 | 0.812 bits |
| ML-KEM-1024 | ζ=0, log2(rop)=262.336 | ζ=32, log2(rop)=261.143, t=82, β=819 | 1.192 bits |

These are the final, fully-resolved gaps — larger than the preliminary step=10-in-t figures from Experiment 5 (0.36 / 0.23 / 1.08), which themselves understated the artifact. Notably, for ML-KEM-1024 the gap grew even at ζ=0 itself: full t-resolution at ζ=0 gives log2(rop)=262.054, 0.28 bits below the tool's own default output — confirming that the coarse t grid distorts the cost surface even before ζ is varied at all.

Total scan time (full ζ×t×p enumeration, ternary-search β at each point): 512 — 6201s, 768 — 7001s, 1024 — 11658s (≈6.9h combined), illustrating why `early_abort_range(step=10)` is likely a deliberate speed/precision tradeoff rather than an oversight.

**Experiment 8 — MATZOV vs. Guo-Johansson(2021) FFT paths are different attacks (2026-07)**

Checked whether `LWE.dual_hybrid` (MATZOV) and `dual_hybrid(fft=True)` — both described as using an "FFT distinguisher" — actually implement the same technique. They do not: `dual_hybrid(fft=True)` uses the Walsh–Hadamard transform-based distinguisher from [AC:GuoJoh21], a different published algorithm from MATZOV's own FFT-based guessing cost model.

For ML-KEM-768 (SageMath, lattice-estimator main):

| Method | log2(rop) |
| --- | --- |
| MATZOV (`LWE.dual_hybrid`), default | 196.37 |
| `dual_hybrid(fft=True)` (Guo-Johansson) | 203.79 (ζ=31, t=120) |
| `dual_hybrid(fft=False)` | 206.36 |

MATZOV finds an attack 8.2 bits cheaper than the generic framework's FFT option, despite both nominally "using FFT." This confirms that the ζ discrepancy discussed in Experiments 1–4 is not just a resolution artifact but reflects genuinely different attack algorithms under the hood — and that MATZOV, not `dual_hybrid(fft=True)`, is the correct baseline for the sawtooth analysis in Experiment 7.

### Conclusion (updated)

The MATZOV ζ/t search's hardcoded `early_abort_range(step=10)` produces a sawtooth cost surface in both ζ and t. Confirmed, fully-resolved gaps between the tool's default output and the true global optimum: **0.60 bits (ML-KEM-512), 0.81 bits (ML-KEM-768), 1.19 bits (ML-KEM-1024)**. This supersedes the preliminary 0.36/0.23/1.08-bit figures, which themselves used a step=10 t-search and therefore understated the artifact. The direct `dual_hybrid` call (no FFT assumption) remains unaffected and robust to `opt_step` (Experiment 4).

**Practical implication:** security-bit estimates produced by `LWE.dual_hybrid` (MATZOV) for ML-KEM-1024-class parameters may underestimate attacker cost by up to ~1.2 bits due to this grid artifact, independent of the separate FFT-assumption question already discussed in Section 5 / Section 11.

### Reproducibility

All experiments use lattice-estimator (main branch, as of 2026-07), SageMath. Scripts:

- `verify_zeta_isolated.sage` — Experiment 1
- `verify_zeta_full_scan.sage` — Experiment 2
- `verify_zeta_direct_scan.sage` — Experiment 3
- `verify_opt_step_direct.sage` — Experiment 4
- `verify_zeta_512_1024.sage` — Experiment 5
- `verify_mechanism_512.sage`, `verify_mechanism_1024.sage` — Experiment 6

### Open follow-up

- Confirm whether NIST/official MATZOV-derived security-bit claims for ML-KEM-1024 are affected in practice (the 1.19-bit gap is within typical safety margins but should be checked against the specific claim being cited).
- Consider filing this as a GitHub issue against malb/lattice-estimator (see draft below).

### Draft GitHub issue (malb/lattice-estimator)

**Title:** MATZOV ζ/t search (early_abort_range, step=10) can miss the global optimum due to sawtooth cost structure

**Body:**

`LWE.dual_hybrid` (which dispatches to `MATZOV.__call__` in `lwe_dual.py`) searches ζ (k_enum) and t (k_fft) using `early_abort_range(..., step=10)`, a greedy search that stops at the first cost increase. Because t is only re-optimized in steps of 10, the cost-vs-ζ curve has a sawtooth shape (cost decreases within a t-block, then jumps up at each block boundary). If the first step from ζ=0 happens to land on a local rise, the search terminates immediately and reports ζ=0, even when much lower-cost points exist at larger ζ in later blocks.

Reproduced for ML-KEM-512, ML-KEM-768, and ML-KEM-1024 (full ζ×t enumeration vs default `LWE.dual_hybrid` output):

| Parameter   | Default (MATZOV) ζ | True optimum ζ (full ζ×t scan) | log2(rop) gap |
|-------------|---------------------|----------------------------------|----------------|
| ML-KEM-512  | 0                   | 14 (t=34, β=385)                 | 0.60           |
| ML-KEM-768  | 20                  | 23 (t=59, β=586)                 | 0.81           |
| ML-KEM-1024 | 0                   | 32 (t=82, β=819)                 | 1.19           |

Minimal repro (ML-KEM-1024): full enumeration of ζ and t (step=1) shows the tool's internal `early_abort_range(step=10)` restricts t to values {120,110,...,70}, producing a sawtooth cost surface; the greedy search terminates at ζ=0→1 and never reaches the true global minimum at ζ=32.

Happy to share the full scan scripts/data if useful.
## Proposed Fix (2026-08, added)

### Why full step=1 enumeration is not a viable permanent fix
Experiment 7's fully-resolved scan took ~6.9h combined for three parameter
sets. This rules out full enumeration as a default estimator behavior — it
would make routine security-margin verification (e.g. at CRYPTREC scale,
across dozens of parameter sets) impractical.

### Design iterations (three attempts before a working design was found)

**Attempt 1 — naive coarse-to-fine on zeta alone, t fully resolved.**
Non-greedy coarse scan over zeta (step=10, fast t-search) to locate a
candidate region, then a fine window around it with t resolved by full
step=1 enumeration (as in Experiment 7). Correct (ML-KEM-512: zeta=14,
t=34, log2(rop)=139.057, exact match), but only ~2x faster than full
enumeration (3159.9s vs 6201s) — the fine window still inherits the large
t_max near zeta=0, so most of the savings from skipping most of the zeta
range were eaten by the remaining exhaustive t-scans.

**Attempt 2 — greedy zeta search with fully-resolved t (t_step=1 via
`early_abort_range` instead of the internal step=10).** Hypothesis: if the
t-sawtooth (Experiment 6) is the only source of non-convexity, resolving t
more finely should make the zeta cost curve smooth enough for a plain
greedy search. Result: **falsified**. The greedy zeta search still
terminated at zeta=0/1 (26.5s, `MISMATCH` — true optimum is zeta=14). The
zeta cost curve has genuine local non-convexity independent of t
resolution (consistent with the ~1.7-bit amplitude non-convexity already
noted in Experiment 2) — a non-greedy exhaustive pass over zeta candidates
is unavoidable.

**Attempt 3 — greedy *t* search (via `early_abort_range`, step=1) inside
an otherwise-correct non-greedy zeta scan.** Hypothesis: t itself, for a
*fixed* zeta, might be unimodal even though zeta is not. Result: **also
falsified**. At zeta=14, greedy t-search reported log2(rop)=139.554 versus
the true 139.057 (0.5-bit error) — t has local non-convexity too, at fixed
zeta. Exhaustive t search cannot be replaced by a greedy search at any
level.

### Final design: recursive non-greedy coarse-to-fine (zeta and t), parallelized

Since neither zeta nor t admits a greedy search, the same "non-greedy
coarse pass + boundary-expanding fine pass" pattern proposed for zeta in
the original draft of this section is applied **twice, recursively**:

1. **Zeta, coarse:** scan zeta over `[0, n)` at step=10, evaluating every
   point (no early-abort) using the fast per-zeta cost (t via
   `early_abort_range`, step=10) — this stage does not need to be exact,
   only good enough to locate the right neighborhood.
2. **Zeta, fine:** for each zeta in a window around the coarse candidate
   (initial width 15, doubling on boundary hits, capped at 3 doublings),
   compute an accurate per-zeta cost using:
3. **t, coarse:** scan t over `[0, n-zeta)` at step=10, evaluating every
   point (no early-abort).
4. **t, fine:** window (width 15, doubling on boundary hits, capped at 3
   doublings) around the t coarse candidate, step=1.
5. **Parallelization:** the zeta-fine stage's per-zeta evaluations are
   independent and are distributed across a `multiprocessing.Pool`. No
   algorithmic change, no accuracy cost — this only reduces wall-clock
   time.

Steps 3–4 mirror steps 1–2 exactly, one level down. This is the same
pattern already validated for zeta, now shown to be necessary (not
optional) for t as well.

### Validation results

All three ML-KEM parameter sets reproduce Experiment 7's confirmed global
optima exactly (zeta, t, and log2(rop) match to the values reported
there):

| Parameter    | zeta | t  | log2(rop) | Match | Time (this design) | Time (Exp. 7, full enum.) | Speedup |
| ------------ | ---- | -- | --------- | ----- | ------------------- | -------------------------- | ------- |
| ML-KEM-512   | 14   | 34 | 139.057   | OK    | 97.2s               | ~6201s                     | ~64x    |
| ML-KEM-768   | 23   | 59 | 195.554   | OK    | 123.0s              | ~7001s                     | ~57x    |
| ML-KEM-1024  | 32   | 82 | 261.143   | OK    | 155.7s              | ~11658s                    | ~75x    |
| **Combined** |      |    |           | **OK**| **376.9s (6m17s)**  | **~24840s (6.9h)**         | **~66x**|

Exact match to Experiment 7 on all three parameter sets, at roughly 1/66th
of the wall-clock cost, on a 32-core machine.

### Complexity comparison (revised)

The original draft of this section (Section 5.0 as originally written)
claimed an "order-of-magnitude" speedup from a single-level coarse+fine
design based on an asymptotic argument alone; that claim did not survive
contact with measurement (Attempt 1 above measured only ~2x). The
recursive two-level design plus parallelization is what actually delivers
the improvement, and only the measured numbers above should be cited going
forward — not the earlier asymptotic estimate.

### Reproducibility
All experiments in this section use lattice-estimator (main branch, as of
2026-08), SageMath, run on a 32-core container. The script auto-detects
available cores via `os.cpu_count()` rather than hardcoding a core count,
so wall-clock time will scale with whatever hardware it runs on. Script:
`results/1_core/verify_parallel_v2.sage` (contains both the fast per-zeta
cost function used for coarse passes and the recursive two-level t-search
used for fine passes; parallelizes the zeta-fine stage via
`multiprocessing.Pool`).

Earlier falsified variants are kept for the record:
`results/1_core/verify_fine_t_hypothesis.sage` (Attempt 2 — greedy zeta,
falsified) and `results/1_core/test_t_two_level.sage` (single-zeta
validation of the working two-level t-search, used to confirm Attempt 3's
failure and motivate the final design).

### Open follow-up
- Confirm this design scales similarly for ML-DSA and NTRU+ parameter
  sets (Section 3's full 42-set list), not just the three ML-KEM sizes
  tested here.
- The initial window (15) and coarse step (10) were chosen to match the
  t-sawtooth period identified in Experiment 6; a smaller coarse step
  might further reduce the (small) residual boundary-doubling overhead
  seen in the ML-KEM-1024 case.
- Consider filing the recursive coarse-to-fine design, alongside the
  original bug report, as a proposed patch to `malb/lattice-estimator`'s
  `MATZOV.__call__`.
## Bug Scope Beyond ML-KEM (2026-08, added)

The Proposed Fix section above validates the recursive coarse-to-fine
design on ML-KEM only. Two natural questions follow: does the underlying
bug (Section 5.0 / Appendix) affect other lattice-based schemes, and does
the fix generalize? Both were tested directly rather than assumed.

### NTRU is not affected

`NTRU.dual_hybrid` and `LWE.dual_hybrid` (the MATZOV instance carrying the
`early_abort_range(step=10)` bug) are **different function objects**
(`NTRU.dual_hybrid is LWE.dual_hybrid` → `False`). Inspecting
`NTRU.dual_hybrid`'s source shows it dispatches to a generic `DH`
optimizer controlled by `opt_step` (default 8), not the MATZOV-specific
greedy ζ/t search. NTRU parameter sets (NTRUHPS/HRSS) were therefore
**not** re-tested with the fix, since there is no bug at that code path
to fix.

### ML-DSA is affected

ML-DSA's MLWE key-recovery parameters are not built into
`lattice-estimator`'s `schemes` module (only `_MSIS_*` forgery parameters
are). Constructed manually per FIPS 204 (k, l, η) triples, following the
standard convention n = l·256, m = k·256, q = 8380417,
Xs = Xe = `ND.UniformMod(η)`:

| Parameter | Default ζ | True ζ | Default log2(rop) | True log2(rop) | Gap (bits) |
| --------- | --------- | ------ | ------------------ | -------------- | ---------- |
| ML-DSA-44 | 10        | 23     | 126.682             | 126.337         | 0.345      |
| ML-DSA-65 | 20        | 21     | 180.456             | 179.992         | 0.464      |
| ML-DSA-87 | 0         | 13     | 225.521             | 225.255         | 0.266      |

All three show a real gap in the same direction as ML-KEM. Notably,
ML-DSA-65's ζ changes by only 1 (20→21) yet has the largest gap of the
three (0.464 bits) — the artifact is driven by t's re-optimization
granularity at least as much as by which ζ is selected.

Script: `results/1_core/verify_dsa_all.sage`.

### NTRU+ (KpqC) is affected — unlike NIST NTRU

NTRU+ differs from NIST NTRU in a way that matters here: its security
paper states that RLWE-based security is estimated using
`lattice-estimator`, in addition to a separate NTRU-problem estimate.
Constructing NTRU+576 (n=576, q=3457, approximated as binary secret +
sparse ternary error) and calling `LWE.dual_hybrid` directly confirms it
routes through the buggy MATZOV path:

| Parameter | Default ζ | True ζ | Default log2(rop) | True log2(rop) | Gap (bits) |
| --------- | --------- | ------ | ------------------ | --------------- | ---------- |
| NTRU+576  | 20        | 26     | 130.716             | 130.588          | 0.128      |

(Distribution is an approximation, not a certified match to the NTRU+
spec's exact sampling — this establishes the code path is affected, not
a certified security number.) Script: `results/1_core/test_ntruplus.sage`.

### FrodoKEM — largest gap found in the entire study

FrodoKEM is built into `lattice-estimator`'s `schemes` module
(`Frodo640`, `Frodo976`, `Frodo1344`), so these results carry the same
confidence as the ML-KEM results — no manual parameter construction.

| Parameter    | Default ζ, t | Default log2(rop) | True ζ, t | True log2(rop) | Gap (bits) |
| ------------ | ------------ | ------------------- | --------- | --------------- | ---------- |
| Frodo640     | 10, 0        | 170.146              | 13, 0     | 169.326          | 0.821      |
| **Frodo976** | **20, 0**    | **231.495**          | **17, 43**| **223.513**      | **7.981**  |
| Frodo1344    | 30, 60       | 281.739              | 28, 69    | 281.271          | 0.468      |

Frodo976 is qualitatively different from every other case in this study:
the default search concludes t=0 is optimal (no benefit from the FFT
distinguisher), while the true optimum uses t=43 and is 7.981 bits
cheaper (~253x in cost). This is a wrong conclusion about which attack
technique applies, not just an imprecise number.

Script: `results/1_core/verify_frodo_all.sage`.

### Updated conclusion
The bug affects every `LWEParameters`-based scheme tested that routes
through `LWE.dual_hybrid` (ML-KEM, ML-DSA, NTRU+, FrodoKEM), with gaps
ranging from 0.128 to **7.981** bits depending on parameters. It does not
affect NTRU (HPS/HRSS), whose dual-hybrid estimation bypasses MATZOV via
a separate function (`NTRU.dual_hybrid`).

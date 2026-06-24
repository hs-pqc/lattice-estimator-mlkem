# Empirical Analysis of Optimal Splitting Dimension ζ in Dual Hybrid Attacks on LWE

**Status:** Research Note (Work in Progress)  
**Author:** Hangshin Jo (hs-pqc)  
**Date:** June 2026

---

## Abstract

## Key Finding

### 1. Cost Model Comparison: LWE.dual_hybrid vs dual_hybrid direct call

Two methods give different results due to different underlying assumptions:

| Method | Model | ζ (ML-KEM-512) | rop (ML-KEM-512) |
|---|---|---|---|
| LWE.dual_hybrid | MATZOV (official standard) | 0 | 139.7 bits |
| dual_hybrid direct | Conservative model | 20 | 145.5 bits |

**LWE.dual_hybrid (MATZOV)** is the standard used in official security
analyses (CRYPTREC, NIST submissions, MATZOV report).
It gives lower security estimates — conservative and safe for parameter design.

**dual_hybrid direct call** uses a different cost model that finds
larger ζ values and gives higher security estimates.

Neither is "wrong" — they reflect different assumptions about attack costs.

### 2. Primal vs Dual under conservative cost model

Under dual_hybrid direct call, primal attack is more efficient
than dual hybrid across all standard parameters:

| Standard | Primal (bits) | Dual (bits) | Primal advantage |
|---|---|---|---|
| ML-KEM-512 | 143.8 | 145.5 | +1.7 |
| ML-KEM-768 | 204.9 | 206.4 | +1.5 |
| ML-KEM-1024 | 275.1 | 277.5 | +2.4 |
| HAETAE-2 | 102.5 | 102.9 | +0.4 |

HAETAE-2 has the smallest gap — dual hybrid is closest to primal
at small n with large q.

### 3. ζ_optimal empirical formula

Despite cost model differences, ζ_optimal follows a predictable
pattern — suggesting the formula reflects an intrinsic structure
of the dual hybrid optimization landscape, independent of cost model choice.

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

**Research Question:**  
Can ζ_optimal be expressed as a closed-form function of (n, q, η)?

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

### 4.2 ζ scales with n/log2(q)

ζ/n ≈ 0.030 + 0.096/√(n/log2(q))

This suggests ζ ≈ 0.030n + 0.096√(n × log2(q)) for small q.

### 4.3 C(η) follows η^(1/3)

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

## 5. Proposed Formulas

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

Key files:
- `all_security_revised.txt`: Primal/dual comparison
- `zeta_verify3.txt`: Formula 1 validation (17 parameter sets)
- `zeta_final.txt`: Formula 1 blind test
- `zeta_cube_verify.txt`: Formula 2 validation (31 parameter sets)
- `C_vs_nlogq.png`: C vs n/log2(q) visualization

---

## 9. References

1. Albrecht et al., "On the Concrete Hardness of LWE", JMC 2015
2. MATZOV, "Report on the Security of LWE", 2022
3. Espitau, Joux, Kharchenko, "On a Dual/Hybrid Approach", INDOCRYPT 2020
4. NIST FIPS 203 (ML-KEM), 2024
5. NIST FIPS 204 (ML-DSA), 2024
6. lattice-estimator, github.com/malb/lattice-estimator
7. Cheon et al., "HAETAE: Shorter Lattice-Based Fiat-Shamir Signatures", CHES 2024
8. Cheon et al., "SMAUG-T: Post-Quantum KEM", KpqC 2025

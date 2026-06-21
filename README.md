# ML-KEM Security Analysis using lattice-estimator

Security bit estimation for ML-KEM parameter sets using the
[lattice-estimator](https://github.com/malb/lattice-estimator).

## Environment

- SageMath 9.x
- lattice-estimator
- Docker

## Parameters

```python
# ML-KEM-512
LWE.Parameters(n=512, q=3329, Xs=ND.CenteredBinomial(3), Xe=ND.CenteredBinomial(3))

# ML-KEM-768
LWE.Parameters(n=768, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2))

# ML-KEM-1024
LWE.Parameters(n=1024, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2))
```

## Results

### Primal Attack (usvp)

| Parameter | n | Security (bits) | β (BKZ block size) |
|-----------|---|-----------------|-------------------|
| ML-KEM-512  | 512  | 143.8 | 406 |
| ML-KEM-768  | 768  | 204.9 | 624 |
| ML-KEM-1024 | 1024 | 275.1 | 874 |

### Dual Hybrid Attack

| Parameter | Security (bits) | β | p | ζ |
|-----------|----------------|---|---|---|
| ML-KEM-512  | 139.7 | 387 | 5 | 0  |
| ML-KEM-768  | 196.4 | 589 | 4 | 20 |
| ML-KEM-1024 | 262.3 | 823 | 4 | 0  |

## Key Observations

**1. Dual hybrid is consistently more efficient than primal**
ML-KEM-512:  Δ = 4.1 bits

ML-KEM-768:  Δ = 8.5 bits

ML-KEM-1024: Δ = 12.8 bits
The gap grows as n increases.

**2. Security scaling**
512  → 768:  +61.1 bits (primal), +56.7 bits (dual)

768  → 1024: +70.2 bits (primal), +65.9 bits (dual)
Not exactly 64 bits per step — the increment grows with n.

**3. ML-KEM-768 anomaly (ζ=20)**  
ζ represents additional dimension reduction beyond the guessing part.  
The optimizer finds different strategies per parameter:

- ML-KEM-512 (η=3): Secret coefficients are relatively large →  
  guessing (p=5) is more efficient than dimension reduction (ζ=0)
- ML-KEM-768 (η=2): Balanced n and small secret →  
  dimension reduction (ζ=20) hits the optimal cost tradeoff
- ML-KEM-1024 (η=2): n is too large →  
  dimension reduction gains are outweighed by guessing cost (ζ=0)

This suggests the dual hybrid optimizer is sensitive to the
interaction between secret distribution width (η) and lattice dimension (n).
## Sparse Secret Analysis

Effect of Hamming weight on security bits (ML-KEM-768, dual hybrid).

| Hamming weight (h) | Security (bits) | Δ from standard | ζ |
|-------------------|----------------|-----------------|---|
| Standard (CBD η=2) | 196.4 | - | 20 |
| h=30 | 164.1 | -32.3 | 0 |
| h=50 | 169.9 | -26.5 | 0 |
| h=100 | 181.1 | -15.3 | 30 |
| h=200 | 190.9 | -5.5 | 20 |

**Key Observation:**  
Security bits drop sharply as Hamming weight decreases.  
h=30 shows the largest gap (-32.3 bits from standard).  
As h increases, security converges back toward the standard parameter.

This directly relates to the key finding of:  
*"From Perfect to Approximate Hints"* (S&P 2026, Changmin Lee et al.)  
— sparse secrets amplify the effectiveness of hints in LWE attacks.

## References

- [NIST FIPS 203](https://csrc.nist.gov/pubs/fips/203/final) — ML-KEM Standard
- [lattice-estimator](https://github.com/malb/lattice-estimator) — Albrecht et al.
- [From Perfect to Approximate Hints](https://eprint.iacr.org/2025/1621) — S&P 2026

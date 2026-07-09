# ML-KEM Security Analysis using lattice-estimator

Security bit estimation for ML-KEM parameter sets using the lattice-estimator.

See `RESEARCH_NOTE.md` for the full research note, including the
Appendix (2026-07) on the MATZOV ζ/t grid-resolution mechanism referenced
below.

---

## Environment

- SageMath 9.x
- lattice-estimator (Albrecht et al.) — not on PyPI, install from source (see below)
- Docker

## 재현 방법

`lattice-estimator`는 PyPI 패키지가 아니라 소스에서 직접 받아야 합니다
(`pip install lattice-estimator`는 존재하지 않는 패키지이므로 동작하지 않습니다).

```bash
# SageMath Docker 컨테이너 실행
docker run --rm -it -v ${PWD}:/home/sage/work sagemath/sagemath bash

# lattice-estimator 소스 클론
git clone https://github.com/malb/lattice-estimator.git
cd lattice-estimator

# 이 레포의 스크립트를 lattice-estimator 루트에 복사한 뒤 실행
sage estimate_mlkem.sage
```

또는 이미 SageMath와 lattice-estimator 소스가 준비된 환경에서:

```bash
cd /path/to/lattice-estimator
sage estimate_mlkem.sage
```

### 검증 스크립트 (results/1_core/, 2026-07 추가)

ζ 탐색 해상도 이슈를 재현하는 7개 스크립트. 각각의 역할은
`RESEARCH_NOTE.md` Appendix / Section 8 참고.

```bash
sage verify_zeta_isolated.sage
sage verify_zeta_full_scan.sage
sage verify_zeta_direct_scan.sage
sage verify_opt_step_direct.sage
sage verify_zeta_512_1024.sage
sage verify_mechanism_512.sage
sage verify_mechanism_1024.sage
```

---

## Parameters

```python
# ML-KEM-512
LWE.Parameters(n=512, q=3329, Xs=ND.CenteredBinomial(3), Xe=ND.CenteredBinomial(3))

# ML-KEM-768
LWE.Parameters(n=768, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2))

# ML-KEM-1024
LWE.Parameters(n=1024, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2))
```

---

## Results

### Primal Attack (usvp)

| Parameter | n | Security (bits) | β (BKZ block size) |
|---|---|---|---|
| ML-KEM-512 | 512 | 143.8 | 406 |
| ML-KEM-768 | 768 | 204.9 | 624 |
| ML-KEM-1024 | 1024 | 275.1 | 874 |

### Dual Hybrid Attack (LWE.dual_hybrid / MATZOV, default output)

> These ζ values are the *default* `LWE.dual_hybrid` (MATZOV) output.
> `RESEARCH_NOTE.md` Appendix (2026-07) shows these are subject to a
> confirmed grid-resolution artifact (see note under "ML-KEM-768 anomaly"
> below) — the true optimum ζ differs from what's shown here by
> 0.23-1.08 bits depending on the parameter set.

| Parameter | Security (bits) | β | p | ζ |
|---|---|---|---|---|
| ML-KEM-512 | 139.7 | 387 | 5 | 0 |
| ML-KEM-768 | 196.4 | 589 | 4 | 20 |
| ML-KEM-1024 | 262.3 | 823 | 4 | 0 |

---

## Key Observations

### 1. Dual hybrid이 primal보다 일관되게 효율적

| 파라미터셋 | Δ (bits) |
|---|---|
| ML-KEM-512 | 4.1 bits |
| ML-KEM-768 | 8.5 bits |
| ML-KEM-1024 | 12.8 bits |

격자 차원 n이 커질수록 dual hybrid의 meet-in-the-middle 단계가
더 효과적으로 작동함. sparse secret 특성이 강해질수록 이 격차는 더 벌어짐.

### 2. 보안 레벨 스케일링

| 구간 | Primal 증가 | Dual 증가 |
|---|---|---|
| 512 → 768 | +61.1 bits | +56.7 bits |
| 768 → 1024 | +70.2 bits | +65.9 bits |

정확히 64-bit씩 증가하지 않음. n이 커질수록 증가폭이 커지는 비선형 구조.

### 3. ML-KEM-768 anomaly (ζ=20) — 2026-07 업데이트

ζ: 추측(guessing) 외에 추가적인 차원 축소(dimension reduction) 크기.

| 파라미터셋 | η | MATZOV 기본 출력 | 실제 전수조사 최적값 |
|---|---|---|---|
| ML-KEM-512 | 3 | ζ=0 | ζ=16 (0.36 bits 차이) |
| ML-KEM-768 | 2 | ζ=20 | ζ=24 (0.23 bits 차이) |
| ML-KEM-1024 | 2 | ζ=0 | ζ=34 (1.08 bits 차이) |

> 원래 이 섹션은 세 값(0, 20, 0)을 각 파라미터셋의 "진짜 최적 전략"으로
> 설명했으나, `RESEARCH_NOTE.md` Appendix에서 전수조사로 확인한 결과
> 셋 다 MATZOV의 `early_abort_range(step=10)` 그리디 탐색이 t(FFT 차원)의
> 10-단위 이산화로 생기는 톱니형 비용 곡선에서 국소 극값에 조기 멈춘
> 결과임이 확인되었다. 특히 ML-KEM-1024는 1.08비트, 무시하기 어려운
> 수준의 격차. 상세 메커니즘은 `RESEARCH_NOTE.md` Appendix Experiment
> 5-6 참고.

→ dual hybrid optimizer가 η(secret 분포 폭)와 n(격자 차원)의
상호작용에 민감하게 반응하는 것은 사실이지만, 그 출력값 자체를
그대로 신뢰하기 전에 그리드 해상도 검증이 필요함을 보여주는 사례.

---

## Sparse Secret Analysis

ML-KEM-768, dual hybrid 기준 — Hamming weight별 보안 비트.

| Hamming weight (h) | Security (bits) | Δ from standard | ζ |
|---|---|---|---|
| Standard (CBD η=2) | 196.4 | — | 20 |
| h=200 | 190.9 | -5.5 | 20 |
| h=100 | 181.1 | -15.3 | 30 |
| h=50 | 169.9 | -26.5 | 0 |
| h=30 | 164.1 | -32.3 | 0 |

h가 낮을수록 보안 비트가 급격히 감소.
h=30에서 표준 대비 -32.3 bits — NIST 128-bit 목표(Level 1)에
근접하는 수준의 약화.

이 결과는 아래 논문의 핵심 발견과 직접 연결됨:

> "From Perfect to Approximate Hints" (S&P 2026, Changmin Lee et al.)
> — sparse secret일수록 hint 하나당 LWE 공격 효과가 증폭됨.

---

## References

- [NIST FIPS 203](https://csrc.nist.gov/pubs/fips/203/final) — ML-KEM Standard
- [lattice-estimator](https://github.com/malb/lattice-estimator) — Albrecht et al.
- [ePrint 2026/1081] Hhan et al., "From Perfect to Approximate Hints", S&P 2026
- 관련 구현: [ml-kem-ntt](https://github.com/hs-pqc/ml-kem-ntt)

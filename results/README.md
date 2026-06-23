# Results

실험 환경: SageMath + lattice-estimator 0.1.0

| 파일 | 내용 |
|---|---|
| mlkem_security.txt | ML-KEM 세 파라미터셋 primal/dual 보안 추정 |
| mldsa_security.txt | ML-DSA 세 파라미터셋 primal/dual 보안 추정 |
| ntruplus_security.txt | NTRU+ 네 파라미터셋 primal/dual 보안 추정 |
| sparse_secret.txt | ML-KEM-768 sparse secret 분석 (h=30~200) |

## ζ Anomaly 비교 분석

세 표준의 dual hybrid attack ζ 값 비교:

| 표준 | 파라미터셋 | ζ |
|---|---|---|
| ML-KEM | 512 | 0 |
| ML-KEM | 768 | **20** |
| ML-KEM | 1024 | 0 |
| ML-DSA | 44 | 10 |
| ML-DSA | 65 | 20 |
| ML-DSA | 87 | 30 |
| NTRU+ | 576 | 20 |
| NTRU+ | 768 | 20 |
| NTRU+ | 864 | 20 |
| NTRU+ | 1152 | 30 |

## 주요 관찰

**1. ML-KEM만 ζ=0인 파라미터셋이 존재**
ML-DSA와 NTRU+는 모든 파라미터셋에서 ζ>0.
ML-KEM은 512와 1024에서 ζ=0, 768에서만 ζ=20.

**2. ML-DSA는 n이 커질수록 ζ도 단조 증가**
44(ζ=10) → 65(ζ=20) → 87(ζ=30)
q=8380417로 고정된 상태에서 n 증가에 따른 패턴.

**3. NTRU+는 q/n 비율이 일정 (≈6.0)**
q가 n에 비례해서 증가하는 구조.
576~864에서 ζ=20 유지, 1152에서 ζ=30으로 증가.

**4. 가설: q/n 비율이 ζ anomaly 발생 조건을 결정**
- ML-KEM: q/n ≈ 13 (512), 4.3 (768), 3.3 (1024)
- NTRU+: q/n ≈ 6.0 (고정)
- ML-DSA: q/n >> 100 (매우 큼)

→ q/n 비율과 ζ 발생 조건의 관계를 추가 실험으로 검증 필요.

## 향후 실험

- [ ] q를 변수로 바꿔가면서 ζ 변화 관찰
- [ ] n을 600~900 구간으로 세밀하게 분석
- [ ] ζ anomaly 발생 조건 수식화

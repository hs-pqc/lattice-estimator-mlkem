# Results

실험 환경: SageMath + lattice-estimator 0.1.0

| 파일 | 내용 |
|---|---|
| mlkem_security.txt | ML-KEM 세 파라미터셋 primal/dual 보안 추정 |
| mldsa_security.txt | ML-DSA 세 파라미터셋 primal/dual 보안 추정 |
| ntruplus_security.txt | NTRU+ 네 파라미터셋 primal/dual 보안 추정 |
| sparse_secret.txt | ML-KEM-768 sparse secret 분석 (h=30~200) |
| zeta_search.txt | n=600~900 구간 ζ 변화 관찰 (q=3329, η=2) |
| zeta_eta.txt | η 변화에 따른 ζ 거동 (n=768, q=3329) |
| zeta_n_eta.txt | n=512에서 η 변화에 따른 ζ 거동 |

## ζ 거동 분석

### 세 표준 비교

| 표준 | 파라미터셋 | ζ |
|---|---|---|
| ML-KEM | 512 (η=3) | 0 |
| ML-KEM | 768 (η=2) | 20 |
| ML-KEM | 1024 (η=2) | 0 |
| ML-DSA | 44 | 10 |
| ML-DSA | 65 | 20 |
| ML-DSA | 87 | 30 |
| NTRU+ | 576 | 20 |
| NTRU+ | 768 | 20 |
| NTRU+ | 864 | 20 |
| NTRU+ | 1152 | 30 |

### n과 η의 상호작용

| | η=1 | η=2 | η=3 | η=4 |
|---|---|---|---|---|
| n=512 | ζ=20 | ζ=10 | ζ=0 | ζ=10 |
| n=768 | ζ=0 | ζ=20 | ζ=20 | ζ=20 |

n=512와 n=768에서 η에 따른 ζ 거동이 반대 방향.
ζ는 n과 η의 상호작용으로 결정되며 어느 하나만으로 설명 불가.

### 핵심 발견

ML-KEM-512는 η=3을 채택해서 결과적으로 ζ=0.
η=2였다면 ζ=10이 나왔을 것.
→ ML-KEM 파라미터 설계에서 η=3 선택이 dual hybrid
차원 축소 효과를 차단한 결과를 낳았음.
의도적 설계인지 여부는 추가 분석 필요.

## 향후 실험

- [ ] q를 변수로 바꿔가면서 ζ 변화 관찰
- [ ] n=512 구간 세밀하게 분석
- [ ] ζ 거동 조건 수식화

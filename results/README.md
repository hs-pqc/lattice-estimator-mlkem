# Results

실험 환경: SageMath + lattice-estimator 0.1.0

| 파일 | 내용 |
|---|---|
| mlkem_security.txt | ML-KEM 세 파라미터셋 primal/dual 보안 추정 |
| sparse_secret.txt | ML-KEM-768 sparse secret 분석 (h=30~200) |

## 주요 발견

**1. Dual hybrid가 primal보다 일관되게 효율적**
- ML-KEM-512: Δ=4.1 bits
- ML-KEM-768: Δ=8.5 bits
- ML-KEM-1024: Δ=12.8 bits

**2. ML-KEM-768 ζ=20 anomaly**
- ML-KEM-512: ζ=0
- ML-KEM-768: ζ=20 ← optimizer가 차원 축소를 선택
- ML-KEM-1024: ζ=0
- η=2와 n=768의 균형점에서 발생

**3. Sparse secret 보안 감소**
- h=30: 표준 대비 -32.3 bits
- h=50: 표준 대비 -26.5 bits
- h=100: 표준 대비 -15.3 bits
- h=200: 표준 대비 -5.5 bits

→ S&P 2026 "From Perfect to Approximate Hints"의 핵심 발견과 직접 연결.

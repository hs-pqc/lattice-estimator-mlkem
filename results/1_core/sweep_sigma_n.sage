"""
sweep_sigma_n.sage

Frodo640만 MATZOV vs Guo-Johansson gap이 뒤집히는 이유를 진단하기 위한 통제 실험.

Frodo640 실제 파라미터: n=640, q=2^15=32768, Xs=Xe=DiscreteGaussian(σ≈2.8)
ML-KEM-512는 n=512로 비슷하지만 q=3329(훨씬 작음), 분포도 CenteredBinomial로 다름.
그래서 1차 진단(fix_and_diagnose.sage)에서 본 "n이 클수록 MATZOV가 유리해진다"는
그룹 내부 패턴이 n 자체 때문인지, σ(에러 폭) 때문인지, q 때문인지 뒤섞여 있었음.

이 스크립트는 n과 σ를 각각 하나씩 고정하고 나머지 하나만 바꿔서, 어느 쪽이
Frodo640의 역전(gap>0)을 실제로 설명하는지 분리한다.

실험 A: n=640, q=32768 고정, σ ∈ {1.0, 1.6, 2.2, 2.8, 3.4, 4.0} 스윕
실험 B: σ=2.8, q=32768 고정, n ∈ {400, 550, 640, 800, 1000, 1200} 스윕
  (n=640, σ=2.8 지점은 두 실험 모두에 포함 — Frodo640 자체이므로 sanity check 겸용)

주의: MATZOV 쪽(two_stage_search, verify_parallel_v2.sage에서 load)은 스킴마다
1.5~3분 걸림. 총 11개 지점(A 6개 + B 6개 - 중복 1개) x 이 정도면 20~35분 예상.
GJ 쪽(dual_hybrid fft=True)은 지점당 1초 내외로 빠름.
"""

import time
from estimator import *

load("verify_parallel_v2.sage")


def gj_log2_rop(params, opt_step=4):
    r = dual_hybrid(params, opt_step=opt_step, fft=True)
    return float(log(r["rop"], 2).n())


def run_point(label, params):
    t0 = time.time()
    matzov_result = two_stage_search(params)
    matzov_v = matzov_result["log2_rop"]
    gj_v = gj_log2_rop(params)
    gap = matzov_v - gj_v
    elapsed = time.time() - t0
    flag = "  <-- GJ 우위 (역전)" if gap > 0 else ""
    print(f"{label:<28} matzov={matzov_v:>8.3f}  gj={gj_v:>8.3f}  gap={gap:>+8.3f}{flag}  ({elapsed:.0f}s)")
    return gap


print("=" * 78)
print("실험 A: n=640, q=32768 고정, sigma 스윕")
print("=" * 78)
Q_FIXED = 32768  # 2^15, Frodo640과 동일
N_FIXED_A = 640

for sigma in [1.0, 1.6, 2.2, 2.8, 3.4, 4.0]:
    params = LWE.Parameters(
        n=N_FIXED_A, q=Q_FIXED,
        Xs=ND.DiscreteGaussian(sigma), Xe=ND.DiscreteGaussian(sigma),
        tag=f"sweep-n640-sigma{sigma}",
    ).normalize()
    run_point(f"n=640, sigma={sigma}", params)

print()
print("=" * 78)
print("실험 B: sigma=2.8, q=32768 고정, n 스윕")
print("=" * 78)
SIGMA_FIXED = 2.8  # Frodo640과 동일

for n in [400, 550, 640, 800, 1000, 1200]:
    params = LWE.Parameters(
        n=n, q=Q_FIXED,
        Xs=ND.DiscreteGaussian(SIGMA_FIXED), Xe=ND.DiscreteGaussian(SIGMA_FIXED),
        tag=f"sweep-sigma2.8-n{n}",
    ).normalize()
    run_point(f"n={n}, sigma=2.8", params)

print()
print("(참고: 위 n=640,sigma=2.8 지점이 두 실험에 각각 한 번씩 나오는데, 값이 서로")
print(" 비슷하게 나오면 계산 자체는 안정적이라는 뜻 — Frodo640 원본 gap=+0.830과도 비교해볼 것)")

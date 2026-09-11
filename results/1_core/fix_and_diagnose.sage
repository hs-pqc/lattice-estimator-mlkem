"""
fix_and_diagnose.sage

compare_dual_attacks.sage 1차 실행 결과에서 나온 두 가지 후속 작업:

1) zeta=0 아티팩트 제거
   ML-KEM-512, ML-KEM-1024, ML-DSA-87에서 matzov_indirect가 zeta=0을 "최적"으로
   보고했는데, 이건 실제 최적점이 아니라 LWE.dual_hybrid(=matzov.cost 내부)가 쓰는
   early_abort_range(step=10) 그리드 서치가 첫 상승점에서 조기 종료하는 아티팩트
   (Issue #219 Appendix에서 이미 규명됨). verify_parallel_v2.sage의 two_stage_search
   (coarse-to-fine, greedy가 아니라 전 구간 스캔 후 window 확장)로 다시 구해서 교체한다.

2) MATZOV vs Guo-Johansson 역전 패턴 진단
   compare_dual_attacks.sage 1차 결과에서 Frodo640만 GJ가 MATZOV보다 더 싼 공격을
   찾아냈다 (169.316 vs 170.146, gap -0.83bit). 나머지 스킴은 전부 MATZOV가 우위.
   이 스크립트는 (n, q, 표준편차) 등 파라미터 특성과 gap의 관계를 표로 뽑아서
   Frodo640이 왜 예외인지 후보 가설을 좁힌다. 인과 증명이 아니라 패턴 탐색 단계.
"""

import os
import time
from estimator import *

# verify_parallel_v2.sage에 이미 구현된 병렬(Pool 기반) two_stage_search를 그대로 불러와 씀.
# load()는 같은 디렉토리(results/1_core)에서 실행한다는 전제. load된 파일의
# `if __name__ == "__main__": main()` 은 load 시점엔 실행되지 않으므로 안전함.
load("verify_parallel_v2.sage")


# ---------------------------------------------------------------------------
# 1. zeta=0 아티팩트 제거 — 1차 실행에서 zeta=0이 나온 4개 스킴만 다시 계산
# ---------------------------------------------------------------------------
AFFECTED_SCHEMES = {
    "ML-KEM-512": LWE.Parameters(n=512, q=3329, Xs=ND.CenteredBinomial(3), Xe=ND.CenteredBinomial(3), tag="ML-KEM-512").normalize(),
    "ML-KEM-1024": LWE.Parameters(n=1024, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="ML-KEM-1024").normalize(),
    "ML-DSA-87": LWE.Parameters(n=7 * 256, q=8380417, Xs=ND.UniformMod(2), Xe=ND.UniformMod(2), m=8 * 256, tag="ML-DSA-87").normalize(),
}

print("=" * 70)
print("1. zeta=0 아티팩트 제거 (coarse-to-fine 재계산)")
print("=" * 70)

corrected_matzov = {}
for name, params in AFFECTED_SCHEMES.items():
    result = two_stage_search(params)
    corrected_matzov[name] = result
    print(f"{name}: zeta={result['zeta']}, t={result['t']}, log2(rop)={result['log2_rop']:.3f}  ({result['elapsed_s']:.1f}s)")

print()
print("(참고용 GROUND_TRUTH: ML-KEM-512 zeta=14 log2_rop=139.057, "
      "ML-KEM-1024 zeta=32 log2_rop=261.143 — 이 값과 비교해서 재계산이 맞는지 확인)")


# ---------------------------------------------------------------------------
# 2. MATZOV vs Guo-Johansson gap과 파라미터 특성 진단
# ---------------------------------------------------------------------------
ALL_SCHEMES_FOR_DIAGNOSIS = {
    "ML-KEM-512": AFFECTED_SCHEMES["ML-KEM-512"],
    "ML-KEM-768": LWE.Parameters(n=768, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="ML-KEM-768").normalize(),
    "ML-KEM-1024": AFFECTED_SCHEMES["ML-KEM-1024"],
    "ML-DSA-44": LWE.Parameters(n=4 * 256, q=8380417, Xs=ND.UniformMod(2), Xe=ND.UniformMod(2), m=4 * 256, tag="ML-DSA-44").normalize(),
    "ML-DSA-65": LWE.Parameters(n=5 * 256, q=8380417, Xs=ND.UniformMod(4), Xe=ND.UniformMod(4), m=6 * 256, tag="ML-DSA-65").normalize(),
    "ML-DSA-87": AFFECTED_SCHEMES["ML-DSA-87"],
    "NTRU+": LWE.Parameters(n=576, q=3457, Xs=ND.Binary, Xe=ND.SparseTernary(576, 192), m=576, tag="NTRU+576").normalize(),
    "Frodo640": schemes.Frodo640,
    "Frodo976": schemes.Frodo976,
    "Frodo1344": schemes.Frodo1344,
}

GJ_LOG2_ROP_FROM_FIRST_RUN = {
    "ML-KEM-512": 143.788, "ML-KEM-768": 203.788, "ML-KEM-1024": 273.817,
    "ML-DSA-44": 128.286, "ML-DSA-65": 183.936, "ML-DSA-87": 232.118,
    "NTRU+": 132.363,
    "Frodo640": 169.316, "Frodo976": 232.696, "Frodo1344": 296.296,
}

MATZOV_LOG2_ROP_FROM_FIRST_RUN = {
    "ML-KEM-768": 196.366, "ML-DSA-44": 126.682, "ML-DSA-65": 180.456,
    "NTRU+": 130.716, "Frodo640": 170.146, "Frodo976": 231.495, "Frodo1344": 281.739,
}
for name, result in corrected_matzov.items():
    MATZOV_LOG2_ROP_FROM_FIRST_RUN[name] = result["log2_rop"]

print()
print("=" * 70)
print("2. gap = MATZOV - GJ (양수: MATZOV가 더 싼 공격 / 음수: GJ가 더 싼 공격)")
print("=" * 70)
print(f"{'scheme':<14} {'n':>5} {'log2(q)':>8} {'gap(bit)':>10}")
for name, params in ALL_SCHEMES_FOR_DIAGNOSIS.items():
    matzov_v = MATZOV_LOG2_ROP_FROM_FIRST_RUN[name]
    gj_v = GJ_LOG2_ROP_FROM_FIRST_RUN[name]
    gap = matzov_v - gj_v
    log2q = float(log(params.q, 2).n())
    flag = "  <-- 역전" if gap < 0 else ""
    print(f"{name:<14} {params.n:>5} {log2q:>8.2f} {gap:>+10.3f}{flag}")

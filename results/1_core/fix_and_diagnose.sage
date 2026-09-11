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
from estimator.lwe_dual import dual_hybrid

# verify_parallel_v2.sage와 같은 폴더에서 실행한다는 전제 (같은 디렉토리에 있어야
# import 없이 그냥 재정의해서 씀 — .sage 파일 간 import가 번거로워서 핵심 함수만 복붙)
from estimator.lwe_dual import matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default


def cost_fixed_zeta_two_level_t(k_enum_val, params, t_coarse_step=10, t_window=15, t_max_doublings=3, red_cost_model=red_cost_model_default):
    for p in early_abort_range(2, params.q):
        t_max = params.n - k_enum_val
        coarse_results = []
        for t in range(0, t_max + 1, t_coarse_step):
            precision = 1
            max_beta = max(min(params.m - k_enum_val - t, max_beta_global), 40 + precision)
            with local_minimum(40, max_beta, precision=precision) as it:
                for beta in it:
                    cost = matzov.cost(beta, params, p=p[0], k_enum=k_enum_val, k_fft=t, red_cost_model=red_cost_model)
                    it.update(cost)
            coarse_results.append((t, it.y))
        t_c, best_coarse = min(coarse_results, key=lambda x: x[1]["rop"])

        window = t_window
        doublings = 0
        p_best = best_coarse
        while True:
            lo = max(0, t_c - window)
            hi = min(t_max, t_c + window)
            fine_results = []
            for t in range(lo, hi + 1):
                precision = 1
                max_beta = max(min(params.m - k_enum_val - t, max_beta_global), 40 + precision)
                with local_minimum(40, max_beta, precision=precision) as it:
                    for beta in it:
                        cost = matzov.cost(beta, params, p=p[0], k_enum=k_enum_val, k_fft=t, red_cost_model=red_cost_model)
                        it.update(cost)
                fine_results.append((t, it.y))
            t_opt, p_best = min(fine_results, key=lambda x: x[1]["rop"])
            hit_boundary = (t_opt == lo or t_opt == hi)
            if not hit_boundary or doublings >= t_max_doublings:
                break
            window *= 2
            doublings += 1

        p[1].update(p_best)
        if p[1].y["t"] == 0 and p[0] > 2:
            break
    return p[1].y


def coarse_scan(params, step=10):
    n = params.n
    results = []
    for zeta in range(0, n, step):
        r = cost_fixed_zeta_two_level_t(zeta, params)
        log2rop = float(log(r["rop"], 2).n())
        results.append((zeta, log2rop, r.get("t"), r.get("beta")))
    zeta_c, rop_c, t_c, beta_c = min(results, key=lambda x: x[1])
    return zeta_c, rop_c, results


def fine_scan(params, center, window):
    n = params.n
    lo = max(0, center - window)
    hi = min(n - 1, center + window)
    results = []
    for zeta in range(lo, hi + 1):
        r = cost_fixed_zeta_two_level_t(zeta, params)
        log2rop = float(log(r["rop"], 2).n())
        results.append((zeta, log2rop, int(r.get("t")), int(r.get("beta"))))
    zeta_opt, rop_opt, t_opt, beta_opt = min(results, key=lambda x: x[1])
    hit_boundary = (zeta_opt == lo or zeta_opt == hi)
    return zeta_opt, rop_opt, t_opt, beta_opt, hit_boundary


def two_stage_search(params, coarse_step=10, initial_window=15, max_doublings=3):
    """verify_parallel_v2.sage의 병렬(Pool) 버전 대신 순차 버전. 4개 스킴만 다시
    돌리는 거라 순차로도 충분히 빠름 — 병렬 워커 초기화 오버헤드를 피하려고 일부러 순차로 씀."""
    t0 = time.time()
    zeta_c, rop_c, coarse_results = coarse_scan(params, step=coarse_step)
    window = initial_window
    doublings = 0
    while True:
        zeta_opt, rop_opt, t_opt, beta_opt, hit_boundary = fine_scan(params, zeta_c, window)
        if not hit_boundary or doublings >= max_doublings:
            break
        window *= 2
        doublings += 1
    elapsed = time.time() - t0
    return dict(zeta=zeta_opt, t=t_opt, beta=beta_opt, log2_rop=rop_opt, elapsed_s=elapsed)


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
# 1차 실행 결과(zeta=0 아티팩트 있는 3곳은 위에서 재계산한 값으로 교체)를 모아서
# gap = log2_rop(matzov) - log2_rop(guo_johansson) 을 계산. 음수면 GJ가 더 싼 공격
# (=MATZOV보다 GJ가 우위) — Frodo640만 이 부호였음.

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

# 1차 실행(compare_dual_attacks.sage)에서 나온 guo_johansson_fft_true 값을 그대로 하드코딩.
# (재계산 안 해도 됨 — GJ 경로는 zeta=0 버그가 없었으므로 1차 값 그대로 신뢰 가능)
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
# zeta=0이었던 3곳은 위에서 재계산한 corrected_matzov로 덮어씀
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

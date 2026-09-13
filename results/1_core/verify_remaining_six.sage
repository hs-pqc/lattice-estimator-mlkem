"""
verify_remaining_six.sage

RESEARCH_NOTE.md Experiment 9 테이블에서 †(raw default 값, 미검증)로 표시했던
6개 스킴을 coarse-to-fine(two_stage_search)으로 재검증한다:
    ML-DSA-44, ML-DSA-65, NTRU+, Frodo640, Frodo976, Frodo1344

목적: default LWE.dual_hybrid()가 이 6곳에서도 ζ=0류 아티팩트나 그보다 작은
sawtooth 오차로 실제보다 나쁜(=너무 큰) MATZOV log2(rop)를 보고하고 있는지 확인.
특히 Frodo976은 문서에 이미 "t=0으로 잘못 착지하는 케이스"로 기록돼 있어
gap 값이 바뀔 가능성이 가장 높음. Frodo640은 이번 연구의 핵심 사례라 재검증이
특히 중요함 — gap=+0.830(GJ 우위)이라는 결론이 재검증 후에도 유지되는지가 관건.

GJ(dual_hybrid fft=True) 쪽은 처음부터 zeta=0 버그가 없었으므로 재검증 없이
Experiment 9의 기존 값을 그대로 사용.
"""

import time
from estimator import *
from estimator.lwe_dual import dual_hybrid, matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default
import os
from multiprocessing import Pool


def cost_fixed_zeta_fast(k_enum_val, params, red_cost_model=red_cost_model_default):
    for p in early_abort_range(2, params.q):
        for k_fft in early_abort_range(0, params.n - k_enum_val, 10):
            precision = 1
            max_beta = max(min(params.m - k_enum_val - k_fft[0], max_beta_global), 40 + precision)
            with local_minimum(40, max_beta, precision=precision) as it:
                for beta in it:
                    cost = matzov.cost(beta, params, p=p[0], k_enum=k_enum_val, k_fft=k_fft[0], red_cost_model=red_cost_model)
                    it.update(cost)
                k_fft[1].update(it.y)
        p[1].update(k_fft[1].y)
        if p[1].y["t"] == 0 and p[0] > 2:
            break
    return p[1].y


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
        r = cost_fixed_zeta_fast(zeta, params)
        log2rop = float(log(r["rop"], 2).n())
        results.append((zeta, log2rop, r.get("t"), r.get("beta")))
    zeta_c, rop_c, t_c, beta_c = min(results, key=lambda x: x[1])
    return zeta_c, rop_c, results


_WORKER_PARAMS = None


def _init_worker(params):
    global _WORKER_PARAMS
    _WORKER_PARAMS = params


def _fine_worker(zeta):
    r = cost_fixed_zeta_two_level_t(zeta, _WORKER_PARAMS)
    log2rop = float(log(r["rop"], 2).n())
    return (zeta, log2rop, int(r.get("t")), int(r.get("beta")))


def fine_scan_parallel(params, center, window, nproc=None):
    if nproc is None:
        nproc = os.cpu_count() or 4
    n = params.n
    lo = max(0, center - window)
    hi = min(n - 1, center + window)
    zetas = list(range(lo, hi + 1))
    with Pool(processes=min(nproc, len(zetas)), initializer=_init_worker, initargs=(params,)) as pool:
        results = pool.map(_fine_worker, zetas)
    zeta_opt, rop_opt, t_opt, beta_opt = min(results, key=lambda x: x[1])
    hit_boundary = (zeta_opt == lo or zeta_opt == hi)
    return zeta_opt, rop_opt, t_opt, beta_opt, hit_boundary


def two_stage_search(params, coarse_step=10, initial_window=15, max_doublings=3, nproc=None):
    zeta_c, rop_c, coarse_results = coarse_scan(params, step=coarse_step)
    window = initial_window
    doublings = 0
    while True:
        zeta_opt, rop_opt, t_opt, beta_opt, hit_boundary = fine_scan_parallel(params, zeta_c, window, nproc=nproc)
        if not hit_boundary or doublings >= max_doublings:
            break
        window *= 2
        doublings += 1
    return dict(zeta=zeta_opt, t=t_opt, beta=beta_opt, log2_rop=rop_opt)


# GJ 값은 이미 검증된(zeta=0 버그 없는) Experiment 9 값을 그대로 사용
GJ_LOG2_ROP = {
    "ML-DSA-44": 128.286, "ML-DSA-65": 183.936,
    "NTRU+": 132.363,
    "Frodo640": 169.316, "Frodo976": 232.696, "Frodo1344": 296.296,
}

# 이전에 raw default로만 나왔던(†) MATZOV 값 — 비교용으로 같이 출력
MATZOV_DEFAULT_RAW = {
    "ML-DSA-44": 126.682, "ML-DSA-65": 180.456,
    "NTRU+": 130.716,
    "Frodo640": 170.146, "Frodo976": 231.495, "Frodo1344": 281.739,
}

TO_VERIFY = {
    "ML-DSA-44": LWE.Parameters(n=4 * 256, q=8380417, Xs=ND.UniformMod(2), Xe=ND.UniformMod(2), m=4 * 256, tag="ML-DSA-44").normalize(),
    "ML-DSA-65": LWE.Parameters(n=5 * 256, q=8380417, Xs=ND.UniformMod(4), Xe=ND.UniformMod(4), m=6 * 256, tag="ML-DSA-65").normalize(),
    "NTRU+": LWE.Parameters(n=576, q=3457, Xs=ND.Binary, Xe=ND.SparseTernary(576, 192), m=576, tag="NTRU+576").normalize(),
    "Frodo640": schemes.Frodo640,
    "Frodo976": schemes.Frodo976,
    "Frodo1344": schemes.Frodo1344,
}

print("=" * 90)
print("6개 미검증 스킴 coarse-to-fine 재검증")
print("=" * 90)
print(f"{'scheme':<12} {'default(raw)':>13} {'corrected':>11} {'diff':>8}  {'GJ':>10} {'new gap':>10}  winner")

for name, params in TO_VERIFY.items():
    t0 = time.time()
    result = two_stage_search(params)
    elapsed = time.time() - t0
    corrected = result["log2_rop"]
    raw = MATZOV_DEFAULT_RAW[name]
    diff = raw - corrected
    gj = GJ_LOG2_ROP[name]
    new_gap = corrected - gj
    winner = "GJ" if new_gap > 0 else "MATZOV"
    changed = "  <-- raw와 다름" if abs(diff) > 0.01 else ""
    print(f"{name:<12} {raw:>13.3f} {corrected:>11.3f} {diff:>+8.3f}  {gj:>10.3f} {new_gap:>+10.3f}  {winner}{changed}  ({elapsed:.0f}s)")

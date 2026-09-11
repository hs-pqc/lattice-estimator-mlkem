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

주의: MATZOV 쪽(two_stage_search)은 스킴마다 1.5~3분 걸림.
총 11개 지점(A 6개 + B 6개 - 중복 1개) x 이 정도면 20~35분 예상.
GJ 쪽(dual_hybrid fft=True)은 지점당 1초 내외로 빠름.
"""

import time
from estimator import *
from estimator.lwe_dual import dual_hybrid, matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default
import os
from multiprocessing import Pool

# NOTE: verify_parallel_v2.sage를 load()로 불러오면 Sage에서는 __name__이 "__main__"으로
# 설정되어 그 파일의 main()(ML-KEM 3종 재검증, 실험당 5~6분)까지 매번 같이 실행돼버림.
# 그래서 함수 본체만 직접 복사해서 씀 (verify_parallel_v2.sage와 100% 동일한 로직).


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

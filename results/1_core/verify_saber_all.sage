import os
import time
from multiprocessing import Pool
from estimator import *
from estimator.lwe_dual import matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default


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
    if nproc is None:
        nproc = os.cpu_count() or 4
    t0 = time.time()
    zeta_c, rop_c, coarse_results = coarse_scan(params, step=coarse_step)
    window = initial_window
    doublings = 0
    while True:
        zeta_opt, rop_opt, t_opt, beta_opt, hit_boundary = fine_scan_parallel(params, zeta_c, window, nproc=nproc)
        if not hit_boundary or doublings >= max_doublings:
            break
        window *= 2
        doublings += 1
    elapsed = time.time() - t0
    return dict(zeta=zeta_opt, t=t_opt, beta=beta_opt, log2_rop=rop_opt,
                zeta_coarse=zeta_c, final_window=window, doublings=doublings, elapsed_s=elapsed)


def main():
    print(f"{'Parameter':<12}{'def zeta':>10}{'def t':>8}{'def log2rop':>14}{'true zeta':>11}{'true t':>8}{'true log2rop':>14}{'gap':>8}{'time':>8}")
    print("-" * 95)
    for name in ["LightSaber", "Saber", "FireSaber"]:
        params = getattr(schemes, name)
        default = LWE.dual_hybrid(params)
        default_log2rop = float(log(default["rop"], 2).n())

        result = two_stage_search(params)

        gap = default_log2rop - result["log2_rop"]
        print(f"{name:<12}{default['zeta']:>10}{default.get('t', '-'):>8}{default_log2rop:>14.3f}"
              f"{result['zeta']:>11}{result['t']:>8}{result['log2_rop']:>14.3f}{gap:>8.3f}{result['elapsed_s']:>8.1f}s")

if __name__ == "__main__":
    main()

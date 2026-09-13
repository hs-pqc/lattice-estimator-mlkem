"""
sweep_2d_grid.sage
"""
import json
import os
import time
from estimator import *
from estimator.lwe_dual import dual_hybrid, matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default
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


def gj_log2_rop(params, opt_step=4):
    r = dual_hybrid(params, opt_step=opt_step, fft=True)
    return float(log(r["rop"], 2).n())


OUT_PATH = "/home/sage/lattice-estimator-mlkem/results/1_core/results/sweep_2d_grid.json"
Q_FIXED = 32768

n_values = [700, 730, 760, 800]
sigma_values = [2.5, 2.65, 2.8, 2.95, 3.1]

results = []
if os.path.exists(OUT_PATH):
    with open(OUT_PATH) as f:
        results = json.load(f)
done_keys = {(r["n"], r["sigma"]) for r in results}

print(f"Grid: n={n_values} x sigma={sigma_values} = {len(n_values)*len(sigma_values)} points")
print(f"Already done: {len(done_keys)}")

for n in n_values:
    for sigma in sigma_values:
        if (n, sigma) in done_keys:
            print(f"skip n={n}, sigma={sigma} (already done)")
            continue
        params = LWE.Parameters(
            n=n, q=Q_FIXED,
            Xs=ND.DiscreteGaussian(sigma), Xe=ND.DiscreteGaussian(sigma),
            tag=f"grid-n{n}-sigma{sigma}",
        ).normalize()
        t0 = time.time()
        matzov_v = two_stage_search(params)["log2_rop"]
        gj_v = gj_log2_rop(params)
        gap = matzov_v - gj_v
        elapsed = time.time() - t0
        flag = "GJ" if gap > 0 else "MATZOV"
        print(f"n={n:4d} sigma={sigma:.1f}  matzov={matzov_v:8.3f}  gj={gj_v:8.3f}  gap={gap:+8.3f}  winner={flag}  ({elapsed:.0f}s)")
        results.append({"n": int(n), "sigma": float(sigma), "matzov": float(matzov_v), "gj": float(gj_v), "gap": float(gap)})
        with open(OUT_PATH, "w") as f:
            json.dump(results, f, indent=2)

print(f"\nDone. Wrote {len(results)} points to {OUT_PATH}")

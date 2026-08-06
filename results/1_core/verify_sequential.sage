import time
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


def fine_scan_sequential(params, center, window):
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


def two_stage_search_sequential(params, coarse_step=10, initial_window=15, max_doublings=3):
    t0 = time.time()
    zeta_c, rop_c, coarse_results = coarse_scan(params, step=coarse_step)
    window = initial_window
    doublings = 0
    while True:
        zeta_opt, rop_opt, t_opt, beta_opt, hit_boundary = fine_scan_sequential(params, zeta_c, window)
        if not hit_boundary or doublings >= max_doublings:
            break
        window *= 2
        doublings += 1
    elapsed = time.time() - t0
    return dict(zeta=zeta_opt, t=t_opt, beta=beta_opt, log2_rop=rop_opt,
                zeta_coarse=zeta_c, final_window=window, doublings=doublings,
                elapsed_s=elapsed, n_fine_evals=(2*window+1))


PARAM_SETS = {
    "ML-KEM-512": LWE.Parameters(n=512, q=3329, Xs=ND.CenteredBinomial(3), Xe=ND.CenteredBinomial(3), tag="ML-KEM-512").normalize(),
    "ML-KEM-768": LWE.Parameters(n=768, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="ML-KEM-768").normalize(),
    "ML-KEM-1024": LWE.Parameters(n=1024, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="ML-KEM-1024").normalize(),
}

GROUND_TRUTH = {
    "ML-KEM-512": dict(zeta=14, t=34, log2_rop=139.057),
    "ML-KEM-768": dict(zeta=23, t=59, log2_rop=195.554),
    "ML-KEM-1024": dict(zeta=32, t=82, log2_rop=261.143),
}

def main():
    print(f"{'Parameter':<14}{'zeta':>6}{'t':>6}{'log2(rop)':>12}{'match':>10}{'time(s)':>10}{'fine evals':>12}")
    print("-" * 70)
    for name, params in PARAM_SETS.items():
        result = two_stage_search_sequential(params)
        gt = GROUND_TRUTH[name]
        match = "OK" if (result["zeta"] == gt["zeta"] and abs(result["log2_rop"] - gt["log2_rop"]) < 0.05) else "MISMATCH"
        print(f"{name:<14}{result['zeta']:>6}{result['t']:>6}{result['log2_rop']:>12.3f}{match:>10}{result['elapsed_s']:>10.1f}{result['n_fine_evals']:>12}")

if __name__ == "__main__":
    main()

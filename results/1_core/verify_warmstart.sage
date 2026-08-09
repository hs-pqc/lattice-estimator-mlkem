import os
import time
from estimator import *
from estimator.lwe_dual import matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default


def cost_at_t(k_enum_val, t, p_val, params, red_cost_model=red_cost_model_default):
    """Single (zeta, t, p) evaluation: BKZ beta local_minimum search only."""
    precision = 1
    max_beta = max(min(params.m - k_enum_val - t, max_beta_global), 40 + precision)
    with local_minimum(40, max_beta, precision=precision) as it:
        for beta in it:
            cost = matzov.cost(beta, params, p=p_val, k_enum=k_enum_val, k_fft=t, red_cost_model=red_cost_model)
            it.update(cost)
    return it.y


def cost_fixed_zeta_fast(k_enum_val, params, red_cost_model=red_cost_model_default):
    """Original coarse-scan cost function (unchanged, used for the outer zeta coarse pass)."""
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


def cost_fixed_zeta_warmstart(k_enum_val, params, t_hint=None, t_coarse_step=10, t_window=15,
                               t_max_doublings=3, red_cost_model=red_cost_model_default, eval_counter=None):
    """
    Same two-level (coarse t -> fine t window) search as cost_fixed_zeta_two_level_t,
    EXCEPT: if t_hint is given, skip the full 0..t_max coarse sweep entirely and open
    the fine window directly around t_hint. Falls back to the full coarse sweep when
    t_hint is None (i.e. the very first zeta processed).
    """
    for p in early_abort_range(2, params.q):
        t_max = params.n - k_enum_val

        if t_hint is None:
            coarse_results = []
            for t in range(0, t_max + 1, t_coarse_step):
                y = cost_at_t(k_enum_val, t, p[0], params, red_cost_model)
                if eval_counter is not None:
                    eval_counter[0] += 1
                coarse_results.append((t, y))
            t_c, best_coarse = min(coarse_results, key=lambda x: x[1]["rop"])
        else:
            t_c = min(max(t_hint, 0), t_max)
            best_coarse = cost_at_t(k_enum_val, t_c, p[0], params, red_cost_model)
            if eval_counter is not None:
                eval_counter[0] += 1

        window = t_window
        doublings = 0
        p_best = best_coarse
        while True:
            lo = max(0, t_c - window)
            hi = min(t_max, t_c + window)
            fine_results = []
            for t in range(lo, hi + 1):
                if t == t_c and t_hint is not None:
                    fine_results.append((t, best_coarse))
                    continue
                y = cost_at_t(k_enum_val, t, p[0], params, red_cost_model)
                if eval_counter is not None:
                    eval_counter[0] += 1
                fine_results.append((t, y))
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
    return zeta_c, rop_c, t_c, results


def fine_scan_warmstart(params, center, window, t_seed):
    """
    Sequential fine scan over zeta in [center-window, center+window], carrying the
    optimal t found at each zeta forward as the t_hint for its neighbor. Starts from
    the coarse scan's own t estimate at the center (t_seed) rather than re-deriving
    it from scratch.
    """
    n = params.n
    lo = max(0, center - window)
    hi = min(n - 1, center + window)
    eval_counter = [0]

    results = {}
    t_hint = t_seed
    for zeta in range(center, hi + 1):
        r = cost_fixed_zeta_warmstart(zeta, params, t_hint=t_hint, eval_counter=eval_counter)
        results[zeta] = r
        t_hint = int(r.get("t"))

    t_hint = t_seed
    for zeta in range(center - 1, lo - 1, -1):
        r = cost_fixed_zeta_warmstart(zeta, params, t_hint=t_hint, eval_counter=eval_counter)
        results[zeta] = r
        t_hint = int(r.get("t"))

    zeta_opt = min(results, key=lambda z: results[z]["rop"])
    r_opt = results[zeta_opt]
    log2rop = float(log(r_opt["rop"], 2).n())
    hit_boundary = (zeta_opt == lo or zeta_opt == hi)
    return zeta_opt, log2rop, int(r_opt.get("t")), int(r_opt.get("beta")), hit_boundary, eval_counter[0]


def two_stage_search_warmstart(params, coarse_step=10, initial_window=15, max_doublings=3):
    t0 = time.time()
    zeta_c, rop_c, t_c, coarse_results = coarse_scan(params, step=coarse_step)
    print(f"  coarse landed at zeta={zeta_c}, t~{t_c}", flush=True)
    window = initial_window
    doublings = 0
    while True:
        zeta_opt, rop_opt, t_opt, beta_opt, hit_boundary, n_fine_evals = fine_scan_warmstart(
            params, zeta_c, window, t_seed=t_c if t_c is not None else 0
        )
        if not hit_boundary or doublings >= max_doublings:
            break
        window *= 2
        doublings += 1
    elapsed = time.time() - t0
    return dict(
        zeta=zeta_opt, t=t_opt, beta=beta_opt, log2_rop=rop_opt,
        zeta_coarse=zeta_c, final_window=window, doublings=doublings,
        elapsed_s=elapsed, n_coarse_evals=len(coarse_results), n_fine_evals=n_fine_evals,
    )


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
    for name, params in PARAM_SETS.items():
        result = two_stage_search_warmstart(params)
        gt = GROUND_TRUTH[name]
        zeta_match = (result["zeta"] == gt["zeta"])
        rop_close = abs(result["log2_rop"] - gt["log2_rop"]) < 0.05
        match = "OK" if (zeta_match and rop_close) else "MISMATCH"
        print(f"{name}: zeta={result['zeta']} t={result['t']} log2(rop)={result['log2_rop']:.3f}  "
              f"GT zeta={gt['zeta']} GT log2(rop)={gt['log2_rop']}  match={match}  "
              f"time={result['elapsed_s']:.1f}s  fine_evals={result['n_fine_evals']}")


if __name__ == "__main__":
    main()

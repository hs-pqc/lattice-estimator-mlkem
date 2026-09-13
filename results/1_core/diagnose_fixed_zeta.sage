import estimator.schemes as schemes
from estimator import *
from estimator.lwe_dual import matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default, dual_hybrid


def matzov_cost_at_fixed_zeta(params, k_enum_val, red_cost_model=red_cost_model_default):
    """cost_fixed_zeta_fast에서 그대로 가져온 것: k_fft(t)와 p는 계속 최적화하되
    k_enum(zeta)만 강제로 고정한다."""
    best = None
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


def guo_johansson_at_fixed_zeta(params, zeta_val):
    """plain dual_hybrid, fft=True, zeta를 opt_step 탐색 범위를 1로 좁혀서 강제 고정에 가깝게."""
    r = dual_hybrid(params, zeta=zeta_val, fft=True)
    return r


targets = {
    "ML-KEM-768": (schemes.Kyber768, [15, 20, 23, 25, 30]),
    "Frodo640": (schemes.Frodo640, [5, 10, 13, 15, 20]),
}

for name, (params, zetas) in targets.items():
    print(f"\n=== {name}: MATZOV vs Guo-Johansson at fixed zeta ===")
    for z in zetas:
        try:
            m = matzov_cost_at_fixed_zeta(params, z)
            m_log2 = float(log(m["rop"], 2).n())
        except Exception as e:
            m_log2 = f"ERR: {e}"
        try:
            g = guo_johansson_at_fixed_zeta(params, z)
            g_log2 = float(log(g["rop"], 2).n())
        except Exception as e:
            g_log2 = f"ERR: {e}"
        print(f"  zeta={z:3d}  MATZOV={m_log2}   GuoJohansson={g_log2}")

from estimator import *
from estimator.lwe_dual import matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default

def cost_for_fixed_zeta_verbose(k_enum_val, params, red_cost_model=red_cost_model_default):
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

MLKEM_512 = LWE.Parameters(n=512, q=3329, Xs=ND.CenteredBinomial(3), Xe=ND.CenteredBinomial(3), tag="ML-KEM-512")
params = MLKEM_512.normalize()

print("zeta,log2_rop,p,t(k_fft),beta")
for zeta in range(0, 21, 1):
    try:
        r = cost_for_fixed_zeta_verbose(zeta, params)
        log2rop = float(log(r["rop"], 2).n())
        p_val = r.get("p", "?")
        t_val = r.get("t", "?")
        beta_val = r.get("beta", "?")
        print(f"{zeta},{log2rop:.3f},{p_val},{t_val},{beta_val}")
    except Exception as e:
        print(f"{zeta},ERROR:{e}")

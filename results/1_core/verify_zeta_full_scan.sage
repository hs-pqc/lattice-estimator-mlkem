from estimator import *
import estimator.lwe_dual as lwe_dual
from estimator.lwe_dual import matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default

MLKEM_768 = LWE.Parameters(n=3*256, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), m=3*256, tag="MLKEM_768")
params = MLKEM_768.normalize()

def cost_for_fixed_zeta(k_enum_val, params, red_cost_model=red_cost_model_default):
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

print("zeta,log2_rop")
results = []
step = 4
for zeta in range(0, params.n + 1, step):
    try:
        r = cost_for_fixed_zeta(zeta, params)
        log2rop = float(log(r["rop"], 2).n())
        results.append((zeta, log2rop))
        print(f"{zeta},{log2rop:.3f}")
    except Exception as e:
        print(f"{zeta},ERROR:{e}")

print()
print("=== summary ===")
best = min(results, key=lambda x: x[1])
print("global minimum found at zeta =", best[0], "log2(rop) =", best[1])
print("default LWE.dual_hybrid gave zeta=20 for comparison")

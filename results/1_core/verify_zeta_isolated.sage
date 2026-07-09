from estimator import *
import estimator.lwe_dual as lwe_dual
from estimator.lwe_dual import matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default

MLKEM_768 = LWE.Parameters(n=3*256, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), m=3*256, tag="MLKEM_768")

print("1. default (k_enum step=10, k_fft step=10)")
r_default = LWE.dual_hybrid(MLKEM_768)
print(r_default)

def custom_call(self, params, red_cost_model=red_cost_model_default, log_level=1, zeta_step=10):
    params = params.normalize()
    for p in early_abort_range(2, params.q):
        for k_enum in early_abort_range(0, params.n, zeta_step):
            for k_fft in early_abort_range(0, params.n - k_enum[0], 10):
                precision = 1
                max_beta = max(min(params.m - k_enum[0] - k_fft[0], max_beta_global), 40 + precision)
                with local_minimum(40, max_beta, precision=precision, log_level=log_level + 4) as it:
                    for beta in it:
                        cost = self.cost(beta, params, p=p[0], k_enum=k_enum[0], k_fft=k_fft[0], red_cost_model=red_cost_model)
                        it.update(cost)
                    k_fft[1].update(it.y)
            k_enum[1].update(k_fft[1].y)
        p[1].update(k_enum[1].y)
        if p[1].y["t"] == 0 and p[0] > 2:
            break
    return p[1].y

print("2. zeta(k_enum) step=10 -> 1 ONLY (k_fft stays at 10)")
r_fine = custom_call(matzov, MLKEM_768, zeta_step=1)
print(r_fine)

print("3. compare")
print("default zeta:", r_default["zeta"], "rop:", r_default["rop"])
print("fine    zeta:", r_fine["zeta"], "rop:", r_fine["rop"])
if r_default["zeta"] != r_fine["zeta"]:
    print(">>> zeta search resolution alone changes the result")
    if r_fine["rop"] < r_default["rop"]:
        print(">>> finer zeta grid found strictly better (lower) rop")
    else:
        print(">>> finer zeta grid found DIFFERENT but not better rop -- greedy path issue, not pure resolution issue")
else:
    print(">>> zeta unchanged when isolated -- original 20 vs 6 discrepancy came from k_fft interaction, not zeta step alone")

import time
from estimator import *
from estimator.lwe_dual import matzov, early_abort_range, local_minimum, max_beta_global, red_cost_model_default

params = LWE.Parameters(n=512, q=3329, Xs=ND.CenteredBinomial(3), Xe=ND.CenteredBinomial(3), tag="test").normalize()

def worker(zeta):
    p = 2
    max_beta = 60
    with local_minimum(40, max_beta, precision=1) as it:
        for beta in it:
            cost = matzov.cost(beta, params, p=p, k_enum=zeta, k_fft=20, red_cost_model=red_cost_model_default)
            it.update(cost)
    return (zeta, float(log(it.y["rop"], 2).n()))

if __name__ == "__main__":
    t0 = time.time()
    results = [worker(z) for z in range(8)]
    print("results (sequential):", results)
    print("elapsed:", time.time() - t0, "s")

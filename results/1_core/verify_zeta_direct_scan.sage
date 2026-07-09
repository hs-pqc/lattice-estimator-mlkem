from estimator import *
from estimator.lwe_dual import DH
from estimator.lwe_guess import exhaustive_search
from estimator.conf import red_cost_model as red_cost_model_default
from estimator.cost import Cost

Cost.register_impermanent(rop=True, mem=False, red=True, beta=False, delta=False, m=True, d=False, zeta=False, t=False)

MLKEM_768 = LWE.Parameters(n=768, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="ML-KEM-768")
params = MLKEM_768.normalize()

print("zeta,log2_rop")
results = []
step = 4
for zeta in range(0, params.n + 1, step):
    try:
        cost = DH.optimize_blocksize(solver=exhaustive_search, params=params, zeta=zeta, red_cost_model=red_cost_model_default, fft=False, log_level=5)
        log2rop = float(log(cost["rop"], 2).n())
        results.append((zeta, log2rop))
        print(f"{zeta},{log2rop:.3f}")
    except Exception as e:
        print(f"{zeta},ERROR:{e}")

print()
print("=== summary (direct call, no FFT) ===")
if results:
    best = min(results, key=lambda x: x[1])
    print("global minimum at zeta =", best[0], "log2(rop) =", best[1])
else:
    print("no results collected")
print("RESEARCH_NOTE.md recorded direct-call zeta=32, rop=206.4")

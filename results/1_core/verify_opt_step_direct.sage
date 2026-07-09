from estimator import *
from estimator.lwe_dual import dual_hybrid

MLKEM_768 = LWE.Parameters(n=768, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="ML-KEM-768")

print("opt_step,zeta,log2_rop")
for step in [1, 2, 4, 8, 16, 32]:
    r = dual_hybrid(MLKEM_768, opt_step=step, fft=False)
    log2rop = float(log(r["rop"], 2).n())
    zeta = int(r.get("zeta", -1))
    print(f"{step},{zeta},{log2rop:.3f}")

print()
print("reference (full scan, step=4): global min at zeta=32, log2(rop)=206.357")

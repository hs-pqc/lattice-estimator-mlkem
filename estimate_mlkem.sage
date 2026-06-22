# estimate_mlkem.sage
# ML-KEM security analysis using lattice-estimator
# Usage: sage estimate_mlkem.sage

from estimator import *
from estimator.nd import CenteredBinomial
import sys

# ── 파라미터 정의 ──────────────────────────────────────────

MLKEM_512 = LWE.Parameters(
    n=512, q=3329,
    Xs=ND.CenteredBinomial(3),
    Xe=ND.CenteredBinomial(3),
    tag="ML-KEM-512"
)

MLKEM_768 = LWE.Parameters(
    n=768, q=3329,
    Xs=ND.CenteredBinomial(2),
    Xe=ND.CenteredBinomial(2),
    tag="ML-KEM-768"
)

MLKEM_1024 = LWE.Parameters(
    n=1024, q=3329,
    Xs=ND.CenteredBinomial(2),
    Xe=ND.CenteredBinomial(2),
    tag="ML-KEM-1024"
)

PARAM_SETS = [MLKEM_512, MLKEM_768, MLKEM_1024]

# ── 유틸 ───────────────────────────────────────────────────

def print_section(title):
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)

def extract_cost(result):
    """결과에서 보안 비트와 β 추출"""
    try:
        rop = float(result.get("rop", 0))
        beta = int(result.get("beta", 0))
        return log(rop, 2).n(), beta
    except Exception:
        return None, None

# ── 1. Primal Attack (usvp) ────────────────────────────────

print_section("1. Primal Attack (usvp)")
print(f"{'Parameter':<15} {'n':<8} {'Security(bits)':<18} {'β':<8}")
print("-" * 50)

for params in PARAM_SETS:
    try:
        result = LWE.primal_usvp(params)
        bits, beta = extract_cost(result)
        if bits:
            print(f"{params.tag:<15} {params.n:<8} {bits:<18.1f} {beta:<8}")
    except Exception as e:
        print(f"{params.tag:<15} ERROR: {e}")

# ── 2. Dual Hybrid Attack ──────────────────────────────────

print_section("2. Dual Hybrid Attack")
print(f"{'Parameter':<15} {'Security(bits)':<18} {'β':<8} {'p':<6} {'ζ':<6}")
print("-" * 55)

for params in PARAM_SETS:
    try:
        result = LWE.dual_hybrid(params)
        bits, beta = extract_cost(result)
        p   = int(result.get("p", 0))
        zeta = int(result.get("zeta", 0))
        if bits:
            print(f"{params.tag:<15} {bits:<18.1f} {beta:<8} {p:<6} {zeta:<6}")
    except Exception as e:
        print(f"{params.tag:<15} ERROR: {e}")

# ── 3. Sparse Secret Analysis (ML-KEM-768) ─────────────────

print_section("3. Sparse Secret Analysis (ML-KEM-768)")
print("Varying Hamming weight h — dual hybrid")
print(f"{'h':<12} {'Security(bits)':<18} {'Δ from std':<14} {'ζ':<6}")
print("-" * 52)

STANDARD_BITS = 196.4  # 표준 파라미터 기준값

hamming_weights = [30, 50, 100, 200]

for h in hamming_weights:
    try:
        params_sparse = LWE.Parameters(
            n=768, q=3329,
            Xs=ND.Uniform(-1, 1, h=h),  # Hamming weight h
            Xe=ND.CenteredBinomial(2),
            tag=f"h={h}"
        )
        result = LWE.dual_hybrid(params_sparse)
        bits, beta = extract_cost(result)
        zeta = int(result.get("zeta", 0))
        if bits:
            delta = bits - STANDARD_BITS
            print(f"h={h:<10} {bits:<18.1f} {delta:<+14.1f} {zeta:<6}")
    except Exception as e:
        print(f"h={h:<10} ERROR: {e}")

# ── 4. ζ anomaly 분석 (ML-KEM-768) ────────────────────────

print_section("4. ML-KEM-768 ζ anomaly 상세 분석")
print("dual hybrid cost around ζ=20")
print(f"{'ζ':<8} {'Security(bits)':<18} {'β':<8}")
print("-" * 36)

for zeta_val in [0, 5, 10, 15, 20, 25, 30]:
    try:
        result = LWE.dual_hybrid(MLKEM_768, zeta=zeta_val)
        bits, beta = extract_cost(result)
        marker = " ← optimal" if zeta_val == 20 else ""
        if bits:
            print(f"ζ={zeta_val:<6} {bits:<18.1f} {beta:<8}{marker}")
    except Exception as e:
        print(f"ζ={zeta_val:<6} ERROR: {e}")

print("\n완료.")

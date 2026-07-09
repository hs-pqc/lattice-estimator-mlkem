# estimate_mlkem.sage
# ML-KEM security analysis using lattice-estimator
# Usage: sage estimate_mlkem.sage
#
# See RESEARCH_NOTE.md Appendix (2026-07) for the full derivation behind
# Sections 4-6 below (MATZOV zeta/t grid resolution mechanism).

from estimator import *
from estimator.nd import CenteredBinomial
from estimator.lwe_dual import matzov, dual_hybrid, early_abort_range, local_minimum, max_beta_global, red_cost_model_default
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

def cost_for_fixed_zeta_matzov(k_enum_val, params, red_cost_model=red_cost_model_default):
    """MATZOV cost model에서 zeta를 고정하고 p, t(k_fft), beta를 재최적화.
    early_abort_range의 그리디 조기종료를 우회해 특정 zeta 지점의
    실제 비용을 정확히 계산한다. (RESEARCH_NOTE.md Appendix Experiment 2/5 참고)
    """
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

# ── 2. Dual Hybrid Attack (MATZOV, official API) ───────────

print_section("2. Dual Hybrid Attack (LWE.dual_hybrid / MATZOV)")
print("NOTE: subject to the grid-resolution artifact documented in")
print("RESEARCH_NOTE.md Appendix. See Section 4 below for the corrected values.")
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

# ── 4. ζ 전수조사로 실제 최적값 계산 (MATZOV) ──────────────
# 이전 버전은 "marker = optimal if zeta_val == 20"으로 하드코딩되어
# 있었음 — 실제 최솟값을 계산한 게 아니라 가정을 표시만 한 것이었음.
# RESEARCH_NOTE.md Appendix Experiment 2에서 이게 부정확함을 확인
# (실제 최적은 ζ=24, 196.133bits, 기본 출력 ζ=20과 0.23bit 차이).
# 아래는 각 파라미터셋에 대해 실제로 argmin을 계산한다.

print_section("4. ζ 전수조사 — MATZOV 실제 최적값 (하드코딩 마커 제거)")

zeta_scan_ranges = {
    "ML-KEM-512": (0, 60, 2),
    "ML-KEM-768": (0, 40, 2),
    "ML-KEM-1024": (0, 60, 2),
}

for params in PARAM_SETS:
    zmin, zmax, zstep = zeta_scan_ranges[params.tag]
    print(f"\n--- {params.tag} (zeta={zmin}..{zmax}, step={zstep}) ---")
    print(f"{'ζ':<8} {'Security(bits)':<18} {'β':<8} {'p':<6} {'t':<6}")
    print("-" * 46)
    normalized = params.normalize()
    results = []
    default_result = LWE.dual_hybrid(params)
    default_zeta = int(default_result.get("zeta", 0))
    for zeta_val in range(zmin, zmax + 1, zstep):
        try:
            r = cost_for_fixed_zeta_matzov(zeta_val, normalized)
            log2rop = float(log(r["rop"], 2).n())
            results.append((zeta_val, log2rop, r))
        except Exception as e:
            print(f"ζ={zeta_val:<6} ERROR: {e}")
    if results:
        best_zeta, best_bits, best_r = min(results, key=lambda x: x[1])
        for zeta_val, bits, r in results:
            marker = " ← true optimum" if zeta_val == best_zeta else ""
            marker += " (default output)" if zeta_val == default_zeta else ""
            print(f"ζ={zeta_val:<6} {bits:<18.3f} {int(r['beta']):<8} {int(r.get('p',0)):<6} {int(r.get('t',0)):<6}{marker}")
        print(f"\n>>> {params.tag}: true optimum zeta={best_zeta} ({best_bits:.3f} bits), "
              f"default LWE.dual_hybrid reported zeta={default_zeta}")

print("\n완료. 상세 메커니즘 분석은 RESEARCH_NOTE.md Appendix 참고.")

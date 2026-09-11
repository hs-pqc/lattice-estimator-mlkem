"""
compare_dual_attacks.sage

Issue #219 후속 개인 연구 스켈레톤.

목적:
    동일 파라미터셋에 대해 세 가지 dual attack 비용 추정 경로를 비교한다.
        1) MATZOV indirect call  -> LWE.dual_hybrid(params)          (자체 그리드 서치)
        2) Guo-Johansson (WHT)   -> dual_hybrid(params, fft=True)    ([AC:GuoJoh21])
        3) plain dual_hybrid     -> dual_hybrid(params, fft=False)

    "어느 게 더 정확한가"가 아니라
    "동일 파라미터에서 어느 경로가 더 싼(=더 강한) 공격을 찾아내는가"를 스킴별/파라미터별로
    테이블화하는 것이 1차 목표. (Experiment 8에서 ML-KEM-768 단일 케이스만 확인했던 것을
    이미 검증된 9개 파라미터셋 전체로 확장)

상태: params/cost 함수 3종 모두 hs-pqc/lattice-estimator-mlkem 리포의 기존 스크립트
(verify_parallel_v2.sage, verify_dsa_all.sage, test_ntruplus.sage, verify_opt_step_direct.sage)
호출 시그니처를 그대로 이식해서 구현 완료. 1차 비교(compare_scheme_at_known_optimum)는
바로 실행 가능. 2차 비교(compare_scheme_own_optimum, 세 알고리즘을 "같은 zeta/t 지점"에서
강제 비교)는 아직 미구현 — 아래 함수 docstring 참고.

주의:
    - 세 함수 모두 "자체 검색"을 돌리는 구조라, zeta/t를 외부에서 고정해 넣는 공통
      인터페이스가 없다. 그래서 1차 비교는 "같은 지점에서 셋이 뭐라 하는가"가 아니라
      "각자 스스로 찾은 최적점이 얼마나 다른가"를 보는 것이다. 이 차이를 논문/노트에
      쓸 때 헷갈리지 않도록 주의.
    - ML-DSA/NTRU+의 known_optimum.log2_rop는 area 파일에 절대값이 기록돼 있지 않아
      None으로 남겨둠 — 이번 실행 결과가 그 기록값이 된다.
"""

from estimator import *
from estimator.lwe_dual import dual_hybrid
import json
import time
import os

# ---------------------------------------------------------------------------
# 1. 파라미터셋 정의
# ---------------------------------------------------------------------------
# hs-pqc/lattice-estimator-mlkem 리포의 기존 스크립트에서 그대로 이식:
#   - ML-KEM: results/1_core/verify_parallel_v2.sage
#   - ML-DSA: results/1_core/verify_dsa_all.sage
#   - NTRU+ : results/1_core/test_ntruplus.sage
#   - FrodoKEM: schemes 모듈 내장
#
# 각 항목의 "known_optimum"은 Issue #219 연구(Experiment 7/8, RESEARCH_NOTE.md)에서
# 이미 확인된 (ζ, t, log2(rop)) 값 — MATZOV(indirect) 기준. 이 스크립트의 1차 목표는
# 이 값을 재현하는 것(sanity check), 2차 목표는 세 알고리즘을 나란히 비교하는 것이다.

SCHEMES = {
    # --- ML-KEM (verify_parallel_v2.sage PARAM_SETS / GROUND_TRUTH 이식) ---
    "ML-KEM-512": {
        "params": LWE.Parameters(n=512, q=3329, Xs=ND.CenteredBinomial(3), Xe=ND.CenteredBinomial(3), tag="ML-KEM-512").normalize(),
        "known_optimum": {"zeta": 14, "t": 34, "log2_rop": 139.057},
    },
    "ML-KEM-768": {
        "params": LWE.Parameters(n=768, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="ML-KEM-768").normalize(),
        # Experiment 8에서 이미 확인된 값 (MATZOV indirect 기준)
        "known_optimum": {"zeta": 23, "t": 59, "log2_rop": 196.37},
    },
    "ML-KEM-1024": {
        "params": LWE.Parameters(n=1024, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="ML-KEM-1024").normalize(),
        "known_optimum": {"zeta": 32, "t": 82, "log2_rop": 261.143},
    },
    # --- ML-DSA (verify_dsa_all.sage PARAM_SETS 이식) ---
    "ML-DSA-44": {
        "params": LWE.Parameters(n=4 * 256, q=8380417, Xs=ND.UniformMod(2), Xe=ND.UniformMod(2), m=4 * 256, tag="ML-DSA-44").normalize(),
        "known_optimum": {"zeta": 23, "t": None, "log2_rop": None},
    },
    "ML-DSA-65": {
        "params": LWE.Parameters(n=5 * 256, q=8380417, Xs=ND.UniformMod(4), Xe=ND.UniformMod(4), m=6 * 256, tag="ML-DSA-65").normalize(),
        "known_optimum": {"zeta": 21, "t": None, "log2_rop": None},
    },
    "ML-DSA-87": {
        "params": LWE.Parameters(n=7 * 256, q=8380417, Xs=ND.UniformMod(2), Xe=ND.UniformMod(2), m=8 * 256, tag="ML-DSA-87").normalize(),
        "known_optimum": {"zeta": 13, "t": None, "log2_rop": None},
    },
    # --- NTRU+ (test_ntruplus.sage 이식) ---
    "NTRU+": {
        "params": LWE.Parameters(n=576, q=3457, Xs=ND.Binary, Xe=ND.SparseTernary(576, 192), m=576, tag="NTRU+576").normalize(),
        "known_optimum": {"zeta": 26, "t": None, "log2_rop": None},
    },
    # --- FrodoKEM (schemes 모듈 내장) ---
    "Frodo640": {
        "params": schemes.Frodo640,
        "known_optimum": {"zeta": 13, "t": 0, "log2_rop": None},
    },
    "Frodo976": {
        "params": schemes.Frodo976,
        # Frodo976은 전체 연구에서 가장 큰 gap이 나온 케이스 (default t=0 vs true t=43)
        "known_optimum": {"zeta": 17, "t": 43, "log2_rop": None},
    },
    "Frodo1344": {
        "params": schemes.Frodo1344,
        "known_optimum": {"zeta": 28, "t": 69, "log2_rop": None},
    },
}

# NOTE: ML-DSA-44/65/87, NTRU+의 known_optimum log2_rop는 아직 기록되어 있지 않음
# (area 파일에는 zeta/gap 비트수만 있고 절대 log2_rop 값이 없었음). None인 채로
# sanity check를 건너뛰고, 일단 세 알고리즘 간 상대 비교(compare_scheme_at_known_optimum)
# 부터 돌려서 채워도 됨 — 그 결과의 MATZOV 열이 곧 log2_rop 기록값이 됨.


# ---------------------------------------------------------------------------
# 2. 세 가지 비용 계산 경로
# ---------------------------------------------------------------------------
def cost_matzov_indirect(params, zeta=None, t=None):
    """
    MATZOV indirect call: LWE.dual_hybrid(params) — 내부적으로 matzov.cost()가 자체
    early_abort_range(step=10) 그리드 서치로 (zeta, t)를 찾아서 반환.
    Experiment 8에서 이 경로가 ML-KEM-768에 대해 Guo-Johansson보다 8.2비트 더 싼
    공격을 찾아낸 바 있음.

    NOTE: zeta/t 인자는 여기서 무시된다 — LWE.dual_hybrid는 zeta/t를 외부에서
    고정해서 넣는 인터페이스를 제공하지 않고, 항상 자체 검색을 돌린다 (이게 바로
    Issue #219의 원인이었던 hardcoded step=10 그리드 서치임). 그래서 이 함수는
    "known_optimum 지점에서의 비용"이 아니라 "이 알고리즘이 스스로 찾은 지점의 비용"을
    반환한다 — 세 알고리즘을 "같은 지점"에서 비교하고 싶다면 이 비대칭을 감안해야 함.
    """
    r = LWE.dual_hybrid(params)
    log2_rop = float(log(r["rop"], 2).n())
    return log2_rop, r


def cost_dual_hybrid_fft_true(params, zeta=None, t=None, opt_step=4):
    """
    Guo-Johansson (Walsh-Hadamard transform 기반) distinguisher.
    dual_hybrid(fft=True) 경로. MATZOV의 자체 FFT 기반 guessing cost model과는
    별개의 알고리즘([AC:GuoJoh21])이므로 혼동 주의 (RESEARCH_NOTE.md 정정 사항 참고).

    verify_opt_step_direct.sage에서 확인된 시그니처: dual_hybrid(params, opt_step=, fft=)
    이 direct 경로는 자체 zeta 검색을 돌리고 결과 dict에 "zeta"를 포함해서 반환한다.
    """
    r = dual_hybrid(params, opt_step=opt_step, fft=True)
    log2_rop = float(log(r["rop"], 2).n())
    return log2_rop, r


def cost_dual_hybrid_fft_false(params, zeta=None, t=None, opt_step=4):
    """
    plain dual_hybrid, FFT/WHT distinguisher 없이 exhaustive guessing만 사용하는 경로.
    """
    r = dual_hybrid(params, opt_step=opt_step, fft=False)
    log2_rop = float(log(r["rop"], 2).n())
    return log2_rop, r


ALGORITHMS = {
    "matzov_indirect": cost_matzov_indirect,
    "guo_johansson_fft_true": cost_dual_hybrid_fft_true,
    "dual_hybrid_fft_false": cost_dual_hybrid_fft_false,
}


# ---------------------------------------------------------------------------
# 3. 비교 실행
# ---------------------------------------------------------------------------
def compare_scheme_at_known_optimum(scheme_name, scheme_entry):
    """
    1차 단계: 세 알고리즘 각각을 스스로의 방식대로 돌려서 결과를 나란히 놓는다.

    중요: 이건 "같은 (zeta,t) 지점에서 세 알고리즘이 뭐라 하는가"가 아니다 — 세 함수
    모두 자체 검색을 돌리는 구조라 외부에서 zeta/t를 강제로 고정해 넣는 공통 인터페이스가
    없다 (MATZOV indirect가 특히 그렇다: Issue #219 자체가 이 자체 검색의 버그였음).
    그래서 여기서는 "각자 스스로 찾은 (zeta, log2_rop)"을 나란히 놓고,
    known_optimum(이미 검증된 참고값)과 얼마나 가까운지를 sanity check로 같이 표시한다.
    진짜 "같은 지점 강제 비교"가 필요하면 2차 단계(compare_scheme_own_optimum)에서
    matzov.cost()/DH.optimize_blocksize() 같은 저수준 함수를 직접 호출해 zeta/t를 고정해야 함
    (verify_zeta_direct_scan.sage 참고).
    """
    params = scheme_entry["params"]
    known = scheme_entry["known_optimum"]

    row = {"scheme": scheme_name, "known_optimum": known}
    for algo_name, algo_fn in ALGORITHMS.items():
        try:
            t0 = time.time()
            log2_rop, raw = algo_fn(params)
            elapsed = time.time() - t0
            row[algo_name] = {
                "log2_rop": round(log2_rop, 3),
                "zeta_found": int(raw.get("zeta")) if raw.get("zeta") is not None else None,
                "seconds": round(elapsed, 3),
            }
        except Exception as e:
            row[algo_name] = {"error": f"{type(e).__name__}: {e}"}
    return row


def compare_scheme_own_optimum(scheme_name, scheme_entry):
    """
    2차 단계 (TODO, 미구현): 각 알고리즘이 자기 자신의 (zeta, t) 검색 공간에서
    독립적으로 최적점을 찾도록 한 뒤 비교. Issue #219에서 밝혀진 것처럼
    step=10 하드코딩 같은 검색 구조 버그가 알고리즘 간 비교를 왜곡시킬 수 있으므로,
    여기서도 coarse-to-fine adaptive search(이미 검증된 verify_parallel_v2.sage 방식)를
    각 알고리즘에 동일하게 적용해야 공정한 비교가 됨.
    """
    raise NotImplementedError(
        "TODO: 각 알고리즘별 독립 (zeta,t) 탐색 — verify_parallel_v2.sage의 "
        "coarse-to-fine 로직을 세 알고리즘에 공통 적용"
    )


def main():
    results = []
    for scheme_name, scheme_entry in SCHEMES.items():
        row = compare_scheme_at_known_optimum(scheme_name, scheme_entry)
        if row is not None:
            results.append(row)

    out_dir = "results"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "dual_attack_comparison.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

        print(f"\n{len(results)}개 스킴 결과:")
    for row in results:
        print(f"\n{row['scheme']}:")
        for algo_name in ALGORITHMS:
            v = row[algo_name]
            if "error" in v:
                print(f"  {algo_name}: ERROR {v['error']}")
            else:
                print(f"  {algo_name}: zeta={v['zeta_found']}, log2(rop)={v['log2_rop']}, {v['seconds']}s")

    out_dir = "results"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "dual_attack_comparison.json")
           def _json_safe(o):
        try:
            return int(o)
        except (TypeError, ValueError):
            try:
                return float(o)
            except (TypeError, ValueError):
                return str(o)

    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=_json_safe)
    print(f"\n저장 완료: {out_path}")
if __name__ == "__main__":
    main()

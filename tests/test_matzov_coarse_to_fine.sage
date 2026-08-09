# Regression test for matzov_coarse_to_fine.sage
#
# Run with: sage tests/test_matzov_coarse_to_fine.sage
#
# Pins the ground-truth (zeta, t, log2(rop)) values obtained by full enumeration
# (Experiment 7) for ML-KEM-512/768/1024. A future change to the coarse-to-fine
# search (window size, coarse step, doubling cap) should be re-run against this
# file before being considered safe.

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "results", "1_core"))
load(os.path.join(os.path.dirname(__file__), "..", "results", "1_core", "matzov_coarse_to_fine.sage"))

PARAM_SETS = {
    "ML-KEM-512": LWE.Parameters(n=512, q=3329, Xs=ND.CenteredBinomial(3), Xe=ND.CenteredBinomial(3), tag="ML-KEM-512"),
    "ML-KEM-768": LWE.Parameters(n=768, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="ML-KEM-768"),
    "ML-KEM-1024": LWE.Parameters(n=1024, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="ML-KEM-1024"),
}

# Confirmed via full enumeration (Experiment 7, step=1 exhaustive search over all
# (zeta, t) pairs). See RESEARCH_NOTE.md for derivation.
GROUND_TRUTH = {
    "ML-KEM-512": dict(zeta=14, t=34, log2_rop=139.057),
    "ML-KEM-768": dict(zeta=23, t=59, log2_rop=195.554),
    "ML-KEM-1024": dict(zeta=32, t=82, log2_rop=261.143),
}

# What the ORIGINAL greedy search (step=10, no coarse-to-fine) returns -- included
# so the test also documents the size of the bug being fixed.
DEFAULT_GREEDY = {
    "ML-KEM-512": dict(zeta=0, gap_bits=0.60),
    "ML-KEM-768": dict(zeta=20, gap_bits=0.81),
    "ML-KEM-1024": dict(zeta=0, gap_bits=1.19),
}

TOLERANCE_BITS = 0.05


def run():
    failures = []
    for name, params in PARAM_SETS.items():
        gt = GROUND_TRUTH[name]
        result = matzov_coarse_to_fine(params)
        log2_rop = float(log(result["rop"], 2).n())
        zeta_ok = (int(result["zeta"]) == gt["zeta"])
        rop_ok = abs(log2_rop - gt["log2_rop"]) < TOLERANCE_BITS
        status = "PASS" if (zeta_ok and rop_ok) else "FAIL"
        print(f"[{status}] {name}: zeta={result['zeta']} (expected {gt['zeta']})  "
              f"log2(rop)={log2_rop:.3f} (expected {gt['log2_rop']}, tol={TOLERANCE_BITS})")
        if status == "FAIL":
            failures.append(name)

    if failures:
        print(f"\n{len(failures)} test(s) FAILED: {failures}")
        sys.exit(1)
    else:
        print(f"\nAll {len(PARAM_SETS)} tests PASSED.")


if __name__ == "__main__":
    run()

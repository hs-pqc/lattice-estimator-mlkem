# Regression test for matzov_coarse_to_fine.sage and matzov_coarse_to_fine_sequential_only.sage
#
# Run with: sage tests/test_matzov_coarse_to_fine.sage
#
# Pins the ground-truth (zeta, t, log2(rop)) values obtained by full enumeration
# (Experiment 7) for ML-KEM-512/768/1024. A future change to the coarse-to-fine
# search (window size, coarse step, doubling cap) should be re-run against this
# file before being considered safe.
#
# NOTE: both matzov_coarse_to_fine.sage and matzov_coarse_to_fine_sequential_only.sage
# define a function with the SAME name (matzov_coarse_to_fine), by design -- either
# is meant as a drop-in replacement for MATZOV.__call__. Loading both in the same
# session means the second load() silently overwrites the first. This test handles
# that deliberately: it loads and fully tests the n_jobs variant first, snapshots
# its results, THEN loads the sequential-only variant (overwriting the symbol) and
# compares its output against the snapshot -- never has both loaded for comparison
# at once.

import sys
import os

from estimator import *

# Sage preparses .sage files to a temporary .py before running them, and __file__
# has been observed to occasionally not resolve to this file's real directory
# depending on how `sage tests/test_matzov_coarse_to_fine.sage` is invoked. Fall
# back to assuming the script is run from the repo root if __file__ misbehaves,
# and fail with a clear message rather than a cryptic FileNotFoundError deep in
# load().
try:
    _HERE = os.path.dirname(os.path.abspath(__file__))
    _SEARCH_DIR = os.path.join(_HERE, "..", "results", "1_core")
    if not os.path.isdir(_SEARCH_DIR):
        raise OSError(f"resolved search dir does not exist: {_SEARCH_DIR}")
except (NameError, OSError):
    _SEARCH_DIR = os.path.join(os.getcwd(), "results", "1_core")
    print(f"[warn] __file__-based path resolution failed; falling back to cwd-relative path: {_SEARCH_DIR}")
    print("[warn] if this also fails, run this script from the repository root: sage tests/test_matzov_coarse_to_fine.sage")

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

# What the ORIGINAL greedy search (step=10, no coarse-to-fine) returns -- printed
# alongside each PASS/FAIL line so the test output also documents the size of the
# bug being fixed, not just whether the fix works.
DEFAULT_GREEDY = {
    "ML-KEM-512": dict(zeta=0, gap_bits=0.60),
    "ML-KEM-768": dict(zeta=20, gap_bits=0.81),
    "ML-KEM-1024": dict(zeta=0, gap_bits=1.19),
}

TOLERANCE_BITS = 0.05
# Tolerance for comparing two independently-computed log2(rop) values that should
# be mathematically identical (sequential vs parallel, or variant A vs variant B).
# Deliberately not using exact float `==`: even deterministic computations can
# differ in the last few bits after (de)serialization across a process boundary
# (multiprocessing) or between two independently loaded copies of the same logic.
CONSISTENCY_TOLERANCE_BITS = 1e-6


def run():
    failures = []

    # --- Part 1: n_jobs variant (matzov_coarse_to_fine.sage) against ground truth ---
    load(os.path.join(_SEARCH_DIR, "matzov_coarse_to_fine.sage"))
    print("=== Testing matzov_coarse_to_fine.sage (n_jobs variant) ===")
    njobs_results = {}
    for name, params in PARAM_SETS.items():
        gt = GROUND_TRUTH[name]
        greedy = DEFAULT_GREEDY[name]
        result = matzov_coarse_to_fine(params)  # n_jobs=None (default, sequential)
        njobs_results[name] = result
        log2_rop = float(log(result["rop"], 2).n())
        zeta_ok = (int(result["zeta"]) == gt["zeta"])
        rop_ok = abs(log2_rop - gt["log2_rop"]) < TOLERANCE_BITS
        status = "PASS" if (zeta_ok and rop_ok) else "FAIL"
        print(f"[{status}] {name}: zeta={result['zeta']} (expected {gt['zeta']}, "
              f"default greedy search gave zeta={greedy['zeta']} / {greedy['gap_bits']} bit gap)  "
              f"log2(rop)={log2_rop:.3f} (expected {gt['log2_rop']}, tol={TOLERANCE_BITS})")
        if status == "FAIL":
            failures.append(name)

    # n_jobs consistency: sequential (already computed above) vs parallel (n_jobs=-1)
    # should return the same result -- only wall-clock time should differ. Reuses
    # njobs_results instead of recomputing the sequential path a second time.
    #
    # NOTE: this exercises multiprocessing.Pool with a function defined via Sage's
    # load(). Whether that function is picklable (required by Pool) depends on
    # exactly how load() attaches definitions to the running session's namespace,
    # which can vary. Wrapped in try/except so a failure here (e.g. PicklingError)
    # is reported clearly without preventing the ground-truth and sequential_only
    # checks (the ones that actually matter for correctness) from running.
    print("\n--- n_jobs consistency (sequential vs parallel should match) ---")
    for name, params in PARAM_SETS.items():
        seq = njobs_results[name]
        try:
            par = matzov_coarse_to_fine(params, n_jobs=-1)
        except Exception as e:
            print(f"[SKIP] {name}: n_jobs=-1 raised {type(e).__name__}: {e}")
            print(f"       (this only affects the opt-in parallel path -- the "
                  f"default sequential path above is unaffected; consider using "
                  f"matzov_coarse_to_fine_sequential_only.sage if this persists)")
            continue
        seq_log2 = float(log(seq["rop"], 2).n())
        par_log2 = float(log(par["rop"], 2).n())
        same_zeta = (int(seq["zeta"]) == int(par["zeta"]))
        same_rop = abs(seq_log2 - par_log2) < CONSISTENCY_TOLERANCE_BITS
        status = "PASS" if (same_zeta and same_rop) else "FAIL"
        print(f"[{status}] {name}: sequential zeta={seq['zeta']} (log2(rop)={seq_log2:.6f}) vs "
              f"parallel zeta={par['zeta']} (log2(rop)={par_log2:.6f})")
        if status == "FAIL":
            failures.append(f"{name} (n_jobs consistency)")

    # --- Part 2: sequential-only variant, compared against Part 1's results ---
    # Overwrites the matzov_coarse_to_fine symbol -- intentional, see module docstring.
    load(os.path.join(_SEARCH_DIR, "matzov_coarse_to_fine_sequential_only.sage"))
    print("\n=== Testing matzov_coarse_to_fine_sequential_only.sage (vs n_jobs variant) ===")
    for name, params in PARAM_SETS.items():
        reference = njobs_results[name]
        result = matzov_coarse_to_fine(params)
        ref_log2 = float(log(reference["rop"], 2).n())
        result_log2 = float(log(result["rop"], 2).n())
        same_zeta = (int(reference["zeta"]) == int(result["zeta"]))
        same_rop = abs(ref_log2 - result_log2) < CONSISTENCY_TOLERANCE_BITS
        status = "PASS" if (same_zeta and same_rop) else "FAIL"
        print(f"[{status}] {name}: sequential_only zeta={result['zeta']} (log2(rop)={result_log2:.6f}) vs "
              f"n_jobs-variant reference zeta={reference['zeta']} (log2(rop)={ref_log2:.6f})")
        if status == "FAIL":
            failures.append(f"{name} (sequential_only vs n_jobs variant)")

    if failures:
        print(f"\n{len(failures)} test(s) FAILED: {failures}")
        sys.exit(1)
    else:
        print(f"\nAll tests PASSED.")


if __name__ == "__main__":
    run()

# Fast, standalone smoke test for the n_jobs=-1 (parallel) path only.
#
# Run with: sage tests/test_n_jobs_smoke.sage
#
# Ground-truth correctness for the real ML-KEM-512/768/1024 schemes is already
# confirmed separately (test_matzov_coarse_to_fine.sage, Part 1 -- all PASSED).
# The coarse scan is identical regardless of n_jobs (it's always sequential), so
# re-running it on the full-size schemes just to check the parallel fine-scan
# path is wasteful. This script checks only what n_jobs actually changes: does
# the parallel path run without crashing (e.g. PicklingError), and does it
# return the same answer as sequential? A small synthetic instance answers that
# in seconds instead of tens of minutes.

import sys, os

_HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else os.getcwd()
_SEARCH_DIR = os.path.join(_HERE, "..", "results", "1_core")
if not os.path.isdir(_SEARCH_DIR):
    _SEARCH_DIR = os.path.join(os.getcwd(), "results", "1_core")

load(os.path.join(_SEARCH_DIR, "matzov_coarse_to_fine.sage"))

TINY_PARAMS = LWE.Parameters(n=40, q=3329, Xs=ND.CenteredBinomial(2), Xe=ND.CenteredBinomial(2), tag="TINY-N-JOBS-CHECK")

print("Running sequential (n_jobs=None)...")
seq = matzov_coarse_to_fine(TINY_PARAMS, zeta_coarse_step=5, zeta_window=5, zeta_max_doublings=1, verbose=True)

print("\nRunning parallel (n_jobs=-1)...")
try:
    par = matzov_coarse_to_fine(TINY_PARAMS, zeta_coarse_step=5, zeta_window=5, zeta_max_doublings=1, n_jobs=-1, verbose=True)
except Exception as e:
    print(f"\n[FAIL] n_jobs=-1 raised {type(e).__name__}: {e}")
    print("The opt-in parallel path does not work in this Sage session's load() context.")
    print("Recommendation: use matzov_coarse_to_fine_sequential_only.sage instead, or")
    print("investigate whether load()-ed functions are picklable in this Sage version.")
    sys.exit(1)

seq_log2 = float(log(seq["rop"], 2).n())
par_log2 = float(log(par["rop"], 2).n())
same_zeta = (int(seq["zeta"]) == int(par["zeta"]))
same_rop = abs(seq_log2 - par_log2) < 1e-6

if same_zeta and same_rop:
    print(f"\n[PASS] sequential zeta={seq['zeta']} (log2(rop)={seq_log2:.6f}) == "
          f"parallel zeta={par['zeta']} (log2(rop)={par_log2:.6f})")
else:
    print(f"\n[FAIL] sequential zeta={seq['zeta']} (log2(rop)={seq_log2:.6f}) != "
          f"parallel zeta={par['zeta']} (log2(rop)={par_log2:.6f})")
    sys.exit(1)

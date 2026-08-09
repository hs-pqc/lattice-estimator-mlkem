# Fix MATZOV's ζ/t search missing the true optimum (Issue #219)

## Problem

`MATZOV.__call__` (aliased as `LWE.dual_hybrid`) searches over the guessing
dimension `ζ` and FFT dimension `t` using `early_abort_range(..., step=10)` — a
fixed-step, greedy, early-aborting search. The hardcoded `step=10` is not exposed
to the caller.

The `(ζ, t) → log2(rop)` cost surface is **non-convex at this resolution**: because
the optimal BKZ block size `β` is integer-valued, small changes in `ζ` or `t` can
shift the optimal `β` by 1, producing a sawtooth-shaped cost surface. A fixed-step
greedy search over a sawtooth surface can converge to a local artifact rather than
the true global optimum — and in several observed cases the search's coarse
estimate lands exactly on a boundary (`ζ=0`), compounding the problem.

## Impact (confirmed via full enumeration, step=1, all `(ζ, t)` pairs)

| Parameter set | Default (greedy) ζ | True optimal ζ | Gap |
|---|---|---|---|
| ML-KEM-512 | 0 | 14 | 0.60 bits |
| ML-KEM-768 | 20 | 23 | 0.81 bits |
| ML-KEM-1024 | 0 | 32 | 1.19 bits |
| ML-DSA-44 | 10 | 23 | 0.345 bits |
| ML-DSA-65 | 20 | 21 | 0.464 bits |
| ML-DSA-87 | 0 | 13 | 0.266 bits |
| NTRU+ (KpqC) | 20 | 26 | 0.128 bits |
| FrodoKEM-640 | ζ=10,t=0 | ζ=13,t=0 | 0.821 bits |
| **FrodoKEM-976** | **ζ=20,t=0** | **ζ=17,t=43** | **7.981 bits** |
| FrodoKEM-1344 | ζ=30,t=60 | ζ=28,t=69 | 0.468 bits |

The FrodoKEM-976 case is qualitatively different: the default search concludes the
FFT distinguisher (`t>0`) isn't worth using at all (`t=0`), when in fact `t=43` is
far stronger — a ~253x cost ratio.

Note: this bug is structural to any scheme routed through `LWE.dual_hybrid`/
`MATZOV` — confirmed to affect ML-KEM, ML-DSA, FrodoKEM, and NTRU+ (KpqC). NIST's
standard NTRU (HPS/HRSS) is **not** affected, since `NTRU.dual_hybrid` is a
different function object.

## Fix

Replace the single fixed-step greedy search with an **adaptive coarse-to-fine**
search (`matzov_coarse_to_fine.sage` in this PR), applied independently to both
`ζ` and `t`:

1. **Coarse pass**: evaluate at `step=10` across the full range (same cost as
   today's search).
2. **Fine pass**: exhaustively evaluate every integer point in a window around the
   coarse optimum (default `window=15`).
3. **Boundary safeguard**: if the fine-pass optimum lands on the window's edge
   (signal that the window was too narrow), double the window and repeat, up to
   3 doublings. This is what correctly recovers `ζ=14/23/32` even when the coarse
   pass lands at the `ζ=0` boundary.

This is verified to reproduce the full-enumeration ground truth for all three
ML-KEM parameter sets to within 0.05 bits (see `tests/test_matzov_coarse_to_fine.sage`).

## Performance

Measured on a 32-core desktop (results/1_core/verify_parallel_v2.sage):

| | Full enumeration (step=1) | Coarse-to-fine, sequential | Coarse-to-fine, parallel (32-core) |
|---|---|---|---|
| ML-KEM-512 | ~6201s | 1473.8s | 97.2s |
| ML-KEM-768 | ~7001s | 2400.4s | 123.0s |
| ML-KEM-1024 | ~11658s | 2904.3s | 155.7s |
| **Combined** | **~24840s (6.9h)** | **6784.8s (113m)** | **376.9s (6m17s)** |

Even single-threaded, coarse-to-fine is **3.66x** faster than full enumeration
while matching it exactly.

Parallelism is **opt-in** via a `n_jobs` parameter (default `None` = fully
sequential, zero new dependency, zero process-pool overhead). Passing
`n_jobs=-1` (all cores) or an explicit integer parallelizes the fine zeta-scan
via `multiprocessing.Pool` and reproduces the 66x figure above. The coarse pass
always runs sequentially regardless of `n_jobs` (it's a small fraction of total
cost — see the coarse/fine timing breakdown in this repo's research notes), so
`n_jobs=None` and `n_jobs>1` are guaranteed to return bit-identical results;
`n_jobs` only changes wall-clock time, never the answer. This keeps the default
path dependency-free for library use while giving callers who want the speedup
a one-line way to opt in.

## What was tried and did not work

A warm-start variant (reusing the neighboring ζ's optimal `t` as a search hint,
to skip most of the coarse `t`-sweep) was explored to reduce runtime further. It
was measured and abandoned:

- On ML-KEM-512, the safeguarded warm-start variant took **1532.9s**, slightly
  *slower* than the plain coarse-to-fine baseline (**1473.8s**, no hint).
- The warm-start hint triggered its boundary-fallback safeguard (falling back to
  the exact full coarse sweep) on **16/16 (100%)** of the ζ values in the fine
  window — meaning the hint was never actually usable, and every ζ paid for both
  the failed hint attempt *and* the full fallback sweep.

The mechanism: near sawtooth-heavy regions of the cost surface (which, per the
Impact section, is most of it), a warm-started narrow window is reliably misled
by a local artifact into a false boundary hit. The neighboring-ζ optimal-`t`
values were in fact close to each other (e.g. ζ=0→t=50, ζ=1→t=49), confirming
the hint's *direction* was correct — but sawtooth noise near the hint still
produced a spurious boundary reading often enough to force fallback essentially
every time, at 100% measured on this parameter set.

We consider plain coarse-to-fine (without warm-starting) the right
accuracy/complexity trade-off for a security-estimation tool: the safeguard that
makes warm-start safe is exactly what makes it not actually save any work.

## Files in this PR

- `results/1_core/matzov_coarse_to_fine.sage` — the core fix. No new dependency
  on the default path (`n_jobs=None`); `multiprocessing` is only imported if the
  caller opts into `n_jobs`.
- `results/1_core/matzov_coarse_to_fine_sequential_only.sage` — functionally
  identical alternative with zero multiprocessing code anywhere in the file, not
  even opt-in. Offered in case the maintainer prefers no process-pool code in
  the core library on principle. Pick one, not both — see the note at the top of
  `tests/test_matzov_coarse_to_fine.sage` about why they share a function name.
- `tests/test_matzov_coarse_to_fine.sage` — regression test against confirmed
  ground truth, covering both files above and the n_jobs sequential/parallel
  consistency check
- `results/1_core/verify_dsa_all.sage`, `verify_frodo_all.sage`,
  `test_ntruplus.sage` — cross-scheme validation scripts

from estimator import *
import inspect

print("=== LWE.dual_hybrid signature ===")
print(inspect.signature(LWE.dual_hybrid))

print()
print("=== plain dual_hybrid signature ===")
print(inspect.signature(dual_hybrid))

print()
print("=== matzov module contents ===")
import estimator.reduction as red
try:
    from estimator.cost import MATZOV
    print(inspect.signature(MATZOV.__call__))
except Exception as e:
    print("MATZOV import/signature failed:", e)

try:
    from estimator import LWE as LWEmod
    print("LWE.dual_hybrid.__wrapped__?", getattr(LWE.dual_hybrid, '__wrapped__', None))
except Exception as e:
    print(e)

"""Smoke-тест Шага 4: коэффициенты распределения ткань/плазма.

Запуск из корня проекта: ``python -m scripts.smoke_test_step4``.
"""

from src.pbpk import k_liver, k_kidney, k_brain

# Reference values, must match the formula.
assert abs(k_liver(1.0) - 10**0.6) < 1e-9
assert abs(k_liver(2.0) - 10**1.1) < 1e-9

# Monotonicity.
assert k_liver(3.0) > k_liver(1.0)              # more lipophilic -> higher K
assert k_brain(2.0, 20) > k_brain(2.0, 80)      # higher TPSA -> lower K

print("[OK] partition.py работает")
print(f"  k_liver(logp=1.0)  = {k_liver(1.0):.3f}")
print(f"  k_kidney(logp=2, tpsa=60) = {k_kidney(2.0, 60.0):.3f}")
print(f"  k_brain(logp=2, tpsa=60)  = {k_brain(2.0, 60.0):.4f}")
"""Smoke-тест Шага 5: проницаемость гематоэнцефалического барьера.

Запуск из корня проекта: ``python -m scripts.smoke_test_step5``.
"""

from src.pbpk import permeability

# Neutral case: score = 0.5 -> sig ~ 0.622.
p_neutral = permeability(2.0, 60.0, 1.0)
print(f"neutral: {p_neutral:.4f}")   # ~ 0.622
assert 0.60 < p_neutral < 0.65

# Lipophilic molecule: score = 4.5 -> sig ~ 0.989 (saturation).
p_lipophilic = permeability(4.0, 20.0, 1.0)
print(f"lipophilic: {p_lipophilic:.4f}")  # ~ 0.989
assert p_lipophilic > 0.9

# Polar molecule: score = -4.5 -> sig ~ 0.011 (close to zero, as expected).
p_polar = permeability(0.0, 120.0, 1.0)
print(f"polar: {p_polar:.4f}")  # ~ 0.011
assert 0.005 < p_polar < 0.02

# Damaged BBB: permeability is exactly halved.
p_damaged = permeability(2.0, 60.0, 0.5)
assert abs(p_damaged - p_neutral * 0.5) < 1e-12
print(f"damaged BBB (integrity=0.5): {p_damaged:.4f}")

# Out-of-range integrity must raise.
try:
    permeability(2.0, 60.0, 1.5)
    print("[FAIL] integrity=1.5 прошёл!")
except ValueError as e:
    print(f"[OK] integrity=1.5 → {e}")

print("[OK] bbb.py работает")
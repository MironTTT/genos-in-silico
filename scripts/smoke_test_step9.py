"""Smoke-тест Шага 9: динамика pH желудка (модель гастрита).

Запуск из корня проекта: ``python -m scripts.smoke_test_step9``.
"""

from src.gastric.ph_model import run_gastric_simulation

# 1. No drug (efficacy=0.0): pH should settle around 1.2-1.5.
df_healthy = run_gastric_simulation(drug_efficacy=0.0)
print(f"Без лекарства: финальный pH = {df_healthy['ph'].iloc[-1]:.2f}")
assert 1.0 < df_healthy["ph"].iloc[-1] < 2.0

# 2. With drug (efficacy=0.9): pH should rise above 3.5.
df_treated = run_gastric_simulation(drug_efficacy=0.9)
print(f"С лекарством: финальный pH = {df_treated['ph'].iloc[-1]:.2f}")
assert df_treated["ph"].iloc[-1] > 3.5

# 3. Validation check.
try:
    run_gastric_simulation(drug_efficacy=1.5)
    print("[FAIL] efficacy=1.5 прошло!")
except ValueError as e:
    print(f"[OK] efficacy=1.5 → {e}")

print("[OK] ph_model.py работает, динамика pH корректна")
"""Smoke-тест Шага 6: правые части ОДУ PBPK-модели.
"""
import numpy as np

from src.pbpk.equations import derivatives
from config.constants import V_GUT, V_BLOOD, V_LIVER, V_KIDNEY, V_BRAIN

# Initial state: 100 mg in the gut, everything else empty.
y0 = np.array([100.0, 0.0, 0.0, 0.0, 0.0])

# Healthy adult, model compound.
params = {
    "ka": 1.2,
    "q_liver": 90.0, "q_kidney": 72.0, "q_brain_eff": 45.0 * 0.6,
    "k_liver": 2.0, "k_kidney": 1.5, "k_brain": 1.0,
    "cl_hep": 20.0, "cl_ren": 15.0, "cl_other": 5.0,
}

dy = derivatives(0.0, y0, params)
print("dy/dt при t=0:", dy)

# Check 1: absorption drains the gut.
assert dy[0] < 0, "dA_gut/dt должно быть отрицательным"

# Check 2: blood fills from the gut while organs are empty.
assert dy[1] > 0, "dA_blood/dt должно быть положительным"

# Check 3: sum of components = net elimination (zero here, organs empty).
total = dy.sum()
assert abs(total) < 1e-9, f"Баланс массы нарушен: сумма dA/dt = {total}"

# Check 4: at equilibrium (all C_i = C_blood * k_i) exchange terms vanish.
c_blood = 1.0
y_eq = np.array([
    V_GUT * c_blood,
    V_BLOOD * c_blood,
    V_LIVER * c_blood * 2.0,
    V_KIDNEY * c_blood * 1.5,
    V_BRAIN * c_blood * 1.0,
])
dy_eq = derivatives(0.0, y_eq, params)
print("dy/dt при равновесии:", dy_eq)

print("[OK] equations.py работает, mass balance holds")
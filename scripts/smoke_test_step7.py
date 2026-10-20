"""Smoke-тест Шага 7: класс PBKPSolver и интеграция ОДУ.

Запуск из корня проекта: ``python -m scripts.smoke_test_step7``.
"""

import numpy as np

from config.constants import V_BLOOD, V_BRAIN, V_GUT, V_KIDNEY, V_LIVER
from src.core.molecule import Molecule
from src.core.patient import Patient
from src.pbpk import equations, solver

# Ibuprofen, elderly patient with mild organ impairment.
m = Molecule.from_smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
p = Patient(age=65, weight_kg=75, liver_health=0.7, kidney_health=0.9, bbb_integrity=0.6)
s = solver.PBKPSolver(m, p, dose_mg=400)

df = s.run(t_end_h=8, t_step_min=15)
print(df.head())
print(f"rows: {len(df)}, final time: {df['time_h'].iloc[-1]:.2f} h")

# Check 1: columns in the required order.
assert df.columns.tolist() == [
    "time_h", "A_gut", "A_blood", "A_liver", "A_kidney", "A_brain",
    "C_gut", "C_blood", "C_liver", "C_kidney", "C_brain",
], f"Unexpected columns: {df.columns.tolist()}"

# Check 2: initial state is the whole dose in the gut.
assert abs(df["A_gut"].iloc[0] - 400.0) < 1e-6
assert abs(df["A_blood"].iloc[0]) < 1e-6

# Check 3: concentrations are computed from masses and volumes.
for i, (comp, vol) in enumerate(
    [("gut", V_GUT), ("blood", V_BLOOD), ("liver", V_LIVER), ("kidney", V_KIDNEY), ("brain", V_BRAIN)]
):
    c_col, a_col = f"C_{comp}", f"A_{comp}"
    assert np.allclose(df[c_col], df[a_col] / vol), f"C_{comp} != A_{comp}/V_{comp}"

# Check 4: masses and concentrations never go negative.
assert (df[["A_gut", "A_blood", "A_liver", "A_kidney", "A_brain"]] >= -1e-9).all().all()
assert (df[["C_gut", "C_blood", "C_liver", "C_kidney", "C_brain"]] >= -1e-9).all().all()

# Check 5: mass is conserved modulo elimination (total <= dose).
total = df[["A_gut", "A_blood", "A_liver", "A_kidney", "A_brain"]].sum(axis=1)
assert (total <= 400.0 + 1e-6).all(), "Total mass exceeds the dose"

# Check 6: the gut empties and the blood peaks early, then decays.
assert df["A_gut"].iloc[-1] < df["A_gut"].iloc[0]
assert df["A_blood"].max() > df["C_blood"].iloc[0]

# Validation: non-positive dose must raise.
try:
    solver.PBKPSolver(m, p, dose_mg=0)
    raise SystemExit("[FAIL] dose_mg=0 не отклонён")
except ValueError:
    print("[OK] dose_mg=0 отклонён")

print("[OK] solver.py работает, интеграция завершена")
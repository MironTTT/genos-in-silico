"""Smoke-тест Шага 10: связка PBPK и pH-моделей (GastricSimulation).

Запуск из корня проекта: ``python -m scripts.smoke_test_step10``.
"""

from src.core.molecule import Molecule
from src.core.patient import Patient
from src.orchestrator.gastric_pipeline import GastricSimulation

# Ibuprofen, elderly patient with mild organ impairment.
_IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
_PATIENT = Patient(
    age=65, weight_kg=75, liver_health=0.7, kidney_health=0.9, bbb_integrity=0.6
)

# 1. Column contract and bounded efficacy.
m_strong = Molecule.from_smiles(_IBUPROFEN, ki=0.1)
df = GastricSimulation(m_strong, _PATIENT, dose_mg=400.0).run(t_end_h=24.0)
assert df.columns.tolist() == ["time_h", "C_blood", "efficacy", "acid_mmol", "ph"]
assert (df["efficacy"] >= 0.0).all() and (df["efficacy"] <= 1.0).all()
assert (df["acid_mmol"] >= 0.0).all()

# 2. A weak inhibitor (large ki) barely shifts pH off the baseline.
m_weak = Molecule.from_smiles(_IBUPROFEN, ki=100.0)
df_weak = GastricSimulation(m_weak, _PATIENT, dose_mg=400.0).run(t_end_h=24.0)
print(f"Слабый ингибитор (ki=100): финальный pH = {df_weak['ph'].iloc[-1]:.2f}, "
      f"пик pH = {df_weak['ph'].max():.2f}")
assert df_weak["ph"].iloc[-1] < 1.15

# 3. A strong inhibitor (small ki) raises pH while the drug is in blood.
print(f"Сильный ингибитор (ki=0.1): финальный pH = {df['ph'].iloc[-1]:.2f}, "
      f"пик pH = {df['ph'].max():.2f}")
assert df["ph"].max() > df_weak["ph"].max() + 1.0

# 4. ki must be strictly positive.
try:
    Molecule.from_smiles(_IBUPROFEN, ki=0.0)
    print("[FAIL] ki=0 прошло!")
except ValueError as e:
    print(f"[OK] ki=0 → {e}")

print("[OK] gastric_pipeline.py работает, PBPK и pH связаны")
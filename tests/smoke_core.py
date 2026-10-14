"""Smoke-тест ядра PBPK-симулятора (модули molecule и patient).

Запуск: ``python tests/smoke_core.py`` из корня проекта.
При успехе всех проверок печатается ``ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ``.
"""

import sys
from pathlib import Path

# Корень проекта добавляется в sys.path, чтобы скрипт работал
# при запуске из любого каталога: python tests/smoke_core.py
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.core.molecule import Molecule
from src.core.patient import Patient


def expect_value_error(label: str, fn) -> None:
    """Проверяет, что переданная функция бросает ValueError.

    :param label: описание проверки для сообщения об ошибке.
    :param fn: функция без аргументов.
    :raises SystemExit: если ValueError не поднят.
    """
    try:
        fn()
    except ValueError:
        return
    raise SystemExit(f"ОШИБКА: {label} не вызвал ValueError")


# --- Molecule ---
ibu = Molecule.from_smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
assert ibu.smiles == "CC(C)Cc1ccc(C(C)C(=O)O)cc1"  # каноническая форма RDKit
assert ibu.mw > 100 and isinstance(ibu.logp, float) and ibu.tpsa > 0
assert isinstance(ibu.solubility_logs, float)
expect_value_error("невалидный SMILES", lambda: Molecule.from_smiles("C1CC1C(INVALID"))
expect_value_error("пустой SMILES", lambda: Molecule.from_smiles("   "))
try:
    ibu.mw = 1.0  # noqa: SLF001 - проверка неизменяемости
    raise SystemExit("ОШИБКА: Molecule изменяема")
except Exception:
    pass

# --- Patient: базовое создание и границы ---
p = Patient(age=25, weight_kg=70)
assert p.liver_health == 1.0 and p.cardiac_output == 6.5
assert Patient(age=120, weight_kg=300, cardiac_output=15.0)  # верхние границы
assert Patient(age=0, weight_kg=1.0, cardiac_output=1.0)  # нижние границы

# --- 5 строк по заданию ---
assert Patient(age=25, weight_kg=70)  # OK
expect_value_error("age=150", lambda: Patient(age=150, weight_kg=70))
expect_value_error("weight_kg=400", lambda: Patient(age=25, weight_kg=400))
expect_value_error("cardiac_output=20", lambda: Patient(age=25, weight_kg=70, cardiac_output=20))

p = Patient(age=25, weight_kg=70)
p2 = p.copy_with(liver_health=0.3)
assert p2 is not p and p2.liver_health == 0.3
assert p.liver_health == 1.0  # исходный p не изменился

# --- to_dict / copy_with ---
d = p.to_dict()
assert d == {"age": 25, "weight_kg": 70.0, "liver_health": 1.0,
             "kidney_health": 1.0, "bbb_integrity": 1.0, "cardiac_output": 6.5}
assert all(isinstance(v, (float, int)) for v in d.values())
p3 = p.copy_with(age=60, kidney_health=0.5, cardiac_output=5.0)
assert p3.age == 60 and p3.kidney_health == 0.5 and p3.cardiac_output == 5.0
try:
    p.copy_with(weight_kg=400)
    raise SystemExit("ОШИБКА: copy_with допустил вес 400")
except ValueError:
    pass

print("ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ")
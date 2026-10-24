"""Pytest-тесты сохранения массы в PBPK-модели.

Проверяется ключевой инвариант системы ОДУ: при отключённых клиренсах
(cl_hep = cl_ren = cl_other = 0) обменные потоки между кровью и органами
взаимно сокращаются, поэтому общая масса препарата обязана строго
сохраняться на всём интервале интегрирования. Тест использует реальный
пайплайн ядра (``Molecule``, ``Patient``, ``PBKPSolver``, ``derivatives``).

Запуск из корня проекта: ``python -m pytest tests/test_pbpk_massbalance.py -v``.
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from src.core.molecule import Molecule
from src.core.patient import Patient
from src.pbpk.equations import derivatives
from src.pbpk.solver import PBKPSolver

# Reference compound: ibuprofen.
_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
_DOSE_MG = 400.0
_T_END_H = 24.0
_T_STEP_H = 0.25


def _zero_clearance_params() -> dict:
    """Параметры модели с обнулёнными клиренсами (замкнутая система).

    Клиренсы зануляются через приватный ``_build_params`` — это часть
    контракта теста: проверить инвариант на реальном наборе параметров.
    """
    molecule = Molecule.from_smiles(_SMILES)
    patient = Patient(
        age=65, weight_kg=75, liver_health=0.7, kidney_health=0.9, bbb_integrity=0.6
    )
    params = PBKPSolver(molecule, patient, dose_mg=_DOSE_MG)._build_params()
    params["cl_hep"] = 0.0
    params["cl_ren"] = 0.0
    params["cl_other"] = 0.0
    return params


def _integrate(params: dict) -> tuple[np.ndarray, np.ndarray]:
    """Интегрирует модель из состояния «вся доза в ЖКТ».

    :param params: параметры модели для ``derivatives``.
    :return: (временные точки, матрица масс формата (5, N)).
    """
    y0 = np.array([_DOSE_MG, 0.0, 0.0, 0.0, 0.0])
    t_eval = np.arange(0.0, _T_END_H + _T_STEP_H, _T_STEP_H)
    sol = solve_ivp(
        derivatives,
        t_span=(0.0, _T_END_H),
        y0=y0,
        t_eval=t_eval,
        args=(params,),
        method="LSODA",
        rtol=1e-6,
        atol=1e-9,
    )
    assert sol.success, sol.message
    return sol.t, sol.y


def test_total_mass_is_conserved_without_clearances() -> None:
    """Доза целиком остаётся в системе: mass(t) == dose на всех точках."""
    _, masses = _integrate(_zero_clearance_params())
    total = masses.sum(axis=0)
    np.testing.assert_allclose(total, _DOSE_MG, rtol=1e-9, atol=1e-8)


def test_derivative_sum_is_zero_without_clearances() -> None:
    """Сумма правых частей ОДУ равна нулю: потоки взаимно сокращаются."""
    params = _zero_clearance_params()
    y = np.array([100.0, 20.0, 10.0, 5.0, 2.0])
    dy = derivatives(0.0, y, params)
    assert dy.sum() == pytest.approx(0.0, abs=1e-12)


def test_mass_decays_when_clearances_are_active() -> None:
    """Анти-тест: с включёнными клиренсами масса обязана убывать.

    Гарантирует, что первый тест не тривиален: сохранение массы — это
    следствие именно нулевых клиренсов, а не ошибки в интеграции.
    """
    molecule = Molecule.from_smiles(_SMILES)
    patient = Patient(
        age=65, weight_kg=75, liver_health=0.7, kidney_health=0.9, bbb_integrity=0.6
    )
    params = PBKPSolver(molecule, patient, dose_mg=_DOSE_MG)._build_params()
    _, masses = _integrate(params)
    assert masses.sum(axis=0)[-1] < _DOSE_MG
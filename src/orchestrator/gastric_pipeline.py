"""Связка PBPK- и pH-моделей: сценарий «лекарство и кислотность».

Модуль берёт кривую концентрации препарата в крови из PBPK-симуляции
(``src.pbpk.solver.PBKPSolver``), пересчитывает её в эффективность
ингибирования протонного насоса по модели насыщения
``efficacy = C_blood / (ki + C_blood)`` и подаёт зависящую от времени
эффективность в модель кислотности желудка (``src.gastric.ph_model``).

Пакет зависит и от Блока 2, и от Блока 3; сами модели не изменяются.
"""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from config.constants import (
    GASTRIC_EVACUATION_RATE,
    GASTRIC_INITIAL_ACID_MMOL,
    GASTRIC_NEUTRALIZATION,
    GASTRIC_SECRETION_BASE,
    V_GUT,
)
from src.core.molecule import Molecule
from src.core.patient import Patient
from src.gastric.ph_model import calculate_ph, gastric_acid_derivatives
from src.pbpk.solver import PBKPSolver

# Output time step for both models, min (matches the PH model default).
_T_STEP_MIN = 15.0


def _gastric_rhs_with_dynamic_efficacy(
    t: float,
    y: np.ndarray,
    params: dict,
    t_grid: np.ndarray,
    efficacy_grid: np.ndarray,
) -> np.ndarray:
    """Правая часть pH-ОДУ с эффективностью, интерполированной по времени.

    Функция-обёртка над ``gastric_acid_derivatives``: вместо фиксированной
    ``drug_efficacy`` из ``params`` подставляется значение из массива
    ``efficacy_grid``, интерполированное в текущий момент ``t`` (np.interp).

    :param t: время, ч.
    :param y: массив [A] — количество кислоты, ммоль.
    :param params: базовый словарь параметров кислотности.
    :param t_grid: временные точки, соответствующие ``efficacy_grid``.
    :param efficacy_grid: эффективность ингибирования в точках ``t_grid``.
    :return: массив ``[dA_dt]``, форма (1,).
    """
    live_params = dict(params)
    live_params["drug_efficacy"] = float(np.interp(t, t_grid, efficacy_grid))
    return gastric_acid_derivatives(t, y, live_params)


class GastricSimulation:
    """Полный сценарий: PBPK-фармакокинетика + pH желудка.

    Объединяет кривую концентрации препарата в крови с моделью
    кислотности: эффективность ингибирования вычисляется в каждой точке
    из ``C_blood`` и ``molecule.ki``, затем pH-модель интегрируется
    с этой динамической эффективностью.
    """

    def __init__(self, molecule: Molecule, patient: Patient, dose_mg: float):
        """Сохраняет входные данные; симуляция не запускается.

        :param molecule: молекула с ``ki`` — константой ингибирования.
        :param patient: физиологический портрет пациента.
        :param dose_mg: разовая доза препарата, мг.
        """
        self.molecule = molecule
        self.patient = patient
        self.dose_mg = dose_mg

    def run(self, t_end_h: float = 24.0) -> pd.DataFrame:
        """Запускает PBPK-симуляцию и интегрирует pH-модель.

        Порядок работы:
        1. ``PBKPSolver`` даёт временной ряд ``C_blood`` по крови.
        2. ``efficacy = C_blood / (ki + C_blood)`` — насыщение, [0, 1).
        3. pH-модель интегрируется c этой эффективностью через
           ``np.interp`` в правой части ОДУ.
        4. Итоговый DataFrame объединяет оба ряда на общей сетке времени.

        :param t_end_h: горизонт симуляции, ч.
        :return: DataFrame с колонками ``time_h, C_blood, efficacy,
            acid_mmol, ph``.
        """
        pbpk_df = PBKPSolver(self.molecule, self.patient, self.dose_mg).run(
            t_end_h=t_end_h, t_step_min=_T_STEP_MIN
        )
        t = pbpk_df["time_h"].to_numpy()
        c_blood = pbpk_df["C_blood"].to_numpy()
        efficacy = c_blood / (self.molecule.ki + c_blood)

        acid_params = {
            "secretion_base": GASTRIC_SECRETION_BASE,
            "drug_efficacy": 0.0,  # replaced per step by the wrapper
            "evacuation_rate": GASTRIC_EVACUATION_RATE,
            "neutralization": GASTRIC_NEUTRALIZATION,
        }
        y0 = np.array([GASTRIC_INITIAL_ACID_MMOL])
        sol = solve_ivp(
            _gastric_rhs_with_dynamic_efficacy,
            t_span=(0.0, t_end_h),
            y0=y0,
            t_eval=t,
            args=(acid_params, t, efficacy),
            method="LSODA",
            rtol=1e-6,
            atol=1e-9,
        )
        if not sol.success:
            raise RuntimeError(f"Integration failed: {sol.message}")

        ph = np.array([calculate_ph(acid, V_GUT) for acid in sol.y[0]])
        return pd.DataFrame(
            {
                "time_h": t,
                "C_blood": c_blood,
                "efficacy": efficacy,
                "acid_mmol": sol.y[0],
                "ph": ph,
            }
        )
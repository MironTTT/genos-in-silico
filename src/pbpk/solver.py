"""Интеграция системы ОДУ PBPK-модели.

Модуль связывает молекулу, пациента и дозу с правой частью системы
ОДУ: строит словарь физиологических параметров (с учётом состояния
здоровья пациента) и запускает численное интегрирование через
``scipy.integrate.solve_ivp``.
"""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from config.constants import (
    CL_HEP_0,
    CL_OTHER_0,
    CL_REN_0,
    KA_BASE,
    Q_BRAIN,
    Q_KIDNEY,
    Q_LIVER,
    V_BLOOD,
    V_BRAIN,
    V_GUT,
    V_KIDNEY,
    V_LIVER,
)
from src.core.molecule import Molecule
from src.core.patient import Patient
from src.pbpk.bbb import permeability
from src.pbpk.equations import derivatives
from src.pbpk.partition import k_brain, k_kidney, k_liver


class PBKPSolver:
    """Класс-обёртка для интеграции PBPK-ОДУ.

    Принимает молекулу, пациента и дозу, формирует параметры модели
    и запускает ``solve_ivp``; результатом является DataFrame с рядами
    масс и концентраций по пяти компартментам.
    """

    def __init__(self, molecule: Molecule, patient: Patient, dose_mg: float):
        """Сохраняет входные данные; интеграция не запускается.

        :param molecule: молекула с физико-химическими дескрипторами.
        :param patient: физиологический портрет пациента.
        :param dose_mg: разовая доза препарата, мг.
        :raises ValueError: если ``dose_mg`` не положительна.
        """
        if dose_mg <= 0:
            raise ValueError(
                f"Поле 'dose_mg' должно быть положительным, получено {dose_mg}."
            )
        self.molecule = molecule
        self.patient = patient
        self.dose_mg = dose_mg

    def _build_params(self) -> dict:
        """Собирает словарь параметров для ``derivatives``.

        Кровотоки и клиренсы масштабируются на состояние здоровья
        пациента; эффективный церебральный кровоток учитывает
        проницаемость ГЭБ; коэффициенты распределения вычисляются
        по QSAR-зависимостям.

        :return: словарь параметров модели (см. ``equations.derivatives``).
        """
        patient = self.patient
        molecule = self.molecule
        return {
            # TODO: tie ka to solubility; fixed for now.
            "ka": KA_BASE,
            "q_liver": Q_LIVER * (0.3 + 0.7 * patient.liver_health),
            "q_kidney": Q_KIDNEY * (0.3 + 0.7 * patient.kidney_health),
            "q_brain_eff": Q_BRAIN
            * permeability(molecule.logp, molecule.tpsa, patient.bbb_integrity),
            "k_liver": k_liver(molecule.logp),
            "k_kidney": k_kidney(molecule.logp, molecule.tpsa),
            "k_brain": k_brain(molecule.logp, molecule.tpsa),
            "cl_hep": CL_HEP_0 * patient.liver_health,
            "cl_ren": CL_REN_0 * patient.kidney_health,
            "cl_other": CL_OTHER_0,
        }

    def run(self, t_end_h: float, t_step_min: float = 1.0) -> pd.DataFrame:
        """Интегрирует систему ОДУ и возвращает временные ряды.

        Начальное состояние — вся доза в ЖКТ. Интеграция выполняется
        методом LSODA; результат разворачивается в DataFrame с массами
        (A_*) и концентрациями (C_*, мг/л) по отсекам.

        :param t_end_h: длительность симуляции, ч.
        :param t_step_min: шаг вывода, мин.
        :return: DataFrame с колонками time_h, A_*, C_*.
        :raises ValueError: если ``t_end_h`` или ``t_step_min`` не положительны.
        :raises RuntimeError: если интегратор не завершился успешно.
        """
        if t_end_h <= 0:
            raise ValueError(
                f"Поле 't_end_h' должно быть положительным, получено {t_end_h}."
            )
        if t_step_min <= 0:
            raise ValueError(
                f"Поле 't_step_min' должно быть положительным, получено {t_step_min}."
            )

        params = self._build_params()
        y0 = np.array([self.dose_mg, 0.0, 0.0, 0.0, 0.0])
        dt_h = t_step_min / 60.0
        t_eval = np.arange(0.0, t_end_h + dt_h, dt_h)

        sol = solve_ivp(
            derivatives,
            t_span=(0.0, t_end_h),
            y0=y0,
            t_eval=t_eval,
            args=(params,),
            method="LSODA",
            rtol=1e-6,
            atol=1e-9,
        )
        if not sol.success:
            raise RuntimeError(f"Integration failed: {sol.message}")

        masses = {"A_gut": sol.y[0], "A_blood": sol.y[1], "A_liver": sol.y[2],
                  "A_kidney": sol.y[3], "A_brain": sol.y[4]}
        volumes = {"C_gut": V_GUT, "C_blood": V_BLOOD, "C_liver": V_LIVER,
                   "C_kidney": V_KIDNEY, "C_brain": V_BRAIN}

        df = pd.DataFrame({"time_h": sol.t})
        for name, mass in masses.items():
            df[name] = mass
        for name, vol in volumes.items():
            df[name] = masses[f"A_{name[2:]}"] / vol
        return df
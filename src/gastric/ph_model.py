"""Динамика кислотности желудка (модель гастрита).

Модуль описывает изменение количества кислоты в желудке (ммоль) под
действием трёх процессов — секреции, эвакуации и нейтрализации — и
пересчёт текущего количества кислоты в pH. Это первая часть модели
гастрита: желудок представлен одним хорошо перемешиваемым отсеком.
"""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from config.constants import (
    GASTRIC_EVACUATION_RATE,
    GASTRIC_INITIAL_ACID_MMOL,
    GASTRIC_NEUTRALIZATION,
    GASTRIC_PH_FLOOR,
    GASTRIC_SECRETION_BASE,
    V_GUT,
)


def gastric_acid_derivatives(t: float, y: np.ndarray, params: dict) -> np.ndarray:
    """Правая часть ОДУ: скорость изменения количества кислоты в желудке.

    Уравнение: dA/dt = secretion - evacuation_rate * A - neutralization,
    где secretion = secretion_base * (1 - drug_efficacy) — лекарство
    подавляет выработку кислоты пропорционально эффективности.

    :param t: время, ч (не используется — система автономна).
    :param y: массив [A], где A — количество кислоты, ммоль.
    :param params: словарь с ключами ``secretion_base``, ``drug_efficacy``,
        ``evacuation_rate`` и ``neutralization`` (см. ``constants``).
    :return: массив ``[dA_dt]``, форма (1,).
    """
    secretion = params["secretion_base"] * (1.0 - params["drug_efficacy"])
    d_acid = secretion - params["evacuation_rate"] * y[0] - params["neutralization"]
    return np.array([d_acid])


def calculate_ph(acid_amount_mmol: float, v_gut: float) -> float:
    """Пересчитывает количество кислоты в pH желудочного содержимого.

    Концентрация: C = acid_amount_mmol / 1000 / v_gut (моль/л). Логарифм
    защищён ``GASTRIC_PH_FLOOR``: если кислота полностью ушла, pH
    упирается в нейтральный предел ≈ 7.

    :param acid_amount_mmol: количество кислоты, ммоль.
    :param v_gut: объём отсека желудка, л.
    :return: pH желудочного содержимого.
    """
    concentration = acid_amount_mmol / 1000.0 / v_gut  # mol/L
    return -np.log10(max(concentration, GASTRIC_PH_FLOOR))


def run_gastric_simulation(drug_efficacy: float, t_end_h: float = 24.0) -> pd.DataFrame:
    """Интегрирует модель кислотности и возвращает временные ряды.

    Начальное условие — ``GASTRIC_INITIAL_ACID_MMOL`` кислоты в желудке.
    Интеграция выполняется методом LSODA (``solve_ivp``); результат —
    DataFrame с кислотной нагрузкой и рассчитанным pH.

    :param drug_efficacy: эффективность лекарства в [0.0, 1.0] — доля
        подавленной секреции (0.0 — без лечения, 1.0 — полная блокада).
    :param t_end_h: длительность симуляции, ч.
    :return: DataFrame с колонками ``time_h``, ``acid_mmol``, ``ph``.
    :raises ValueError: если ``drug_efficacy`` вне [0.0, 1.0] или
        ``t_end_h`` не положителен.
    """
    if not 0.0 <= drug_efficacy <= 1.0:
        raise ValueError(
            f"Поле 'drug_efficacy' должно быть в диапазоне [0.0, 1.0], "
            f"получено {drug_efficacy}."
        )
    if t_end_h <= 0:
        raise ValueError(
            f"Поле 't_end_h' должно быть положительным, получено {t_end_h}."
        )

    params = {
        "secretion_base": GASTRIC_SECRETION_BASE,
        "drug_efficacy": drug_efficacy,
        "evacuation_rate": GASTRIC_EVACUATION_RATE,
        "neutralization": GASTRIC_NEUTRALIZATION,
    }
    y0 = np.array([GASTRIC_INITIAL_ACID_MMOL])
    t_eval = np.arange(0.0, t_end_h + 0.25, 0.25)

    sol = solve_ivp(
        gastric_acid_derivatives,
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

    df = pd.DataFrame({"time_h": sol.t, "acid_mmol": sol.y[0]})
    df["ph"] = df["acid_mmol"].apply(lambda acid: calculate_ph(acid, V_GUT))
    return df
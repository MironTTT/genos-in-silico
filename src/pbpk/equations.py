"""Правая часть системы ОДУ пятикомпартментной PBPK-модели.

Модель перфузионно-ограниченная: ткань считается мгновенно
уравновешенной с отходящей кровью, поэтому концентрация в венозной
крови органа равна C_орган / K_орган. Вектор состояния — массы
препарата (мг) в отсеках ЖКТ, крови, печени, почек и мозга.
"""

import numpy as np

from config.constants import V_BLOOD, V_BRAIN, V_GUT, V_KIDNEY, V_LIVER


def derivatives(t: float, y: np.ndarray, params: dict) -> np.ndarray:
    """Правая часть системы ОДУ: скорости изменения масс по отсекам.

    Вектор состояния ``y = [A_gut, A_blood, A_liver, A_kidney, A_brain]`` —
    массы препарата в отсеках, мг. Концентрации вычисляются как
    ``C_i = A_i / V_i`` (мг/л), объёмы берутся из ``config.constants``.
    Система автономна, поэтому ``t`` в уравнения явно не входит.

    :param t: время, ч (не используется — система автономна).
    :param y: массив масс по пяти компартментам, форма (5,).
    :param params: словарь параметров модели (ka, q_liver, q_kidney,
        q_brain_eff, k_liver, k_kidney, k_brain, cl_hep, cl_ren, cl_other).
    :return: производные ``[dA_gut, dA_blood, dA_liver, dA_kidney, dA_brain]``,
        массив формы (5,).
    """
    c_blood = y[1] / V_BLOOD
    c_liver = y[2] / V_LIVER
    c_kidney = y[3] / V_KIDNEY
    c_brain = y[4] / V_BRAIN

    # Gut: first-order absorption into the blood.
    dy = np.zeros(5)
    dy[0] = -params["ka"] * y[0]

    # Blood: absorption in, venous return from the organs, other clearance.
    dy[1] = (
        params["ka"] * y[0]
        + params["q_liver"] * (c_liver / params["k_liver"] - c_blood)
        + params["q_kidney"] * (c_kidney / params["k_kidney"] - c_blood)
        + params["q_brain_eff"] * (c_brain / params["k_brain"] - c_blood)
        - params["cl_other"] * c_blood
    )

    # Liver: uptake from blood, elimination by hepatic clearance.
    dy[2] = (
        params["q_liver"] * (c_blood - c_liver / params["k_liver"])
        - params["cl_hep"] * (c_liver / params["k_liver"])
    )

    # Kidney: uptake from blood, elimination by renal clearance.
    dy[3] = (
        params["q_kidney"] * (c_blood - c_kidney / params["k_kidney"])
        - params["cl_ren"] * (c_kidney / params["k_kidney"])
    )

    # Brain: exchange with blood only, no elimination.
    dy[4] = params["q_brain_eff"] * (c_blood - c_brain / params["k_brain"])

    return dy


# Mass-balance sanity check: organ exchange terms cancel out exactly
# (whatever leaves the blood enters an organ and vice versa), so the sum
# of all dA_i/dt equals the net elimination:
#   - cl_hep * (C_liver / k_liver) - cl_ren * (C_kidney / k_kidney) - cl_other * C_blood
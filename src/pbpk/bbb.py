"""Оценка проницаемости гематоэнцефалического барьера (ГЭБ).

Модуль вычисляет безразмерный коэффициент ``P_BBB ∈ [0, 1]`` — долю
проникновения молекулы из плазмы в мозг через ГЭБ. В PBPK-модели
(Шаг 6) ``P_BBB`` умножается на церебральный кровоток ``Q_BRAIN``,
т.е. описывает кинетическую доступность мозга, а не равновесное
распределение (``k_brain`` из Шага 4).

Модуль не зависит от ``src.core``: функция принимает скалярные
аргументы, а не объект ``Molecule``.
"""

from math import exp

# Lipophilicity threshold: penetration grows above it, falls below it.
_LOGP_THRESHOLD = 1.5

# Reference TPSA (Å²) at which polarity is neutral.
_TPSA_REFERENCE = 60.0

# Polarity penalty per TPSA unit above the reference.
_TPSA_PENALTY = 0.05


def permeability(logp: float, tpsa: float, bbb_integrity: float = 1.0) -> float:
    """Доля проникновения молекулы через ГЭБ из плазмы в мозг.

    Модельная сигмоида Clark-типа с аддитивным score:

        score = (logp - 1.5) - 0.05 * (tpsa - 60)
        P_BBB = 1 / (1 + exp(-score))
        P_BBB = max(0.0, min(1.0, P_BBB)) * bbb_integrity

    Аддитивная форма выбрана намеренно: штрафы за низкую липофильность
    и высокую полярность складываются независимо, а не перемножаются.
    Мультипликативная форма ``(2.0 - 0.02*tpsa) * (logp - 1.5)`` ломается
    при двух одновременных штрафах: два отрицательных множителя дают
    положительное произведение, и полярная молекула с низким LogP
    ошибочно «проходит» ГЭБ.

    Логика: чем выше LogP (до разумных пределов) — тем лучше
    проникновение; чем выше TPSA сверх 60 Å² — тем сильнее штраф
    (полярные молекулы плохо проходят ГЭБ). ``bbb_integrity ∈ [0, 1]`` —
    параметр состояния пациента из ``Patient``: 1.0 — здоровая ГЭБ,
    0.0 — полностью нарушенная (клинически: нейровоспаление, травма,
    поздняя стадия нейродегенерации).

    Санити-чек: ``permeability(2.0, 60.0, 1.0)`` ≈ 0.622 (score = 0.5,
    нейтральная полярность).

    :param logp: коэффициент распределения октанол/вода (cLogP).
    :param tpsa: топологическая полярная площадь поверхности, Å².
    :param bbb_integrity: целостность ГЭБ пациента в диапазоне [0.0, 1.0].
    :return: коэффициент P_BBB в диапазоне [0.0, 1.0].
    :raises ValueError: если ``bbb_integrity`` вне диапазона [0.0, 1.0].
    """
    if not 0.0 <= bbb_integrity <= 1.0:
        raise ValueError(
            "Поле 'bbb_integrity' должно быть в диапазоне [0.0, 1.0], "
            f"получено {bbb_integrity}."
        )

    score = (logp - _LOGP_THRESHOLD) - _TPSA_PENALTY * (tpsa - _TPSA_REFERENCE)
    sig = 1.0 / (1.0 + exp(-score))
    # Clamp guards against float noise; bbb_integrity <= 1 keeps the result in [0, 1].
    return max(0.0, min(1.0, sig)) * bbb_integrity
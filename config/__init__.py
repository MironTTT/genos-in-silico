"""Пакет ``config`` — конфигурация и физические константы PBPK-модели.

Реэкспортирует все константы эталонного взрослого из ``config.constants``,
чтобы потребители могли импортировать их одной строкой:
``from config import V_BLOOD, Q_LIVER``.
"""

from config.constants import *  # noqa: F401,F403
from config.constants import __all__  # noqa: F401
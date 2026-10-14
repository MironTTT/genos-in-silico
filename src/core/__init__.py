"""Пакет ``src.core`` — ядро PBPK-симулятора.

Содержит доменные модели и алгоритмы симуляции, независимые
от генератора молекул и пользовательского интерфейса.
"""

from src.core.molecule import Molecule
from src.core.patient import Patient

__all__ = ["Molecule", "Patient"]
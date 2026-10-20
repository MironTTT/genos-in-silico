"""Пакет ``src.pbpk`` — PBPK-моделирование.

Содержит QSAR-оценки коэффициентов распределения ткань/плазма
и последующие компоненты физиологически обоснованной
фармакокинетической модели.
"""

from src.pbpk.bbb import permeability
from src.pbpk.partition import k_brain, k_kidney, k_liver

__all__ = ["k_liver", "k_kidney", "k_brain", "permeability"]
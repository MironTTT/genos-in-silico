"""Пакет ``src.gastric`` — модели заболеваний желудка (гастрит).

Первый модуль — ``ph_model``: динамика кислотности желудка (pH) как
первая часть модели гастрита. Пакет не зависит от ``src.pbpk`` и
``src.core``, использует только ``config.constants``.
"""

from src.gastric.ph_model import (
    calculate_ph,
    gastric_acid_derivatives,
    run_gastric_simulation,
)

__all__ = ["gastric_acid_derivatives", "calculate_ph", "run_gastric_simulation"]
"""Smoke-тест Шага 3: проверка загрузки и согласованности констант PBPK-модели.

Запуск из корня проекта: ``python -m scripts.smoke_test_step3``.
"""

from config.constants import *

print(f"V_BLOOD={V_BLOOD}, Q_LIVER={Q_LIVER}, CL_HEP_0={CL_HEP_0}")

# Сумма органных кровотоков не должна превышать сердечный выброс,
# пересчитанный в л/ч (согласовано с assert в config/constants.py).
assert Q_LIVER + Q_KIDNEY + Q_BRAIN <= CARDIAC_OUTPUT_BASE * 60
print("[OK] Сумма кровотоков не превышает сердечный выброс")
print("[OK] Все константы импортируются")
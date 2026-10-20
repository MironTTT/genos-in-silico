"""Модель пациента ядра PBPK-симулятора.

Модуль определяет ``Patient`` — физиологический портрет субъекта
симуляции: антропометрия, возраст и параметры состояния органов,
влияющие на фармакокинетику.
"""

from dataclasses import dataclass, fields

_HEALTH_FIELDS = ("liver_health", "kidney_health", "bbb_integrity")

_HEALTH_MIN = 0.0
_HEALTH_MAX = 1.0
_AGE_MIN = 0
_AGE_MAX = 120
_WEIGHT_MIN = 1.0
_WEIGHT_MAX = 300.0
_CARDIAC_MIN = 1.0
_CARDIAC_MAX = 15.0


def _check_number(value: object, field_name: str) -> float:
    """Проверяет, что значение является конечным числом, и приводит к float.

    :param value: проверяемое значение.
    :param field_name: имя поля для сообщения об ошибке.
    :return: значение как ``float``.
    :raises ValueError: если значение не является числом (включая bool).
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(
            f"Поле {field_name!r} должно быть числом, получено {type(value).__name__}."
        )
    return float(value)


@dataclass
class Patient:
    """Физиологическая модель пациента для PBPK-симуляции.

    Коэффициенты ``liver_health``, ``kidney_health`` и ``bbb_integrity``
    (целостность гематоэнцефалического барьера) нормированы на
    диапазон [0.0, 1.0], где 1.0 — полностью здоровый орган.
    """

    age: int
    """Возраст, годы (целое в диапазоне [0, 120])."""

    weight_kg: float
    """Масса тела, кг (в диапазоне [1, 300])."""

    liver_health: float = 1.0
    """Функциональное состояние печени в диапазоне [0.0, 1.0]."""

    kidney_health: float = 1.0
    """Функциональное состояние почек в диапазоне [0.0, 1.0]."""

    bbb_integrity: float = 1.0
    """Целостность гематоэнцефалического барьера в диапазоне [0.0, 1.0]."""

    cardiac_output: float = 6.5
    """Сердечный выброс, л/мин (в диапазоне [1.0, 15.0])."""

    def __post_init__(self) -> None:
        """Валидирует поля пациента после инициализации.

        Проверяются типы и допустимые диапазоны значений; health-поля
        обязаны лежать в [0.0, 1.0].

        :raises ValueError: если хотя бы одно поле вне допустимого диапазона.
        """
        if isinstance(self.age, bool) or not isinstance(self.age, int):
            raise ValueError(
                f"Поле 'age' должно быть целым числом в диапазоне "
                f"[{_AGE_MIN}, {_AGE_MAX}], получено {self.age!r} "
                f"(тип {type(self.age).__name__})."
            )
        if not _AGE_MIN <= self.age <= _AGE_MAX:
            raise ValueError(
                f"Поле 'age' должно быть в диапазоне "
                f"[{_AGE_MIN}, {_AGE_MAX}], получено {self.age}."
            )

        self.weight_kg = _check_number(self.weight_kg, "weight_kg")
        if not _WEIGHT_MIN <= self.weight_kg <= _WEIGHT_MAX:
            raise ValueError(
                f"Поле 'weight_kg' должно быть в диапазоне "
                f"[{_WEIGHT_MIN}, {_WEIGHT_MAX}], получено {self.weight_kg}."
            )

        for field_name in _HEALTH_FIELDS:
            value = _check_number(getattr(self, field_name), field_name)
            if not _HEALTH_MIN <= value <= _HEALTH_MAX:
                raise ValueError(
                    f"Поле {field_name!r} должно быть в диапазоне "
                    f"[{_HEALTH_MIN}, {_HEALTH_MAX}], получено {value}."
                )
            setattr(self, field_name, value)

        self.cardiac_output = _check_number(self.cardiac_output, "cardiac_output")
        if not _CARDIAC_MIN <= self.cardiac_output <= _CARDIAC_MAX:
            raise ValueError(
                f"Поле 'cardiac_output' должно быть в диапазоне "
                f"[{_CARDIAC_MIN}, {_CARDIAC_MAX}], получено {self.cardiac_output}."
            )

    def to_dict(self) -> dict[str, float | int]:
        """Возвращает все поля пациента как словарь.

        :return: словарь с именами полей в качестве ключей.
        """
        return {field.name: getattr(self, field.name) for field in fields(self)}

    def copy_with(self, **kwargs) -> "Patient":
        """Создаёт копию пациента с изменёнными полями.

        Базовые значения берутся из ``to_dict``, затем переопределяются
        переданными аргументами; новая копия проходит полную валидацию.

        :param kwargs: поля пациента и их новые значения.
        :return: новый экземпляр ``Patient``.
        :raises ValueError: если новое значение поля вне допустимого диапазона.
        """
        values = self.to_dict()
        values.update(kwargs)
        return Patient(**values)
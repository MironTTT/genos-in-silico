"""Молекулярная модель ядра PBPK-симулятора.

Модуль определяет неизменяемое представление молекулы ``Molecule``
и фабричный метод ``from_smiles``, вычисляющий физико-химические
дескрипторы через RDKit.
"""

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

# Cap on SMILES length so pathological inputs can't hang the RDKit parser.
_MAX_SMILES_LENGTH = 10_000

# ESOL (Delaney, 2004) solubility model constants.
_ESOL_INTERCEPT = 0.16
_ESOL_LOGP_COEF = -0.63
_ESOL_MW_COEF = -0.0062
_ESOL_RB_COEF = 0.066
_ESOL_AROM_COEF = -0.74


def _esol_log_s(mol: Chem.Mol, logp: float, mw: float, rotatable_bonds: int) -> float:
    """Оценка логарифма водорастворимости (logS, моль/л) по модели ESOL.

    Модель Delaney использует logP, молекулярную массу, число
    вращающихся связей и долю ароматических атомов.

    :param mol: молекула RDKit, для которой считается доля ароматических атомов.
    :param logp: коэффициент распределения октанол/вода (cLogP).
    :param mw: молекулярная масса, Да.
    :param rotatable_bonds: число вращающихся связей.
    :return: прогноз logS в моль/л.
    """
    heavy_atoms = mol.GetNumHeavyAtoms()
    aromatic_atoms = sum(
        1 for atom in mol.GetAtoms() if atom.GetIsAromatic()
    )
    aromatic_proportion = aromatic_atoms / heavy_atoms if heavy_atoms else 0.0
    return (
        _ESOL_INTERCEPT
        + _ESOL_LOGP_COEF * logp
        + _ESOL_MW_COEF * mw
        + _ESOL_RB_COEF * rotatable_bonds
        + _ESOL_AROM_COEF * aromatic_proportion
    )


@dataclass(frozen=True)
class Molecule:
    """Неизменяемое описание молекулы с физико-химическими дескрипторами.

    Все поля вычисляются из SMILES при создании через ``from_smiles``
    и не могут быть изменены после инициализации.
    """

    smiles: str
    """Каноническая SMILES-запись молекулы."""

    mw: float
    """Молекулярная масса, Да."""

    logp: float
    """Коэффициент распределения октанол/вода (cLogP)."""

    tpsa: float
    """Топологическая полярная площадь поверхности, Å²."""

    hbd: int
    """Число доноров водородной связи."""

    hba: int
    """Число акцепторов водородной связи."""

    rotatable_bonds: int
    """Число вращающихся связей."""

    solubility_logs: float
    """Логарифм водорастворимости (logS, моль/л) по модели ESOL."""

    @classmethod
    def from_smiles(cls, smiles: str) -> "Molecule":
        """Создаёт ``Molecule`` из строки SMILES.

        Дескрипторы вычисляются через RDKit: ``MolFromSmiles`` для
        парсинга, ``Descriptors.MolWt`` — молекулярная масса,
        ``Crippen.MolLogP`` — logP, ``Descriptors.TPSA`` — полярная
        площадь поверхности; доноры/акцепторы водородной связи и число
        вращающихся связей берутся из ``Descriptors``, растворимость —
        по модели ESOL.

        :param smiles: SMILES-строка молекулы.
        :return: экземпляр ``Molecule`` с заполненными дескрипторами.
        :raises ValueError: если SMILES пуст, превышает допустимую длину
            или не может быть разобран RDKit.
        """
        if not isinstance(smiles, str) or not smiles.strip():
            raise ValueError("SMILES должен быть непустой строкой.")

        if len(smiles) > _MAX_SMILES_LENGTH:
            raise ValueError(
                f"SMILES слишком длинный ({len(smiles)} символов, "
                f"лимит {_MAX_SMILES_LENGTH})."
            )

        try:
            mol = Chem.MolFromSmiles(smiles)
        except Exception as exc:  # noqa: BLE001 - RDKit raises various exception types
            raise ValueError(f"Не удалось разобрать SMILES {smiles!r}: {exc}") from exc

        if mol is None:
            raise ValueError(f"Невалидный SMILES: {smiles!r}.")

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)

        return cls(
            smiles=Chem.MolToSmiles(mol),
            mw=mw,
            logp=logp,
            tpsa=Descriptors.TPSA(mol),
            hbd=Descriptors.NumHDonors(mol),
            hba=Descriptors.NumHAcceptors(mol),
            rotatable_bonds=rotatable_bonds,
            solubility_logs=_esol_log_s(
                mol, logp=logp, mw=mw, rotatable_bonds=rotatable_bonds
            ),
        )
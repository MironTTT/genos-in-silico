"""Проверка: значения из agent-кода совпадают с «сырым» RDKit."""
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, Crippen

print(f"RDKit version: {rdBase.rdkitVersion}")
print(f"Python: {__import__('sys').version.split()[0]}")

def probe(smiles: str, label: str) -> None:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        print(f"[{label}] ❌ SMILES не парсится")
        return
    print(f"\n[{label}]")
    print(f"  Input:        {smiles}")
    print(f"  Canonical:    {Chem.MolToSmiles(m)}")
    print(f"  MW:           {Descriptors.MolWt(m):.4f}")
    print(f"  LogP:         {Crippen.MolLogP(m):.4f}")
    print(f"  TPSA:         {Descriptors.TPSA(m):.4f}")
    print(f"  HBD:          {Descriptors.NumHDonors(m)}")
    print(f"  HBA:          {Descriptors.NumHAcceptors(m)}")
    print(f"  RotB:         {Descriptors.NumRotatableBonds(m)}")

probe("CC(=O)Oc1ccccc1C(=O)O",              "Aspirin")
probe("CN1C=NC2=C1C(=O)N(C(=O)N2C)C",       "Caffeine")
probe("CC(=O)O",                            "Acetic acid (sanity)")
probe("c1ccccc1",                           "Benzene (sanity)")
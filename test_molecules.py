"""Run:  python test_molecules.py   (needs rdkit installed)"""
from rdkit import Chem
from molecules import MOLECULE_LIBRARY

seen, bad, total = set(), [], 0
for cat, items in MOLECULE_LIBRARY.items():
    for name, smi in items.items():
        total += 1
        if name in seen:
            bad.append((name, "DUPLICATE NAME"))
        seen.add(name)
        if Chem.MolFromSmiles(smi) is None:
            bad.append((name, f"INVALID SMILES: {smi}"))

print(f"Checked {total} molecules in {len(MOLECULE_LIBRARY)} categories")
for name, why in bad:
    print("  ❌", name, "->", why)
print("✅ All good" if not bad else f"{len(bad)} problem(s) found")

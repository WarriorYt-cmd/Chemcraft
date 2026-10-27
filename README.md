# 🧪 Cheminformatics Virtual Laboratory

## Aim
To develop a foundational understanding of Cheminformatics by using
computational tools to represent, visualize and analyze chemical structures.

## Learning Objectives
Students will be able to
- Convert molecular representations using SMILES.
- Visualize molecular structures in 2D and 3D.
- Calculate molecular descriptors.
- Understand structure–property relationships.

## Molecule Library
100+ molecules in 12 categories (drugs, amino acids, sugars, vitamins,
natural products, industrial chemicals, ...). Add your own in `molecules.py`.

## Molecular Descriptors
Molecular Weight, LogP, TPSA, H-bond donors/acceptors, rotatable bonds,
ring count, σ/π bonds, Fsp3, molar refractivity and more.

## Run locally
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Check the molecule library
```bash
python test_molecules.py
```

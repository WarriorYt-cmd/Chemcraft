"""
Molecule library for the Cheminformatics Virtual Laboratory.

Structure:  { "Category name": { "Molecule name": "SMILES", ... }, ... }

To add a molecule, just add a new  "Name": "SMILES"  line to any category
(or create a new category). Molecule names must be unique across the library.
The app validates every SMILES at start-up and silently skips invalid ones.
"""

MOLECULE_LIBRARY = {

    "Simple Molecules & Gases": {
        "Water": "O",
        "Ammonia": "N",
        "Carbon Dioxide": "O=C=O",
        "Carbon Monoxide": "[C-]#[O+]",
        "Hydrogen Cyanide": "C#N",
        "Hydrogen Peroxide": "OO",
        "Ozone": "[O-][O+]=O",
        "Sulfuric Acid": "OS(=O)(=O)O",
        "Nitric Acid": "O[N+](=O)[O-]",
        "Sodium Chloride": "[Na+].[Cl-]",
    },

    "Alkanes, Alkenes & Alkynes": {
        "Methane": "C",
        "Ethane": "CC",
        "Propane": "CCC",
        "Butane": "CCCC",
        "Hexane": "CCCCCC",
        "Cyclohexane": "C1CCCCC1",
        "Ethene": "C=C",
        "Isoprene": "CC(=C)C=C",
        "Acetylene": "C#C",
    },

    "Alcohols, Ethers & Halides": {
        "Methanol": "CO",
        "Ethanol": "CCO",
        "Isopropanol": "CC(C)O",
        "Ethylene Glycol": "OCCO",
        "Glycerol": "OCC(O)CO",
        "Diethyl Ether": "CCOCC",
        "Tetrahydrofuran": "C1CCOC1",
        "Dichloromethane": "ClCCl",
        "Chloroform": "ClC(Cl)Cl",
        "Menthol": "CC(C)C1CCC(C)CC1O",
    },

    "Aldehydes, Ketones, Acids & Amides": {
        "Formaldehyde": "C=O",
        "Acetaldehyde": "CC=O",
        "Acetone": "CC(C)=O",
        "Formic Acid": "OC=O",
        "Acetic Acid": "CC(=O)O",
        "Lactic Acid": "CC(O)C(=O)O",
        "Citric Acid": "OC(=O)CC(O)(CC(=O)O)C(=O)O",
        "Ethyl Acetate": "CCOC(C)=O",
        "Acetamide": "CC(N)=O",
        "Urea": "NC(N)=O",
        "Acetonitrile": "CC#N",
        "DMSO": "CS(C)=O",
    },

    "Aromatic & Heterocyclic Compounds": {
        "Benzene": "c1ccccc1",
        "Toluene": "Cc1ccccc1",
        "Phenol": "Oc1ccccc1",
        "Aniline": "Nc1ccccc1",
        "Benzoic Acid": "OC(=O)c1ccccc1",
        "Nitrobenzene": "[O-][N+](=O)c1ccccc1",
        "Styrene": "C=Cc1ccccc1",
        "Naphthalene": "c1ccc2ccccc2c1",
        "Anthracene": "c1ccc2cc3ccccc3cc2c1",
        "Pyridine": "c1ccncc1",
        "Pyrrole": "c1cc[nH]c1",
        "Furan": "c1ccoc1",
        "Thiophene": "c1ccsc1",
        "Imidazole": "c1c[nH]cn1",
        "Indole": "c1ccc2[nH]ccc2c1",
    },

    "Amino Acids": {
        "Glycine": "NCC(=O)O",
        "Alanine": "CC(N)C(=O)O",
        "Serine": "NC(CO)C(=O)O",
        "Cysteine": "NC(CS)C(=O)O",
        "Lysine": "NCCCCC(N)C(=O)O",
        "Glutamic Acid": "NC(CCC(=O)O)C(=O)O",
        "Histidine": "NC(Cc1cnc[nH]1)C(=O)O",
        "Phenylalanine": "NC(Cc1ccccc1)C(=O)O",
        "Tyrosine": "NC(Cc1ccc(O)cc1)C(=O)O",
        "Tryptophan": "NC(Cc1c[nH]c2ccccc12)C(=O)O",
    },

    "Sugars, Lipids & Nucleobases": {
        "Glucose": "OCC1OC(O)C(O)C(O)C1O",
        "Fructose": "OCC1(O)OCC(O)C(O)C1O",
        "Sucrose": "OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O",
        "Palmitic Acid": "CCCCCCCCCCCCCCCC(=O)O",
        "Oleic Acid": "CCCCCCCC/C=C\\CCCCCCCC(=O)O",
        "Cholesterol": "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C",
        "Adenine": "Nc1ncnc2[nH]cnc12",
        "Cytosine": "Nc1cc[nH]c(=O)n1",
        "Thymine": "Cc1c[nH]c(=O)[nH]c1=O",
        "Uracil": "O=c1cc[nH]c(=O)[nH]1",
    },

    "Hormones, Neurotransmitters & Vitamins": {
        "Dopamine": "NCCc1ccc(O)c(O)c1",
        "Serotonin": "NCCc1c[nH]c2ccc(O)cc12",
        "Adrenaline": "CNCC(O)c1ccc(O)c(O)c1",
        "Melatonin": "COc1ccc2[nH]cc(CCNC(C)=O)c2c1",
        "Testosterone": "CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O",
        "Estradiol": "CC12CCC3c4ccc(O)cc4CCC3C2CCC1O",
        "Vitamin C": "OCC(O)C1OC(=O)C(O)=C1O",
        "Vitamin B3 (Niacin)": "OC(=O)c1cccnc1",
        "Vitamin A (Retinol)": "CC1=C(C(C)(C)CCC1)C=CC(C)=CC=CC(C)=CCO",
    },

    "Common Drugs": {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
        "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "Naproxen": "COc1ccc2cc(C(C)C(=O)O)ccc2c1",
        "Metformin": "CN(C)C(=N)NC(N)=N",
        "Amoxicillin": "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O",
        "Penicillin G": "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",
        "Ciprofloxacin": "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
        "Diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
        "Fluoxetine": "CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1",
        "Omeprazole": "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C",
        "Lidocaine": "CCN(CC)CC(=O)Nc1c(C)cccc1C",
        "Propranolol": "CC(C)NCC(O)COc1cccc2ccccc12",
        "Salbutamol": "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",
        "Cetirizine": "OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1",
        "Warfarin": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
        "Tamoxifen": "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
        "Sildenafil": "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1",
        "Atorvastatin": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
        "Morphine": "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O",
    },

    "Stimulants & Food Chemistry": {
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Theobromine": "Cn1cnc2c1c(=O)[nH]c(=O)n2C",
        "Nicotine": "CN1CCCC1c1cccnc1",
        "Capsaicin": "COc1cc(CNC(=O)CCCCC=CC(C)C)ccc1O",
        "Vanillin": "COc1cc(C=O)ccc1O",
        "Saccharin": "O=C1NS(=O)(=O)c2ccccc12",
        "Aspartame": "COC(=O)C(Cc1ccccc1)NC(=O)C(N)CC(=O)O",
    },

    "Natural Products": {
        "Limonene": "CC1=CCC(CC1)C(C)=C",
        "Camphor": "CC1(C)C2CCC1(C)C(=O)C2",
        "Eugenol": "C=CCc1ccc(O)c(OC)c1",
        "Resveratrol": "Oc1ccc(C=Cc2cc(O)cc(O)c2)cc1",
        "Curcumin": "COc1cc(C=CC(=O)CC(=O)C=Cc2ccc(O)c(OC)c2)ccc1O",
        "Quercetin": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
    },

    "Industrial & Environmental": {
        "Ethylene Oxide": "C1CO1",
        "Vinyl Chloride": "ClC=C",
        "Acrylonitrile": "C=CC#N",
        "Caprolactam": "O=C1CCCCCN1",
        "Bisphenol A": "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1",
        "MTBE": "COC(C)(C)C",
        "Freon-12 (CFC-12)": "FC(F)(Cl)Cl",
        "DDT": "ClC(Cl)(Cl)C(c1ccc(Cl)cc1)c1ccc(Cl)cc1",
        "TNT": "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",
    },
}


def flat_library():
    """Return a single {name: smiles} dictionary of every molecule."""
    flat = {}
    for category in MOLECULE_LIBRARY.values():
        flat.update(category)
    return flat

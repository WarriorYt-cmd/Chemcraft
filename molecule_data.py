"""
Experimental physical properties for molecules in molecules.py and automated
structure-property intermolecular forces analysis.

IMPORTANT
---------
RDKit can calculate theoretical descriptors like molecular weight, LogP, TPSA,
and H-bond donors/acceptors, but experimental melting point (MP), boiling point (BP),
density, and room-temperature physical state are literature-derived values
(CRC Handbook, Merck Index, and PubChem at 1 atm unless noted).
"""

import re

D = "—"   # "no data"

EXPERIMENTAL = {
    # ---- Simple molecules & gases ----
    "Water":               ("0", "100", "1.00", "Liquid", "Miscible"),
    "Ammonia":             ("-77.7", "-33.3", D, "Gas", "Very soluble"),
    "Carbon Dioxide":      ("-56.6 (5.2 atm)", "-78.5 (sublimes)", D, "Gas", "Slightly soluble"),
    "Carbon Monoxide":     ("-205", "-191.5", D, "Gas", "Slightly soluble"),
    "Hydrogen Cyanide":    ("-13.4", "25.6", "0.69", "Liquid", "Miscible"),
    "Hydrogen Peroxide":   ("-0.4", "150 (dec.)", "1.45", "Liquid", "Miscible"),
    "Ozone":               ("-192", "-112", D, "Gas", "Slightly soluble"),
    "Sulfuric Acid":       ("10.3", "337 (dec.)", "1.83", "Liquid", "Miscible"),
    "Nitric Acid":         ("-42", "83", "1.51", "Liquid", "Miscible"),
    "Sodium Chloride":     ("801", "1413", "2.16", "Solid", "Soluble"),

    # ---- Alkanes, alkenes, alkynes ----
    "Methane":             ("-182.5", "-161.5", D, "Gas", "Very slightly soluble"),
    "Ethane":              ("-182.8", "-88.6", D, "Gas", "Very slightly soluble"),
    "Propane":             ("-187.7", "-42.1", D, "Gas", "Very slightly soluble"),
    "Butane":              ("-138.4", "-0.5", D, "Gas", "Very slightly soluble"),
    "Hexane":              ("-95", "69", "0.66", "Liquid", "Insoluble"),
    "Cyclohexane":         ("6.5", "80.7", "0.78", "Liquid", "Insoluble"),
    "Ethene":              ("-169.2", "-103.7", D, "Gas", "Slightly soluble"),
    "Isoprene":            ("-146", "34", "0.68", "Liquid", "Insoluble"),
    "Acetylene":           ("-80.8 (1.2 atm)", "-84 (sublimes)", D, "Gas", "Slightly soluble"),

    # ---- Alcohols, ethers, halides ----
    "Methanol":            ("-97.6", "64.7", "0.79", "Liquid", "Miscible"),
    "Ethanol":             ("-114.1", "78.4", "0.79", "Liquid", "Miscible"),
    "Isopropanol":         ("-89", "82.6", "0.79", "Liquid", "Miscible"),
    "Ethylene Glycol":     ("-12.9", "197.3", "1.11", "Liquid", "Miscible"),
    "Glycerol":            ("17.8", "290 (dec.)", "1.26", "Liquid", "Miscible"),
    "Diethyl Ether":       ("-116.3", "34.6", "0.71", "Liquid", "Slightly soluble"),
    "Tetrahydrofuran":     ("-108.4", "66", "0.89", "Liquid", "Miscible"),
    "Dichloromethane":     ("-96.7", "39.6", "1.33", "Liquid", "Slightly soluble"),
    "Chloroform":          ("-63.5", "61.2", "1.49", "Liquid", "Slightly soluble"),
    "Menthol":             ("42", "212", "0.89", "Solid", "Slightly soluble"),

    # ---- Aldehydes, ketones, acids, amides ----
    "Formaldehyde":        ("-92", "-19", D, "Gas", "Very soluble"),
    "Acetaldehyde":        ("-123.4", "20.2", "0.78", "Liquid", "Miscible"),
    "Acetone":             ("-94.7", "56", "0.78", "Liquid", "Miscible"),
    "Formic Acid":         ("8.4", "100.8", "1.22", "Liquid", "Miscible"),
    "Acetic Acid":         ("16.6", "118.1", "1.05", "Liquid", "Miscible"),
    "Lactic Acid":         ("17", "122 (15 mmHg)", "1.21", "Liquid", "Miscible"),
    "Citric Acid":         ("153", "175 (dec.)", "1.66", "Solid", "Very soluble"),
    "Ethyl Acetate":       ("-83.6", "77.1", "0.90", "Liquid", "Soluble"),
    "Acetamide":           ("80", "221", "1.16", "Solid", "Very soluble"),
    "Urea":                ("133", "dec.", "1.32", "Solid", "Very soluble"),
    "Acetonitrile":        ("-45.7", "81.6", "0.79", "Liquid", "Miscible"),
    "DMSO":                ("18.5", "189", "1.10", "Liquid", "Miscible"),

    # ---- Aromatic & heterocyclic ----
    "Benzene":             ("5.5", "80.1", "0.88", "Liquid", "Slightly soluble"),
    "Toluene":             ("-95", "110.6", "0.87", "Liquid", "Insoluble"),
    "Phenol":              ("40.5", "181.7", "1.07", "Solid", "Soluble"),
    "Aniline":             ("-6", "184", "1.02", "Liquid", "Slightly soluble"),
    "Benzoic Acid":        ("122.4", "249", "1.27", "Solid", "Slightly soluble"),
    "Nitrobenzene":        ("5.7", "211", "1.20", "Liquid", "Slightly soluble"),
    "Styrene":             ("-31", "145", "0.91", "Liquid", "Insoluble"),
    "Naphthalene":         ("80.3", "218", "1.14", "Solid", "Insoluble"),
    "Anthracene":          ("216", "340", "1.28", "Solid", "Insoluble"),
    "Pyridine":            ("-41.6", "115.2", "0.98", "Liquid", "Miscible"),
    "Pyrrole":             ("-23", "130", "0.97", "Liquid", "Slightly soluble"),
    "Furan":               ("-85.6", "31.4", "0.94", "Liquid", "Insoluble"),
    "Thiophene":           ("-38", "84", "1.06", "Liquid", "Insoluble"),
    "Imidazole":           ("90", "256", "1.23", "Solid", "Soluble"),
    "Indole":              ("52.5", "253", "1.22", "Solid", "Slightly soluble"),

    # ---- Amino acids (L-form, decompose on heating) ----
    "Glycine":             ("233 (dec.)", D, "1.61", "Solid", "Soluble"),
    "Alanine":             ("258 (dec.)", D, "1.42", "Solid", "Soluble"),
    "Serine":              ("228 (dec.)", D, "1.60", "Solid", "Soluble"),
    "Cysteine":            ("240 (dec.)", D, D, "Solid", "Soluble"),
    "Lysine":              ("224 (dec.)", D, D, "Solid", "Very soluble"),
    "Glutamic Acid":       ("199 (dec.)", D, "1.54", "Solid", "Slightly soluble"),
    "Histidine":           ("287 (dec.)", D, D, "Solid", "Soluble"),
    "Phenylalanine":       ("283 (dec.)", D, "1.29", "Solid", "Soluble"),
    "Tyrosine":            ("343 (dec.)", D, D, "Solid", "Slightly soluble"),
    "Tryptophan":          ("289 (dec.)", D, D, "Solid", "Soluble"),

    # ---- Sugars, lipids, nucleobases ----
    "Glucose":             ("146", "dec.", "1.54", "Solid", "Very soluble"),
    "Fructose":            ("103 (dec.)", D, "1.60", "Solid", "Very soluble"),
    "Sucrose":             ("186 (dec.)", D, "1.59", "Solid", "Very soluble"),
    "Palmitic Acid":       ("63", "351", "0.85", "Solid", "Insoluble"),
    "Oleic Acid":          ("13", "360 (dec.)", "0.90", "Liquid", "Insoluble"),
    "Cholesterol":         ("148.5", "360 (dec.)", "1.07", "Solid", "Insoluble"),
    "Adenine":             ("360 (dec.)", D, "1.60", "Solid", "Slightly soluble"),
    "Cytosine":            ("320 (dec.)", D, D, "Solid", "Slightly soluble"),
    "Thymine":             ("316", D, D, "Solid", "Slightly soluble"),
    "Uracil":              ("335 (dec.)", D, "1.32", "Solid", "Slightly soluble"),

    # ---- Hormones, neurotransmitters, vitamins ----
    "Dopamine":            ("128", D, D, "Solid", "Soluble"),
    "Serotonin":           ("167.5", "dec.", "1.25", "Solid", "Soluble"),
    "Adrenaline":          ("212 (dec.)", D, D, "Solid", "Slightly soluble"),
    "Melatonin":           ("117", D, D, "Solid", "Slightly soluble"),
    "Testosterone":        ("155", D, "1.17", "Solid", "Insoluble"),
    "Estradiol":           ("178", D, D, "Solid", "Insoluble"),
    "Vitamin C":           ("191 (dec.)", D, "1.65", "Solid", "Very soluble"),
    "Vitamin B3 (Niacin)": ("237", "sublimes", "1.47", "Solid", "Soluble"),
    "Vitamin A (Retinol)": ("63", D, D, "Solid", "Insoluble"),

    # ---- Common drugs ----
    "Aspirin":             ("135", "140 (dec.)", "1.40", "Solid", "Slightly soluble"),
    "Paracetamol":         ("169", D, "1.26", "Solid", "Slightly soluble"),
    "Ibuprofen":           ("76", "157 (4 mmHg)", "1.03", "Solid", "Practically insoluble"),
    "Naproxen":            ("153", "403.9", "1.26", "Solid", "Practically insoluble"),
    "Metformin":           ("224", "dec.", "1.36", "Solid", "Very soluble"),
    "Amoxicillin":         ("194 (dec.)", D, D, "Solid", "Slightly soluble"),
    "Penicillin G":        ("215 (dec.)", D, "1.41", "Solid", "Soluble"),
    "Ciprofloxacin":       ("256 (dec.)", D, D, "Solid", "Slightly soluble"),
    "Diazepam":            ("133", D, D, "Solid", "Practically insoluble"),
    "Fluoxetine":          ("158", "395", "1.20", "Solid", "Sparingly soluble"),
    "Omeprazole":          ("155 (dec.)", D, D, "Solid", "Slightly soluble"),
    "Lidocaine":           ("68", "181 (4 mmHg)", D, "Solid", "Slightly soluble"),
    "Propranolol":         ("96", "434.5", "1.11", "Solid", "Slightly soluble"),
    "Salbutamol":          ("157.5", "433.5", "1.23", "Solid", "Sparingly soluble"),
    "Cetirizine":          ("112", "dec.", "1.24", "Solid", "Soluble"),
    "Warfarin":            ("161", D, D, "Solid", "Practically insoluble"),
    "Tamoxifen":           ("97", "495", "1.09", "Solid", "Insoluble"),
    "Sildenafil":          ("188", "dec.", "1.39", "Solid", "Insoluble"),
    "Atorvastatin":        ("176", "dec.", "1.28", "Solid", "Very slightly soluble"),
    "Morphine":            ("255 (dec.)", D, "1.32", "Solid", "Slightly soluble"),

    # ---- Stimulants & food chemistry ----
    "Caffeine":            ("236", "178 (sublimes)", "1.23", "Solid", "Soluble"),
    "Theobromine":         ("357", "sublimes", D, "Solid", "Slightly soluble"),
    "Nicotine":            ("-79", "247", "1.01", "Liquid", "Miscible"),
    "Capsaicin":           ("63", D, D, "Solid", "Insoluble"),
    "Vanillin":            ("82", "285", "1.06", "Solid", "Slightly soluble"),
    "Saccharin":           ("229", D, "1.83", "Solid", "Slightly soluble"),
    "Aspartame":           ("246", D, D, "Solid", "Slightly soluble"),

    # ---- Natural products ----
    "Limonene":            ("-74", "176", "0.84", "Liquid", "Insoluble"),
    "Camphor":             ("176", "204 (sublimes)", "0.99", "Solid", "Slightly soluble"),
    "Eugenol":             ("-9", "254", "1.07", "Liquid", "Slightly soluble"),
    "Resveratrol":         ("254", D, D, "Solid", "Insoluble"),
    "Curcumin":            ("183", D, D, "Solid", "Insoluble"),
    "Quercetin":           ("316 (dec.)", D, D, "Solid", "Practically insoluble"),

    # ---- Industrial & environmental ----
    "Ethylene Oxide":      ("-112.5", "10.7", D, "Gas", "Miscible"),
    "Vinyl Chloride":      ("-153.8", "-13.4", D, "Gas", "Slightly soluble"),
    "Acrylonitrile":       ("-83.5", "77.3", "0.81", "Liquid", "Soluble"),
    "Caprolactam":         ("69", "270", "1.02", "Solid", "Very soluble"),
    "Bisphenol A":         ("158", "220 (4 mmHg)", "1.20", "Solid", "Slightly soluble"),
    "MTBE":                ("-109", "55.2", "0.74", "Liquid", "Slightly soluble"),
    "Freon-12 (CFC-12)":   ("-157.7", "-29.8", D, "Gas", "Very slightly soluble"),
    "DDT":                 ("108.5", "260 (dec.)", "0.99", "Solid", "Insoluble"),
    "TNT":                 ("80.5", "240 (dec.)", "1.65", "Solid", "Slightly soluble"),
}

# Cache canonical SMILES mapping to common names
_CANONICAL_INDEX = None


def _get_canonical_index():
    global _CANONICAL_INDEX
    if _CANONICAL_INDEX is None:
        index = {}
        try:
            from rdkit import Chem
            from molecules import MOLECULE_LIBRARY
            for cat, items in MOLECULE_LIBRARY.items():
                for name, smi in items.items():
                    try:
                        m = Chem.MolFromSmiles(smi)
                        if m:
                            can = Chem.MolToSmiles(m)
                            index[can] = name
                    except Exception:
                        pass
        except Exception:
            pass
        _CANONICAL_INDEX = index
    return _CANONICAL_INDEX


def get_experimental(name):
    """Return a dict of experimental properties for a molecule name, or None."""
    row = EXPERIMENTAL.get(name)
    if not row:
        return None
    mp, bp, density, state, solubility = row
    return {
        "name": name,
        "mp": mp,
        "bp": bp,
        "density": density,
        "state": state,
        "solubility": solubility,
        "source": "Curated Literature (CRC / PubChem)"
    }


def get_experimental_by_smiles(smiles, mol=None):
    """
    Given a SMILES or RDKit molecule, resolve its canonical SMILES and retrieve
    experimental data if present in our library. If not found, estimates basic
    physical trends from cheminformatics descriptors.
    """
    if mol is None and smiles:
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles.strip())
        except Exception:
            mol = None

    if mol is None:
        return None

    try:
        from rdkit import Chem
        can = Chem.MolToSmiles(mol)
        idx = _get_canonical_index()
        name = idx.get(can)
        if name and name in EXPERIMENTAL:
            data = get_experimental(name)
            return data
    except Exception:
        pass

    # Fallback heuristic prediction for custom molecules
    mw = mol.GetNumHeavyAtoms() * 14.0  # rough proxy if Descriptors isn't pre-imported
    return {
        "name": "Custom Structure",
        "mp": "— (PubChem required)",
        "bp": "— (PubChem required)",
        "density": "—",
        "state": "Solid (Predicted)" if mol.GetNumHeavyAtoms() > 10 else "Liquid / Gas",
        "solubility": "Estimate via LogP",
        "source": "Custom SMILES (Not in offline library)"
    }


_NUM = re.compile(r"-?\d+(?:\.\d+)?")


def as_number(text):
    """Extract the first numeric value from text like '233 (dec.)' or '-114.1'."""
    if not text:
        return None
    m = _NUM.search(str(text))
    return float(m.group()) if m else None


def analyze_intermolecular_forces(mol, properties, exp_data=None):
    """
    Analyze intermolecular forces (IMF) and structure-property relationships:
    - Hydrogen Bonding (HBD, HBA)
    - Dipole-Dipole interactions (TPSA, polar heteroatoms)
    - London Dispersion Forces (LDF, driven by MW and polarizability)
    - π-π Stacking (Aromatic rings)
    - Thermal physical state & boiling/melting point rationale
    """
    forces = []
    explanations = []

    hbd = properties.get("Hydrogen Bond Donors", 0)
    hba = properties.get("Hydrogen Bond Acceptors", 0)
    tpsa = properties.get("TPSA", 0.0)
    mw = properties.get("Molecular Weight", 0.0)
    aromatic_rings = properties.get("Aromatic Rings", 0)
    logp = properties.get("LogP", 0.0)

    # 1. Hydrogen Bonding
    if hbd > 0 and hba > 0:
        forces.append("Strong Hydrogen Bonding (Donor & Acceptor)")
        explanations.append(
            f"Contains {hbd} H-bond donor(s) (–OH / –NH) and {hba} acceptor(s) (O / N lone pairs). "
            "This creates strong intermolecular network attractions, significantly elevating the boiling point."
        )
    elif hba > 0:
        forces.append("H-Bond Acceptor Only (Can accept H-bonds from protic solvents like water)")
        explanations.append(
            f"Contains {hba} acceptor site(s) but no internal donor hydrogens. It cannot hydrogen-bond with itself, "
            "keeping neat boiling points lower than alcohols/acids of similar molecular weight."
        )

    # 2. Dipole-Dipole
    if tpsa >= 25.0:
        forces.append("Dipole–Dipole Attractions")
        explanations.append(
            f"The substantial polar surface area ({tpsa} Å²) gives rise to permanent molecular dipoles "
            "that align attractively in liquid and solid phases."
        )

    # 3. π-π Stacking
    if aromatic_rings > 0:
        forces.append("Aromatic π–π Stacking")
        explanations.append(
            f"Features {aromatic_rings} planar aromatic ring(s) that permit π-electron cloud overlap "
            "and crystal lattice stabilization (higher melting points)."
        )

    # 4. London Dispersion
    if mw > 150:
        forces.append("Extensive London Dispersion Forces (LDF)")
        explanations.append(
            f"Higher molecular mass ({mw} g/mol) and larger electron cloud create strong temporary dipole attractions."
        )
    else:
        forces.append("London Dispersion Forces (LDF)")
        explanations.append("Universal temporary dipole attractions proportional to molecular surface area.")

    # Solubility insight
    if logp < 0.5:
        solubility_note = "High water solubility favored by hydrophilic polar interactions."
    elif logp <= 2.5:
        solubility_note = "Balanced amphiphilic profile; moderate solubility in both polar and organic solvents."
    else:
        solubility_note = "Lipophilic/hydrophobic character; poor water solubility, prefers organic solvents."

    return {
        "forces": forces,
        "explanations": explanations,
        "solubility_note": solubility_note,
        "primary_force": forces[0] if forces else "London Dispersion Forces"
    }

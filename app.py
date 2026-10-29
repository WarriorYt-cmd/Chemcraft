import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO

# RDKit
from rdkit import Chem
from rdkit.Chem import (
    Draw,
    AllChem,
    Descriptors,
    Crippen,
    Lipinski,
    rdMolDescriptors
)

# 3D Visualization
import py3Dmol
import streamlit.components.v1 as components

# PDF
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from molecule_data import get_experimental
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.colors import HexColor
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    PageBreak,
    Image
)
from reportlab.lib.units import inch
from reportlab.lib.utils import ImageReader
import tempfile
import os


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="ChemCraft - Cheminformatics Lab",
    page_icon="⛏️",
    layout="wide"
)


# ============================================================
# CUSTOM CSS & ENHANCED TYPOGRAPHY
# ============================================================

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Press+Start+2P&family=VT323&display=swap');

/* Make the Main Page Title bigger and glowing */
.main-title {
    text-align: center;
    font-size: 46px !important;
    font-weight: 800 !important;
    color: #fcfc54 !important;
    text-shadow: 3px 3px 0px #3f3f15, 0px 0px 15px rgba(252, 252, 84, 0.5);
    margin-bottom: 10px;
}

/* Make Sidebar Title (ChemCraft Lab) much bigger and more noticeable */
[data-testid="stSidebar"] h1 {
    font-size: 28px !important;
    color: #55ff55 !important;
    text-shadow: 3px 3px 0px #153f15;
    font-family: 'Press Start 2P', monospace !important;
    letter-spacing: 1px;
    padding-bottom: 5px;
}

/* Enhance Sidebar headers / text readability */
[data-testid="stSidebar"] .stMarkdown h3 {
    font-size: 1.1rem !important;
    color: #ffaa00 !important;
    text-shadow: 2px 2px 0px #3f2a00;
}

/* Smooth Metric Card Styling */
[data-testid="stMetric"] {
    background: rgba(45, 45, 45, 0.8);
    border: 2px solid rgba(255, 255, 255, 0.2);
    border-radius: 8px;
    padding: 14px;
    box-shadow: 0 4px 12px rgba(0,0,0,0.3);
}

/* ============================================================
   ENHANCED SIDEBAR COLLAPSE / EXPAND BUTTON STYLING
   ============================================================ */
[data-testid="stSidebarCollapseButton"] button, 
[data-testid="collapsedControl"] button {
    background-color: #5d9b3a !important;
    border: 3px solid #3b6b22 !important;
    border-radius: 8px !important;
    box-shadow: 0 0 10px rgba(93, 161, 48, 0.6);
    transition: all 0.2s ease-in-out;
}

[data-testid="stSidebarCollapseButton"] svg, 
[data-testid="collapsedControl"] svg {
    fill: #fcfc54 !important;
    width: 24px !important;
    height: 24px !important;
}

[data-testid="stSidebarCollapseButton"] button:hover, 
[data-testid="collapsedControl"] button:hover {
    background-color: #7a4f2a !important;
    border-color: #fcfc54 !important;
    transform: scale(1.1);
}
</style>
""", unsafe_allow_html=True)

# ---- Minecraft theme (set to False for the plain look) ----
MINECRAFT_THEME = True

# ============================================================
# YOUR DETAILS
# ============================================================
APP_AUTHOR = "PRIYA SAHU (RA26011026050173)"
APP_ROLE = "Student"
APP_DEPARTMENT = "BTECH CSE AI ML"
APP_ORGANISATION = "SRMIST TRICHY"


if MINECRAFT_THEME:
    from minecraft_theme import apply_minecraft_theme
    apply_minecraft_theme()


# ============================================================
# SHARED MOLECULE INPUT STATE & HELPERS
# ============================================================

from molecules import MOLECULE_LIBRARY


@st.cache_data
def load_valid_library():
    """Keep only the molecules whose SMILES RDKit can parse."""
    valid = {}
    for category, items in MOLECULE_LIBRARY.items():
        good = {
            name: smi for name, smi in items.items()
            if Chem.MolFromSmiles(smi) is not None
        }
        if good:
            valid[category] = good
    return valid


VALID_LIBRARY = load_valid_library()

EXAMPLE_MOLECULES = {
    name: smi
    for category in VALID_LIBRARY.values()
    for name, smi in category.items()
}
if not EXAMPLE_MOLECULES:
    EXAMPLE_MOLECULES = {"Ethanol": "CCO"}

if "active_smiles" not in st.session_state:
    st.session_state.active_smiles = "CCO"

if "selected_example" not in st.session_state:
    st.session_state.selected_example = "Ethanol"

def load_shared_example():
    """Load the selected example into the shared SMILES field."""
    st.session_state.active_smiles = EXAMPLE_MOLECULES[
        st.session_state.selected_example
    ]

def get_molecule(smiles):
    """Convert SMILES into an RDKit molecule."""
    if not smiles or not smiles.strip():
        return None
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        return mol
    except Exception:
        return None

def get_molecule_name_from_smiles(target_smiles):
    """Find the library name of a molecule given its SMILES string."""
    for name, smi in EXAMPLE_MOLECULES.items():
        if smi.strip() == target_smiles.strip():
            return name
    return None


# ============================================================
# CALCULATE SIGMA AND PI BONDS
# ============================================================

def calculate_bond_information(mol):
    heavy_atom_bonds = mol.GetNumBonds()
    kekule_mol = Chem.Mol(mol)

    try:
        Chem.Kekulize(kekule_mol, clearAromaticFlags=True)
    except Exception:
        pass

    pi_bonds = 0
    for bond in kekule_mol.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type == Chem.BondType.DOUBLE:
            pi_bonds += 1
        elif bond_type == Chem.BondType.TRIPLE:
            pi_bonds += 2
        elif bond_type == Chem.BondType.QUADRUPLE:
            pi_bonds += 3
        elif bond.GetIsAromatic():
            pi_bonds += 0.5

    mol_with_h = Chem.AddHs(Chem.Mol(mol))
    sigma_bonds = mol_with_h.GetNumBonds()
    total_bonds = sigma_bonds + pi_bonds

    return {
        "Heavy Atom Bonds": heavy_atom_bonds,
        "Sigma Bonds": sigma_bonds,
        "Pi Bonds": pi_bonds,
        "Total Bonds": total_bonds,
        "Total Atoms Including H": mol_with_h.GetNumAtoms(),
    }


# ============================================================
# CALCULATE COMPLETE MOLECULAR PROPERTIES
# ============================================================
def calculate_properties(mol):
    bond_info = calculate_bond_information(mol)
    molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
    canonical_smi = Chem.MolToSmiles(mol)

    exp_data = None
    try:
        from molecule_data import EXPERIMENTAL, get_experimental
        from molecules import flat_library
        flat_lib = flat_library()
        matched_name = None
        for name, smi in flat_lib.items():
            if Chem.MolToSmiles(Chem.MolFromSmiles(smi)) == canonical_smi:
                matched_name = name
                break
        if matched_name:
            exp_data = get_experimental(matched_name)
    except Exception:
        pass

    properties = {
        "Molecular Formula": molecular_formula,
        "Molecular Weight": round(Descriptors.MolWt(mol), 3),
        "Exact Molecular Weight": round(Descriptors.ExactMolWt(mol), 4),
        "LogP": round(Crippen.MolLogP(mol), 3),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "Hydrogen Bond Donors": Lipinski.NumHDonors(mol),
        "Hydrogen Bond Acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Ring Count": Lipinski.RingCount(mol),
        "Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "Aliphatic Rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "Number of Atoms": mol.GetNumAtoms(),
        "Heavy Atoms": mol.GetNumHeavyAtoms(),
        "Heavy Atom Bonds": bond_info["Heavy Atom Bonds"],
        "Sigma (σ) Bonds": bond_info["Sigma Bonds"],
        "Pi (π) Bonds": bond_info["Pi Bonds"],
        "Total Bonds (σ + π)": bond_info["Total Bonds"],
        "Fraction Csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "Molar Refractivity": round(Crippen.MolMR(mol), 3),
        "Physical State (25°C)": exp_data["state"] if exp_data else "—",
        "Melting Point": f"{exp_data['mp']} °C" if exp_data and exp_data['mp'] != "—" else "—",
        "Boiling Point": f"{exp_data['bp']} °C" if exp_data and exp_data['bp'] != "—" else "—",
        "Density": f"{exp_data['density']} g/cm³" if exp_data and exp_data['density'] != "—" else "—",
        "Water Solubility": exp_data["solubility"] if exp_data else "—"
    }

    return properties


# ============================================================
# GENERATE 2D & 3D STRUCTURES
# ============================================================

def generate_2d_image(mol):
    image = Draw.MolToImage(mol, size=(700, 500))
    return image


def generate_3d_molecule(mol):
    try:
        mol3d = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(mol3d, params)
        if status != 0:
            return None
        try:
            AllChem.MMFFOptimizeMolecule(mol3d)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol3d)
            except Exception:
                pass
        return mol3d
    except Exception:
        return None


def display_3d_structure(mol3d, style, show_surface=False, enable_spin=False):
    mol_block = Chem.MolToMolBlock(mol3d)
    view = py3Dmol.view(width=900, height=600)
    view.addModel(mol_block, "mol")

    if style == "Ball and Stick":
        view.setStyle({"stick": {"radius": 0.15}, "sphere": {"scale": 0.30}})
    elif style == "Stick":
        view.setStyle({"stick": {"radius": 0.20}})
    elif style == "Space Filling":
        view.setStyle({"sphere": {"scale": 1.0}})
    elif style == "Wireframe":
        view.setStyle({"line": {"linewidth": 2}})
    elif style == "Cartoon":
        view.setStyle({"stick": {}, "sphere": {"scale": 0.20}})

    if show_surface:
        view.addSurface(py3Dmol.VDW, {"opacity": 0.7})

    view.setBackgroundColor("black" if MINECRAFT_THEME else "white")
    view.zoomTo()

    if enable_spin:
        view.spin(True)

    html = view._make_html()
    components.html(html, height=620, scrolling=False)


# ============================================================
# CREATE MOLECULAR PDF REPORT
# ============================================================

def generate_molecule_pdf(input_smiles, mol, properties):
    buffer = BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        leftMargin=45,
        rightMargin=45,
        topMargin=45,
        bottomMargin=45
    )

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "TitleStyle",
        parent=styles["Title"],
        fontSize=20,
        alignment=TA_CENTER,
        textColor=HexColor("#17365D")
    )
    heading_style = ParagraphStyle(
        "HeadingStyle",
        parent=styles["Heading1"],
        fontSize=14,
        textColor=HexColor("#17365D")
    )
    body_style = ParagraphStyle(
        "BodyStyle",
        parent=styles["BodyText"],
        fontSize=10,
        leading=14,
        alignment=TA_JUSTIFY
    )

    story = []
    story.append(Paragraph("MOLECULAR PROPERTIES REPORT", title_style))
    story.append(Spacer(1, 0.3 * inch))

    story.append(Paragraph("1. Molecular Information", heading_style))
    canonical_smiles = Chem.MolToSmiles(mol)
    basic_data = [
        ["Property", "Value"],
        ["Input SMILES", input_smiles],
        ["Canonical SMILES", canonical_smiles],
        ["Molecular Formula", properties["Molecular Formula"]]
    ]

    basic_table = Table(basic_data, colWidths=[2 * inch, 4 * inch])
    basic_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#17365D")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("PADDING", (0, 0), (-1, -1), 7)
        ])
    )
    story.append(basic_table)
    story.append(Spacer(1, 0.3 * inch))

    story.append(Paragraph("2. Molecular Properties", heading_style))
    property_data = [["Property", "Value"]]
    for key, value in properties.items():
        property_data.append([key, str(value)])

    property_table = Table(property_data, colWidths=[3 * inch, 3 * inch], repeatRows=1)
    property_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#2F5597")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("PADDING", (0, 0), (-1, -1), 6)
        ])
    )
    story.append(property_table)
    story.append(Spacer(1, 0.3 * inch))

    story.append(Paragraph("3. Interpretation", heading_style))
    interpretation = """
    This report was generated automatically using the
    Cheminformatics Virtual Laboratory. The molecular properties
    were calculated from the chemical structure represented by
    the submitted SMILES notation using the RDKit cheminformatics
    toolkit.
    """
    story.append(Paragraph(interpretation, body_style))

    doc.build(story)
    pdf_data = buffer.getvalue()
    buffer.close()
    return pdf_data


# ============================================================
# REPORT IMAGE HELPERS
# ============================================================

def save_2d_structure_png(mol):
    img = Draw.MolToImage(mol, size=(700, 500))
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    tmp.close()
    img.save(tmp.name)
    return tmp.name


def save_3d_structure_png(mol):
    import matplotlib.pyplot as plt
    import numpy as np

    mol3d = generate_3d_molecule(mol)
    if mol3d is None:
        return None

    conf = mol3d.GetConformer()
    n_atoms = mol3d.GetNumAtoms()

    coords = np.array([
        [
            conf.GetAtomPosition(i).x,
            conf.GetAtomPosition(i).y,
            conf.GetAtomPosition(i).z,
        ]
        for i in range(n_atoms)
    ], dtype=float)

    centered = coords - coords.mean(axis=0)
    if n_atoms >= 3 and np.linalg.matrix_rank(centered) >= 2:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        projected = centered @ vh[:2].T
    else:
        projected = centered[:, :2]

    cpk = {
        "H": "#FFFFFF", "C": "#3F3F3F", "N": "#2F5BFF",
        "O": "#E53935", "F": "#65C466", "P": "#F39C12",
        "S": "#F4D03F", "Cl": "#27AE60", "Br": "#8E3B2F",
        "I": "#8E44AD",
    }
    sizes = {
        "H": 70, "C": 240, "N": 270, "O": 285, "F": 260,
        "P": 300, "S": 310, "Cl": 310, "Br": 330, "I": 350,
    }

    fig, ax = plt.subplots(figsize=(8.8, 6.4), facecolor="white")
    ax.set_facecolor("white")

    for bond in mol3d.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        x1, y1 = projected[i]
        x2, y2 = projected[j]

        dx, dy = x2 - x1, y2 - y1
        length = float(np.hypot(dx, dy))
        if length == 0:
            continue

        px, py = -dy / length, dx / length
        btype = bond.GetBondType()
        if btype == Chem.BondType.DOUBLE:
            offsets = (-0.035, 0.035)
        elif btype == Chem.BondType.TRIPLE:
            offsets = (-0.055, 0.0, 0.055)
        else:
            offsets = (0.0,)

        for off in offsets:
            ax.plot(
                [x1 + px * off, x2 + px * off],
                [y1 + py * off, y2 + py * off],
                color="#707070", linewidth=3.0,
                solid_capstyle="round", zorder=1
            )

    for i, atom in enumerate(mol3d.GetAtoms()):
        symbol = atom.GetSymbol()
        x, y = projected[i]
        ax.scatter(
            x, y,
            s=sizes.get(symbol, 250),
            c=cpk.get(symbol, "#A0A0A0"),
            edgecolors="#222222",
            linewidths=1.0,
            zorder=3
        )

        if symbol not in ("H", "C"):
            ax.text(
                x, y, symbol,
                ha="center", va="center",
                fontsize=8, fontweight="bold",
                color="black", zorder=4
            )

    xmin, ymin = projected.min(axis=0)
    xmax, ymax = projected.max(axis=0)
    width = max(xmax - xmin, 1.0)
    height = max(ymax - ymin, 1.0)
    pad = 0.22 * max(width, height)

    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ax.set_title("3D Molecular Structure (Ball-and-Stick)", fontsize=18, fontweight="bold", pad=18)

    fig.tight_layout(pad=1.2)

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    tmp.close()
    fig.savefig(
        tmp.name,
        dpi=320,
        bbox_inches="tight",
        facecolor="white",
        pad_inches=0.12
    )
    plt.close(fig)
    return tmp.name


# ============================================================
# PROPERTY INTERPRETATION HELPERS
# ============================================================

def get_property_interpretation(property_name, value, properties):
    if property_name == "Molecular Formula":
        return "Shows the elemental composition of the molecule and is useful for understanding its basic chemical composition."
    if property_name == "Molecular Weight":
        return "Represents the average molecular mass. Molecular size and mass can influence diffusion, formulation and many physicochemical behaviours."
    if property_name == "Exact Molecular Weight":
        return "The monoisotopic mass calculated from exact isotopic masses; this is especially useful when comparing mass-spectrometry data."
    if property_name == "LogP":
        if value < 0:
            return "A negative value indicates a preference toward the aqueous phase and relatively high hydrophilicity."
        elif value < 1:
            return "Indicates low to moderate lipophilicity, suggesting a reasonable preference for polar environments."
        elif value < 3:
            return "Indicates moderate lipophilicity, often reflecting a balance between hydrophilic and hydrophobic character."
        else:
            return "Indicates relatively high lipophilicity and a stronger preference for non-polar environments."
    if property_name == "TPSA":
        if value < 60:
            return "A relatively low polar surface area, indicating limited overall molecular polarity."
        elif value < 140:
            return "A moderate polar surface area, reflecting the contribution of heteroatoms to molecular polarity."
        else:
            return "A high polar surface area, indicating substantial polarity and strong hydrogen-bonding potential."
    if property_name == "Hydrogen Bond Donors":
        return f"The molecule has {value} donor site(s) capable of donating hydrogen in hydrogen-bond interactions."
    if property_name == "Hydrogen Bond Acceptors":
        return f"The molecule has {value} acceptor site(s) capable of accepting hydrogen in hydrogen-bond interactions."
    if property_name == "Rotatable Bonds":
        if value == 0:
            return "No freely rotatable bonds are detected, suggesting a relatively rigid molecular framework."
        return f"{value} rotatable bond(s) indicate conformational flexibility in the molecular structure."
    if property_name == "Ring Count":
        return f"The molecule contains {value} ring(s), which contribute to its overall shape and structural rigidity."
    if property_name == "Aromatic Rings":
        return f"{value} aromatic ring(s) contribute to aromaticity, planarity and delocalized π-electron character."
    if property_name == "Aliphatic Rings":
        return f"{value} non-aromatic ring(s) contribute to cyclic structure without aromatic π-electron delocalization."
    if property_name == "Number of Atoms":
        return "Counts atoms explicitly present in the RDKit molecular representation; implicit hydrogens are not included here."
    if property_name == "Heavy Atoms":
        return "Counts all non-hydrogen atoms and provides a useful measure of molecular structural size."
    if property_name == "Heavy Atom Bonds":
        return "Counts connections between non-hydrogen atoms in the molecular graph."
    if property_name == "Sigma (σ) Bonds":
        return "Every bonded atom pair contains one σ bond; this value includes bonds involving explicit or added hydrogens."
    if property_name == "Pi (π) Bonds":
        return "π bonds arise from multiple bonds and aromatic systems and are associated with electron delocalization."
    if property_name == "Total Bonds (σ + π)":
        return "Represents the total bond components obtained by adding the calculated σ and π bond counts."
    if property_name == "Fraction Csp3":
        if value < 0.25:
            return "A low sp³ fraction suggests a comparatively planar, unsaturated or aromatic structural character."
        elif value < 0.60:
            return "An intermediate sp³ fraction indicates a balance between saturated and unsaturated structural features."
        return "A high sp³ fraction indicates a more saturated, three-dimensional molecular framework."
    if property_name == "Molar Refractivity":
        return "Relates to molecular volume and electronic polarizability and reflects how the electron cloud responds to an electric field."
    return "This descriptor provides additional quantitative information about the molecular structure and physicochemical behaviour."


def build_key_interpretation_points(properties):
    points = []
    points.append(
        f"<b>• Molecular identity:</b> The submitted structure has molecular formula "
        f"<b>{properties['Molecular Formula']}</b> and molecular weight "
        f"<b>{properties['Molecular Weight']}</b>."
    )

    logp = properties["LogP"]
    if logp < 0:
        lipophilic_text = "more hydrophilic than lipophilic"
    elif logp < 1:
        lipophilic_text = "low in lipophilicity"
    elif logp < 3:
        lipophilic_text = "moderately lipophilic"
    else:
        lipophilic_text = "strongly lipophilic"

    points.append(
        f"<b>• Lipophilicity:</b> LogP = <b>{logp}</b>, indicating that the molecule is "
        f"<b>{lipophilic_text}</b>. This descriptor describes the balance between "
        f"preference for non-polar and polar environments."
    )

    points.append(
        f"<b>• Molecular polarity:</b> TPSA = <b>{properties['TPSA']}</b> Å². "
        f"TPSA reflects the contribution of polar atoms and functional groups to the "
        f"overall surface polarity of the molecule."
    )

    points.append(
        f"<b>• Hydrogen bonding:</b> The structure contains "
        f"<b>{properties['Hydrogen Bond Donors']}</b> hydrogen-bond donor(s) and "
        f"<b>{properties['Hydrogen Bond Acceptors']}</b> acceptor(s), which helps "
        f"describe possible intermolecular interactions."
    )

    points.append(
        f"<b>• Flexibility:</b> {properties['Rotatable Bonds']} rotatable bond(s) "
        f"indicate the degree of conformational freedom available to the molecule."
    )

    points.append(
        f"<b>• Cyclic and aromatic character:</b> The molecule contains "
        f"{properties['Ring Count']} ring(s), including "
        f"{properties['Aromatic Rings']} aromatic ring(s). Aromatic rings are "
        f"associated with delocalized π-electron systems."
    )

    points.append(
        f"<b>• Bonding pattern:</b> The calculated structure contains "
        f"<b>{properties['Sigma (σ) Bonds']}</b> σ bonds and "
        f"<b>{properties['Pi (π) Bonds']}</b> π bonds. This helps distinguish the "
        f"framework-forming σ bonds from multiple-bond and delocalized π bonding."
    )

    points.append(
        f"<b>• Three-dimensional character:</b> Fraction Csp³ = "
        f"<b>{properties['Fraction Csp3']}</b>. This descriptor provides a simple "
        f"indication of the balance between saturated 3D carbon environments and "
        f"planar/unsaturated environments."
    )

    points.append(
        "<b>• Overall conclusion:</b> The calculated descriptors should be interpreted "
        "together with the 2D connectivity and generated 3D geometry to understand "
        "the molecule's composition, bonding, polarity, and flexibility."
    )
    return points


# ============================================================
# CHEMBOT AI ASSISTANT SIDEBAR WIDGET
# ============================================================

st.sidebar.divider()
st.sidebar.markdown("### 🤖 ChemBot AI Assistant")

from ai_assistant import ask_ai

if "ai_messages" not in st.session_state:
    st.session_state.ai_messages = []

with st.sidebar.expander("💬 Chat with ChemBot", expanded=False):
    st.caption("Ask chemistry or app-related questions!")
    
    for msg in st.session_state.ai_messages:
        if msg["role"] == "user":
            st.markdown(f"**You:** {msg['content']}")
        else:
            st.markdown(f"**ChemBot:** {msg['content']}")

    user_query = st.text_input("Ask ChemBot anything...", key="chembot_input_field")
    
    if st.button("Send to ChemBot"):
        if user_query.strip():
            st.session_state.ai_messages.append({"role": "user", "content": user_query})
            active_mol_context = ""
            if "latest_smiles" in st.session_state:
                active_mol_context = f"SMILES: {st.session_state.latest_smiles}"
            
            success, reply = ask_ai(st.session_state.ai_messages, molecule_context=active_mol_context)
            if success:
                st.session_state.ai_messages.append({"role": "assistant", "content": reply})
                st.rerun()
            else:
                st.error(reply)


# ============================================================
# SIDEBAR NAVIGATION
# ============================================================

st.sidebar.markdown("---")
st.sidebar.markdown("### 🧭 Navigation")

page = st.sidebar.radio(
    "Select Lab Section",
    [
        "📖 Theory",
        "🧬 Molecular Visualization",
        "📊 Complete Properties Report",
        "📥 Export Reports"
    ],
    key="lab_navigation"
)

st.sidebar.markdown("---")
st.sidebar.markdown(f"**Author:** {APP_AUTHOR}")
st.sidebar.markdown(f"**Role:** {APP_ROLE}")
st.sidebar.markdown(f"**Dept:** {APP_DEPARTMENT}")
st.sidebar.markdown(f"**Inst:** {APP_ORGANISATION}")


# ============================================================
# THEORY PAGE
# ============================================================

if page == "📖 Theory":

    st.title("📖 Introduction to Cheminformatics")

    st.markdown("""
## What is Cheminformatics?

Cheminformatics combines:

- Chemistry ⚗️
- Computer Science 💻
- Data Analysis 📊

It uses computational techniques to store, process,
visualize and analyze chemical information.
""")

    st.markdown("""
## SMILES

**SMILES** means:

### Simplified Molecular Input Line Entry System

It represents a molecular structure using a text string.

Examples:

| Molecule | SMILES |
|---|---|
| Ethanol | CCO |
| Benzene | c1ccccc1 |
| Acetic Acid | CC(=O)O |
""")

    st.markdown("""
## Molecular Descriptors

Molecular descriptors provide numerical information
about a chemical structure.

Examples include:

- Molecular Weight
- LogP
- TPSA
- Hydrogen Bond Donors
- Hydrogen Bond Acceptors
- Rotatable Bonds
""")

    st.markdown("""
## 3D Molecular Models

Different visualization models help students understand
the molecular geometry.

### 🟢 Ball and Stick

Shows atoms as spheres and bonds as sticks.

### 🟡 Stick Model

Emphasizes chemical bonding.

### 🔵 Space Filling

Shows the approximate molecular volume.

### ⚫ Wireframe

Provides a simplified representation of the molecule.
""")


# ============================================================
# MOLECULAR VISUALIZATION PAGE
# ============================================================

elif page == "🧬 Molecular Visualization":

    st.title("🧬 Molecular Visualization Laboratory")
    st.info("Enter any SMILES notation or load an example. The interactive 3D rotating model and 2D structure will generate automatically.")

    col1, col2 = st.columns([2, 1])

    with col1:
        smiles = st.text_input(
            "Enter SMILES",
            key="active_smiles"
        )

    with col2:
        st.selectbox(
            "Choose Example Molecule",
            list(EXAMPLE_MOLECULES.keys()),
            key="selected_example"
        )

        st.button(
            "Load Example",
            on_click=load_shared_example,
            key="load_visualization_example"
        )

    mol = get_molecule(smiles)

    if mol is None:
        st.error("❌ Invalid SMILES notation.")
    else:
        st.success("✅ Valid Molecular Structure!")
        st.session_state.latest_smiles = smiles
        st.session_state.latest_properties = calculate_properties(mol)

        col_1, col_2 = st.columns([1.2, 1])

        with col_1:
            st.subheader("2D Molecular Structure")
            image = generate_2d_image(mol)
            st.image(image)

        with col_2:
            st.subheader("Molecular Information")
            st.write("**Input SMILES:**", smiles)
            st.write("**Canonical SMILES:**", Chem.MolToSmiles(mol))
            st.write("**Formula:**", rdMolDescriptors.CalcMolFormula(mol))

            bond_info = calculate_bond_information(mol)
            st.write("**Heavy Atoms:**", mol.GetNumHeavyAtoms())
            st.write("**Total Atoms (incl. H):**", bond_info["Total Atoms Including H"])
            st.write("**Sigma (σ) Bonds:**", bond_info["Sigma Bonds"])
            st.write("**Pi (π) Bonds:**", bond_info["Pi Bonds"])

        st.divider()

        with st.spinner("Generating interactive 3D model..."):
            mol3d = generate_3d_molecule(mol)

        if mol3d is None:
            st.error("Unable to generate the 3D structure for this molecule.")
        else:
            st.subheader("🌐 Interactive 3D Rotating Molecular Model")
            st.caption("🖱️ Click and drag your mouse over the 3D viewer to rotate manually, or check the box below to spin it automatically!")

            col_3d_1, col_3d_2, col_3d_3 = st.columns(3)
            with col_3d_1:
                visualization_style = st.selectbox(
                    "Select Model Style",
                    ["Ball and Stick", "Stick", "Space Filling", "Wireframe", "Cartoon"],
                    key="vis_3d_style"
                )
            with col_3d_2:
                show_surface = st.checkbox("Show van der Waals Surface", key="vis_surface")
            with col_3d_3:
                enable_spin = st.checkbox("🔄 Auto-Rotate / Spin", value=True, key="vis_spin")

            display_3d_structure(mol3d, visualization_style, show_surface, enable_spin)


# ============================================================
# COMPLETE PROPERTIES REPORT PAGE
# ============================================================

elif page == "📊 Complete Properties Report":

    st.title("📊 Complete Molecular Properties Report")

    if "latest_smiles" not in st.session_state or not st.session_state.latest_smiles:
        st.warning("⚠️ Please analyze a molecule first in the 'Molecular Visualization' page.")
    else:
        smiles = st.session_state.latest_smiles
        mol = get_molecule(smiles)

        if mol is not None:
            properties = calculate_properties(mol)
            st.session_state.latest_properties = properties

            st.write(f"**Current Molecule SMILES:** `{smiles}`")
            st.write(f"**Canonical SMILES:** `{Chem.MolToSmiles(mol)}`")

            df_props = pd.DataFrame(
                list(properties.items()),
                columns=["Property Descriptor", "Calculated Value"]
            )
            st.dataframe(df_props, use_container_width=True, hide_index=True)

            st.markdown("### 📝 Detailed Property Explanations")
            for prop_name, prop_val in properties.items():
                explanation = get_property_interpretation(prop_name, prop_val, properties)
                st.markdown(f"- **{prop_name} (`{prop_val}`):** {explanation}")


# ============================================================
# EXPORT REPORTS PAGE
# ============================================================

elif page == "📥 Export Reports":

    st.title("📥 Export Molecular Reports")

    if "latest_smiles" not in st.session_state or not st.session_state.latest_smiles:
        st.warning("⚠️ Please analyze a molecule first in the 'Molecular Visualization' page.")
    else:
        smiles = st.session_state.latest_smiles
        mol = get_molecule(smiles)

        if mol is not None:
            properties = calculate_properties(mol)

            col_exp_1, col_exp_2 = st.columns(2)

            with col_exp_1:
                st.subheader("📄 Download PDF Report")
                st.write("Generate a formal professional PDF document containing summary tables and calculated descriptors.")

                if st.button("🚀 Build PDF Report"):
                    with st.spinner("Compiling PDF document..."):
                        pdf_bytes = generate_molecule_pdf(smiles, mol, properties)
                        st.download_button(
                            label="📥 Download PDF File",
                            data=pdf_bytes,
                            file_name=f"ChemCraft_Report_{properties['Molecular Formula']}.pdf",
                            mime="application/pdf"
                        )

            with col_exp_2:
                st.subheader("📊 Download CSV Data")
                st.write("Export all calculated numerical properties into a comma-separated values file.")

                df_export = pd.DataFrame(
                    list(properties.items()),
                    columns=["Property", "Value"]
                )
                csv_data = df_export.to_csv(index=False).encode("utf-8")

                st.download_button(
                    label="📥 Download CSV File",
                    data=csv_data,
                    file_name=f"ChemCraft_Data_{properties['Molecular Formula']}.csv",
                    mime="text/csv"
                )

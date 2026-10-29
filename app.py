# ============================================================
# SIDEBAR
# ============================================================

st.sidebar.title("⛏️ ChemCraft Lab")

page = st.sidebar.radio(
    "Navigation",
    [
        "🏠 Home",
        "📖 Theory",
        "🧬 Molecular Visualization",
        "📊 Molecular Properties Report",
        "📈 Structure–Property Analysis",
        "📚 100+ Molecule Master Table",
    ]
)

st.sidebar.divider()

st.sidebar.markdown("""
### WHAT WILL YOU GET :

- SMILES notation
- Molecular representation
- 2D visualization
- 3D visualization
- Molecular descriptors
- Structure-property relationships
""")


# ============================================================
# HOME PAGE
# ============================================================

if page == "🏠 Home":

    logo_path = os.path.join(os.path.dirname(__file__), "logo.jpeg")

    if os.path.exists(logo_path):
        logo_col1, logo_col2, logo_col3 = st.columns([1, 4, 1])
        with logo_col2:
            st.image(logo_path, use_container_width=True)

    st.markdown(
        '<div class="main-title">⛏️ ChemCraft: Cheminformatics Lab</div>',
        unsafe_allow_html=True
    )

    st.markdown(
        f'<div class="mc-subtitle">Mine molecules. Craft knowledge.</div>'
        f'<div style="text-align:center; margin-top:8px;">'
        f'Built by {APP_AUTHOR} | {APP_ORGANISATION}</div>',
        unsafe_allow_html=True
    )

    st.divider()

    st.markdown("## Welcome, Crafter!")

    st.write("""
    This interactive virtual laboratory helps students understand
    important concepts in cheminformatics and computational chemistry.
    """)

    st.markdown("### 🧪 What can students do?")

    st.markdown("""
    ✅ Enter any SMILES notation

    ✅ Validate molecular structures

    ✅ Generate 2D structures

    ✅ Generate interactive 3D structures

    ✅ Select different molecular visualization styles

    ✅ Calculate complete molecular properties

    ✅ Download molecular reports
    """)

    st.markdown("### Example Molecules")

    example_df = pd.DataFrame({
        "Molecule": [
            "Methanol",
            "Toluene",
            "Aniline",
            "Benzoic Acid",
            "Paracetamol",
            "Glucose"
        ],
        "SMILES": [
            "CO",
            "Cc1ccccc1",
            "Nc1ccccc1",
            "CC(=O)O",
            "c1ccc(cc1)C(=O)O",
            "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1OC"
        ]
    })

    st.dataframe(
        example_df,
        use_container_width=True,
        hide_index=True
    )

    st.divider()

    st.markdown("### 👨‍💻 Developed by")
    st.markdown(f"**{APP_AUTHOR}**")
    st.markdown(f"{APP_ROLE}  ")
    st.markdown(f"{APP_DEPARTMENT}  ")
    st.markdown(f"{APP_ORGANISATION}")


# ============================================================
# THEORY PAGE
# ============================================================

elif page == "📖 Theory":

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
# COMPLETE MOLECULAR PROPERTIES REPORT
# ============================================================

elif page == "📊 Molecular Properties Report":

    st.title("📊 Complete Molecular Properties Report")

    st.info("""
Enter any valid SMILES notation. The system will automatically
generate a complete molecular properties report.
""")

    st.write("### Select or Enter Your Molecule")

    col1, col2 = st.columns([2, 1])

    with col1:
        smiles = st.text_input(
            "Enter New SMILES",
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
            key="load_properties_example"
        )

    st.caption("Once you enter a valid SMILES or load an example, the complete 2D/3D structures and molecular properties below are generated for the same molecule.")

    if smiles:
        mol = get_molecule(smiles)

        if mol is None:
            st.error("❌ Invalid SMILES notation. Please check the structure.")
        else:
            st.success("✅ Valid Molecular Structure Detected!")

            properties = calculate_properties(mol)
            canonical_smiles = Chem.MolToSmiles(mol)

            st.session_state.latest_smiles = smiles
            st.session_state.latest_properties = properties

            col1, col2 = st.columns([1.2, 1])

            with col1:
                st.subheader("🧬 2D Molecular Structure")
                image = generate_2d_image(mol)
                st.image(image)

            with col2:
                st.subheader("📋 Basic Information")
                st.write("**Input SMILES:**", smiles)
                st.write("**Canonical SMILES:**", canonical_smiles)
                st.write("**Molecular Formula:**", properties["Molecular Formula"])
                st.write("**Atoms:**", properties["Number of Atoms"])
                st.write("**Heavy Atoms:**", properties["Heavy Atoms"])
                st.write("**Heavy-Atom Bonds:**", properties["Heavy Atom Bonds"])
                st.write("**Sigma (σ) Bonds:**", properties["Sigma (σ) Bonds"])
                st.write("**Pi (π) Bonds:**", properties["Pi (π) Bonds"])
                st.write("**Total Bonds (σ + π):**", properties["Total Bonds (σ + π)"])

            st.divider()

            st.subheader("📊 Important Molecular Properties")

            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Molecular Weight", f"{properties['Molecular Weight']} g/mol")
            col2.metric("LogP", properties["LogP"])
            col3.metric("TPSA", f"{properties['TPSA']} Å²")
            col4.metric("Ring Count", properties["Ring Count"])

            col1, col2, col3, col4 = st.columns(4)
            col1.metric("H-Bond Donors", properties["Hydrogen Bond Donors"])
            col2.metric("H-Bond Acceptors", properties["Hydrogen Bond Acceptors"])
            col3.metric("Rotatable Bonds", properties["Rotatable Bonds"])
            col4.metric("Fraction Csp3", properties["Fraction Csp3"])

            st.divider()

            st.subheader("📑 Complete Molecular Property Table")

            report_data = [
                ["Input SMILES", smiles],
                ["Canonical SMILES", canonical_smiles]
            ]
            for key, value in properties.items():
                report_data.append([key, value])

            report_df = pd.DataFrame(
                report_data,
                columns=["Property", "Value"]
            )

            st.dataframe(
                report_df,
                use_container_width=True,
                hide_index=True
            )

            csv = report_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "⬇️ Download Properties Report (CSV)",
                data=csv,
                file_name="Molecular_Properties_Report.csv",
                mime="text/csv"
            )

            if st.button("📄 Generate Molecular PDF Report"):
                with st.spinner("Generating molecular PDF report..."):
                    pdf = generate_molecule_pdf(smiles, mol, properties)

                st.download_button(
                    "⬇️ Download Molecular Properties PDF",
                    data=pdf,
                    file_name="Molecular_Properties_Report.pdf",
                    mime="application/pdf"
                )

            st.divider()

            st.subheader("🌐 3D Molecular Visualization")
            model_style = st.selectbox(
                "Select 3D Model",
                ["Ball and Stick", "Stick", "Space Filling", "Wireframe", "Cartoon"],
                key="properties_3d_style"
            )
            surface = st.checkbox("Show van der Waals Surface", key="properties_surface")

            if st.button("🚀 Generate 3D Molecular Model", key="properties_generate_3d"):
                with st.spinner("Generating 3D molecular model..."):
                    mol3d = generate_3d_molecule(mol)

                if mol3d is None:
                    st.error("Unable to generate the 3D structure.")
                else:
                    display_3d_structure(mol3d, model_style, surface)


# ============================================================
# STRUCTURE PROPERTY ANALYSIS
# ============================================================

elif page == "📈 Structure–Property Analysis":

    st.title("📈 Structure–Property Analysis")

    st.write("""
Compare the properties of multiple molecules.
Enter one molecule per line in the format:

**Molecule Name, SMILES**
""")

    STARTER = {
        "Ethanol": "CCO",
        "Benzene": "c1ccccc1",
        "Phenol": "Oc1ccccc1",
        "Acetic Acid": "CC(=O)O",
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    }

    preset = st.selectbox(
        "Quick-load a molecule set",
        ["Starter set"] + list(VALID_LIBRARY.keys()) + ["Blank"],
        key="analysis_preset"
    )

    if preset == "Starter set":
        preset_items = STARTER
    elif preset == "Blank":
        preset_items = {}
    else:
        preset_items = VALID_LIBRARY[preset]

    default_data = "\n".join(f"{n},{s}" for n, s in preset_items.items())

    molecule_input = st.text_area(
        "Enter Molecules",
        value=default_data,
        height=200
    )

    rows = []
    for line in molecule_input.splitlines():
        if "," not in line:
            continue
        name, smiles = line.split(",", 1)
        mol = get_molecule(smiles)
        if mol:
            prop = calculate_properties(mol)
            rows.append({
                "Molecule": name.strip(),
                "MW": prop["Molecular Weight"],
                "LogP": prop["LogP"],
                "TPSA": prop["TPSA"],
                "HBD": prop["Hydrogen Bond Donors"],
                "HBA": prop["Hydrogen Bond Acceptors"],
                "Rotatable Bonds": prop["Rotatable Bonds"],
                "Ring Count": prop["Ring Count"]
            })

    if rows:
        df = pd.DataFrame(rows)
        st.dataframe(
            df,
            use_container_width=True,
            hide_index=True
        )

        st.subheader("📋 Comparison Table")

        csv_data = df.to_csv(index=False).encode("utf-8")
        st.download_button(
            "⬇️ Download Comparison Table (CSV)",
            data=csv_data,
            file_name="structure_property_comparison.csv",
            mime="text/csv",
            key="download_structure_property_csv"
        )


# ============================================================
# 100+ MOLECULE MASTER TABLE & PRINTABLE VIEW
# ============================================================

elif page == "📚 100+ Molecule Master Table":

    st.title("📚 100+ Molecule Master Reference Table")
    st.info("This comprehensive table compiles all 100+ molecules from your library, combining RDKit calculated descriptors with experimental physical properties. Use the buttons below to export or print the entire 100+ molecule dataset in one go.")

    @st.cache_data
    def build_master_table_data():
        master_rows = []
        for category, items in VALID_LIBRARY.items():
            for name, smiles in items.items():
                mol = get_molecule(smiles)
                if mol:
                    props = calculate_properties(mol)
                    exp = get_experimental(name) or {}
                    master_rows.append({
                        "Category": category,
                        "Molecule": name,
                        "SMILES": smiles,
                        "Formula": props["Molecular Formula"],
                        "MW (g/mol)": props["Molecular Weight"],
                        "LogP": props["LogP"],
                        "TPSA (Å²)": props["TPSA"],
                        "HBD": props["Hydrogen Bond Donors"],
                        "HBA": props["Hydrogen Bond Acceptors"],
                        "State (25°C)": exp.get("state", "—"),
                        "Melting Pt (°C)": exp.get("mp", "—"),
                        "Boiling Pt (°C)": exp.get("bp", "—"),
                    })
        return pd.DataFrame(master_rows)

    with st.spinner("Compiling 100+ molecule library data..."):
        master_df = build_master_table_data()

    st.success(f"Successfully loaded {len(master_df)} molecules across {len(master_df['Category'].unique())} categories!")

    search_query = st.text_input("🔍 Search Master Table (by Molecule Name or Category)", "")
    if search_query:
        filtered_df = master_df[
            master_df["Molecule"].str.contains(search_query, case=False, na=False) |
            master_df["Category"].str.contains(search_query, case=False, na=False)
        ]
    else:
        filtered_df = master_df

    st.dataframe(filtered_df, use_container_width=True, hide_index=True)

    col1, col2 = st.columns(2)
    with col1:
        csv_master = master_df.to_csv(index=False).encode("utf-8")
        st.download_button(
            "⬇️ Download Complete 100+ Molecule Master Table (CSV)",
            data=csv_master,
            file_name="ChemCraft_100_Molecule_Master_Table.csv",
            mime="text/csv"
        )
    with col2:
        if st.button("🖨️ Print Complete Table"):
            st.markdown(
                """
                <script>
                window.print();
                </script>
                """,
                unsafe_allow_html=True
            )
            st.info("Print dialog triggered! You can save as a PDF or print all 100+ rows directly.")


# ============================================================
# FOOTER
# ============================================================

st.divider()

st.caption(
    "⛏️ ChemCraft Cheminformatics Lab | "
    "Developed using Python, Streamlit and RDKit"
)

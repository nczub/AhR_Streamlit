import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from mordred import Calculator, descriptors
from supervised.automl import AutoML
from streamlit_ketcher import st_ketcher
import matplotlib.pyplot as plt


st.markdown("""
<style>
    [data-testid=stSidebar] {
        background-color: #a2c7ab  ;
    }
</style>
""", unsafe_allow_html=True)

success_style = """
    background-color: #a2c7ab;
    color: #525354;
    border-radius: 10px;
    padding: 10px;
    width: 80px;
    fontSize: 25px;
    animation-name: fadeOut;
    animation-duration: 5s;
"""

#background of webpage
page_bg_img = f"""
<style>
[data-testid="stAppViewContainer"] > .main {{
background-image: url("https://raw.githubusercontent.com/nczub/AhR_Streamlit/main/AhR_background.svg");
background-size: cover;
background-position: top;
background-repeat: repeat;
background-attachment: local;
backgroud-color: #45745c;
}}
[data-testid="stHeader"] {{
background: rgba(0,0,0,0);

}}
</style>
"""
st.markdown(page_bg_img, unsafe_allow_html=True)

custom_css = """
<style>
    :root {
        font-size: 20px;
        text-align: justify;
    }
    .centered-image {
        display: flex;
        justify-content: center;
    }
    </style>

</style>
"""

st.markdown(custom_css, unsafe_allow_html=True)


footer_style = """
    position: fixed;
    left: 0;
    z-index: 3;
    bottom: 0;
    width: 100%;
    color: #525354;
    font-style: italic;
    text-align: center;
    padding: 10px;
    font-size: 15px;
"""
st.markdown(
    f'<div style="{footer_style}">Copyright (C) Natalia Łapińska 2024</div>',
    unsafe_allow_html=True
)


#loading models
#classification - threshold 1000 nM EC50
classification_model_path = 'mljar_AutoML_Compete_2024_08_29_23_38_44_1000' 

classification_model = AutoML(classification_model_path)
calc = Calculator(descriptors, ignore_3D=True)

#clearning SMILES input button
def clear_text():
    st.session_state["text"] = ""

st.title('Aryl hydrogen receptor QSAR model')
st.write("If you want to compare your experimental results with AhR QSAR model, provide test type and cell line. If your experiment was conducted according to different methodology, pick 'not selected'. Also for virtual screening choose 'not selected'.")
test_type_mapping = {
    "Not selected": 0,
    "Luminescence": 0,
    "Fluorescence_EROD": 1,
    "Fluorescence_NT": 2}
test_type = st.selectbox("Input test type:", list(test_type_mapping.keys()))
test_type_value = test_type_mapping[test_type]
cell_type_mapping = {
    "Not selected": 0,
    "HepG2": 0, 
    "HT29 ": 1, 
    "Huh-7": 2,
    "U937": 3, 
    "MCF7": 4,
    "HEK293": 5}
cell_line = st.selectbox("Input cell line:", list(cell_type_mapping.keys()))
cell_line_value = cell_type_mapping[cell_line]

smiles_input = st.text_input("Input SMILES", key="text")
col1, col2 = st.columns(2)
if smiles_input:
    try:
        molecule = Chem.MolFromSmiles(smiles_input)
        if molecule:
            img = Draw.MolToImage(molecule, size=(600, 600))
            with col1:
                st.image(img, caption='Chemical structure', use_column_width=True)
        else:
            pass
    except Exception as e:
        st.error(f"Wystąpił błąd: {str(e)}")
if smiles_input:
    try:
        molecule = Chem.MolFromSmiles(smiles_input)
        if molecule is not None:
            descriptors_value = calc.pandas([molecule])
            descriptors_value_df = pd.DataFrame(descriptors_value)
            for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                descriptors_value_df[column] = 0
            descriptors_value_df['TEST_TYPE'] = test_type_value
            descriptors_value_df['CELL_LINE'] = cell_line_value
            
            with col2:
                with st.spinner('Calculation in progress'):
                    prediction = classification_model.predict(descriptors_value_df)
                    st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                    if prediction == 0:
                        prediction_class = "inactive"
                    else:
                        prediction_class = "active"
                    st.write("Tested molecule towards AhR is ", f'<span style="color: #465e4b;">{prediction_class}</span>', unsafe_allow_html=True)
                    list_of_important_descriptors = ['Xch-5dv', 'AATSC1d', 'AATSC0p', 'MATS6v', 'AMID_X', 'PEOE_VSA10', 'GATS3s', 'MPC9', 'ATSC0c', 'AMID_O']
                    min_values = {'Xch-5dv': 0.0, 'AATSC1d': -0.1296285890880486, 'AATSC0p': 0.19011554804521, 'MATS6v': -0.3380648814398859,
 'AMID_X': 0.0, 'PEOE_VSA10': 0.0, 'GATS3s': 0.3614357159689514, 'MPC9': 0, 'ATSC0c': 0.0778819250415749, 'AMID_O': 0.0}
                    max_values = {'Xch-5dv': 0.8288787866419043, 'AATSC1d': 0.3693213296398892, 'AATSC0p': 0.4309255777726231, 'MATS6v': 1.2521308648311122,
 'AMID_X': 0.4319398031872989, 'PEOE_VSA10': 36.30663711745475, 'GATS3s': 1.6285321476764898, 'MPC9': 350, 'ATSC0c': 3.4150107712893245, 'AMID_O': 0.5190087882828514}
                    normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                    values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                    values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                    values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                    labels = normalized_descriptors_df[list_of_important_descriptors].columns
                    desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                    values_1 += values_1[:1]
                    values_2 += values_2[:1]
                    values_2 = [-1 if value < -1 else value for value in values_2]
                    values_2 = [1.5 if value > 1.5 else value for value in values_2]
                    values_3 += values_3[:1]
                    num_labels = len(labels)
                    angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                    angles += angles[:1]
                    fig = plt.figure(figsize=(7,7))
                    ax = fig.add_subplot(111, polar=True)
                    color_1 = '#A6A6A6'
                    color_2 = '#809f87'
                    ax.plot(angles, values_1, color=color_1, label="training set")
                    ax.fill(angles, values_1, alpha=0.25, color=color_1)
                    ax.plot(angles, values_3, color="white")
                    ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                    ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                    ax.fill(angles, values_2, alpha=0)
                    ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                    ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                    ax.legend()
                    plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                    plt.tight_layout()
                    st.pyplot(fig)
                    if desc_condition > 6:
                        st.write("The compound is under applicability domain, prediction is accurate!")
                    else:
                        st.write("The compound is not under applicability domain, prediction may be inaccurate.")
                    st.button("Clear SMILES", on_click=clear_text)
                def normalize_data(data, max_values):
                    return [x / max_val for x, max_val in zip(data, max_values)]
                def plot_radar_chart_one_condition(data1, data3, labels, max_values, label, title):
                    num_vars = len(labels)
                    data1_norm = normalize_data(data1, max_values)
                    data3_norm = normalize_data(data3, max_values)
                    data1_norm += data1_norm[:1]
                    data3_norm += data3_norm[:1]
                    angles = [n / float(num_vars) * 2 * np.pi for n in range(num_vars)]
                    angles += angles[:1]
                    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
                    line1, = ax.plot(angles, data1_norm, linewidth=1.5, linestyle='solid', label=label, color='#658c6d')
                    ax.fill(angles, data1_norm, '#658c6d', alpha=0.2)
                    line3, = ax.plot(angles, data3_norm, linewidth=1.5, linestyle='solid', label='Molecule', color='black')
                    ax.set_yticklabels([])
                    ax.set_xticks(angles[:-1])
                    ax.set_xticklabels(labels, color='black')
                    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
                    plt.title(title)
                    plt.text(0.1, -0.2, "Values were normalized", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 8, color="gray")
                    return fig
                st.write(' ')
            with st.expander("**See druglikeness and lead-like properties of a tested molecule**"):
                rotatable_bonds = Calculator(descriptors.RotatableBond.RotatableBondsCount())(Chem.MolFromSmiles(smiles_input))[descriptors.RotatableBond.RotatableBondsCount()]
                logP = Calculator(descriptors.SLogP.SLogP())(Chem.MolFromSmiles(smiles_input))[descriptors.SLogP.SLogP()]
                mw = Calculator(descriptors.Weight.Weight())(Chem.MolFromSmiles(smiles_input))[descriptors.Weight.Weight()]
                hdonors = Calculator(descriptors.HydrogenBond.HBondDonor())(Chem.MolFromSmiles(smiles_input))[descriptors.HydrogenBond.HBondDonor()]
                hacceptors = Calculator(descriptors.HydrogenBond.HBondAcceptor())(Chem.MolFromSmiles(smiles_input))[descriptors.HydrogenBond.HBondAcceptor()]
                labels_RO3 = ['Hydrogen donors', 'Hydrogen acceptors', 'Molecular weight', 'logP', 'Rotatable bonds']
                labels_RO5 = ['Hydrogen donors', 'Hydrogen acceptors', 'Molecular weight', 'logP']
                data1 = [5, 10, 500, 5]
                data2 = [3, 3, 300, 3, 3]
                data3 = [hdonors, hacceptors, mw, logP]
                data4 = [hdonors, hacceptors, mw, logP, rotatable_bonds]
                max_values_RO5 = [5, 10, 500, 5]
                max_values_RO3 = [3, 3, 300, 3, 3]
                fig_RO5 = plot_radar_chart_one_condition(data1, data3, labels_RO5, max_values_RO5, label = "RO5", title = "Lipiński's rule of five - RO5")
                fig_RO3 = plot_radar_chart_one_condition(data2, data4, labels_RO3, max_values_RO3, label = "RO3", title = "Rule of three - RO3")
                font_size_style = "font-size:13px; text-align:center;"
                col1, col2, col3, col4, col5 = st.columns([1,1,1,1,1])
                with col1:
                    st.write(f'<span style="{font_size_style}">Molecular weight: <span style="color: #5e7563;">{round(mw, 1)} Da</span></span>', unsafe_allow_html=True)
                with col2:
                    st.write(f'<span style="{font_size_style}">CLogP: <span style="color: #5e7563;">{round(logP, 2)}</span></span>', unsafe_allow_html=True)
                with col3:
                    st.write(f'<span style="{font_size_style}">Rotatable bonds: <span style="color: #5e7563;">{rotatable_bonds}</span></span>', unsafe_allow_html=True)
                with col4:
                    st.write(f'<span style="{font_size_style}">Hydrogen donors: <span style="color: #5e7563;">{hdonors}</span></span>', unsafe_allow_html=True)
                with col5:
                    st.write(f'<span style="{font_size_style}">Hydrogen acceptors: <span style="color: #5e7563;">{hacceptors}</span></span>', unsafe_allow_html=True)
                col1, col2 = st.columns(2)
                with col1:
                    st.pyplot(fig_RO5)
                    results_RO5 = [data3[i] <= data1[i] for i in range(len(data3))]
                    if all(results_RO5):
                        st.write("Lipinski's rule is fulfilled")
                    elif not any(results_RO5):
                        st.write("Lipinski's rule is not fulfilled")
                    else:
                        st.write("Lipinski's rule is fulfilled partially")
                with col2:
                    st.pyplot(fig_RO3)
                    results_RO3 = [data4[i] <= data2[i] for i in range(len(data4))]
                    if all(results_RO3):
                        st.write("The rule of three is fulfilled")
                    elif not any(results_RO3):
                        st.write("The rule of three is not fulfilled")
                    else:
                        st.write("The rule of three is fulfilled partially")

        else:
            st.write("Invalid SMILES")
    except Exception as e:
        st.write("Error:", e)
agree_draw_smiles = st.checkbox("Draw chemical structure")
if agree_draw_smiles:
    smile_code = st_ketcher()
    st.markdown(f"SMILES: {smile_code}")
        

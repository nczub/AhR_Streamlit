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
    text-align: right;
    padding: 10px;
    font-size: 20px;
"""
st.markdown(
    f'<div style="{footer_style}">Copyright (C) Natalia Łapińska 2024</div>',
    unsafe_allow_html=True
)


#loading models
#loading models
#classification - threshold 1000 nM EC50
# classification_model_path = 'mljar_AutoML_Compete_2024_08_29_23_38_44_1000' 
#classification - threshold 10000 nM EC50
classification_model_path = 'mljar_AutoML_Compete_2024_08_29_21_37_05_10000'
classification_model = AutoML(classification_model_path)
calc = Calculator(descriptors, ignore_3D=True)

#clearning SMILES input button
def clear_text():
    st.session_state["text"] = ""

st.title('AhR QSAR model - batch mode')
st.write('Predictions based on uploaded csv file')
uploaded_file = st.file_uploader('CSV file')

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


if uploaded_file is not None:
    data_file = pd.read_csv(uploaded_file)
    st.write(data_file.head())
    descriptors = []
    smiles = []
    for smi in data_file["smiles"]:
        mols = [Chem.MolFromSmiles(smi)]
        smiles.append(smi)
        descriptors.append(calc.pandas(mols))
    descriptors_df = pd.DataFrame(np.vstack(descriptors))
    descriptors_df.columns = descriptors[0].columns
    descriptors_df = descriptors_df.astype(float)
    descriptors_df = descriptors_df.fillna(0)
    descriptors_df['TEST_TYPE'] = test_type_value
    descriptors_df['CELL_LINE'] = cell_line_value
    list_of_important_descriptors = ['BCUTc-1l', 'MAXdO', 'MAXaaaC', 'PEOE_VSA9', 'SMR_VSA3', 'SssCH2', 
                                 'AATS4i', 'SpAbs_DzZ', 'AATS6dv', 'VSA_EState5']
    min_values = {'BCUTc-1l': -0.7373620145398293, 'MAXdO': 9.099537037037038, 'MAXaaaC': -0.1340347251287001,
                        'PEOE_VSA9': 0.0, 'SMR_VSA3': 0.0, 'SssCH2': -0.4661989795918362, 'AATS4i': 147.56632501478904,
                        'SpAbs_DzZ': 42.05895519992669, 'AATS6dv': 0.1538461538461538, 'VSA_EState5': -7.181078813682171}
    max_values = {'BCUTc-1l': -0.292146392415571, 'MAXdO': 14.629783178882205, 'MAXaaaC': 1.5381250000000002,
                        'PEOE_VSA9': 78.6625871223764, 'SMR_VSA3': 39.8909626482546, 'SssCH2': 23.225582666887828,
                        'AATS4i': 175.1107976137481, 'SpAbs_DzZ': 1265.278990098867, 'AATS6dv': 7.298507462686567,
                        'VSA_EState5': 8.521368790538302}
    normalized_descriptors_df = (descriptors_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
    values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    data = pd.DataFrame()
    if st.button("Calculate Predictions"):
        data = pd.DataFrame(columns = ["AhR_class"])
        for i in range(len(descriptors_df)):
            compound = pd.DataFrame(descriptors_df.iloc[i]).transpose()
            predictions = classification_model.predict(compound)
            if predictions == 0:
                data.loc[i, "AhR_class"] = "inactive"
            else:
                data.loc[i, "AhR_class"] = "active"
            values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
            values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
            values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
            if desc_condition > 6:
                data.loc[i, "applicability_domain"] = True  
            else:
                data.loc[i, "applicability_domain"] = False
        data_new = pd.concat([data_file, data], axis = 1)
        if data.empty:
            st.write("No predictions to download.")
        else:
            st.download_button(label="Download predictions as csv file", data=data_new.to_csv(index=False), file_name="AhR_predictions.csv") 

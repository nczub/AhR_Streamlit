#Natalia Lapinska - creation webapp started 26 August 2024


#load packages
import streamlit as st
import pandas as pd
import numpy as np
import html
from streamlit.components.v1 import html
from PIL import Image
import matplotlib.pyplot as plt

st.set_page_config(page_title="Aryl hydrocarbon receptor")


st.markdown("""
<style>
    [data-testid=stSidebar] {
        background-color: #a2c7ab  ;
    }
</style>
""", unsafe_allow_html=True)

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

st.title('Aryl hydrogen receptor')
st.write('The aryl hydrocarbon receptor (AhR) is part of a family of essential helix-loop-helix transcription factors. This receptor is critical in determining host physiology and various pathophysiologies, from inflammation and metabolism to cancer. AhR is a ligand-controlled receptor with complex activation pharmacology depending on the type and amount of ligand present.')
st.write('')
st.write('---')
st.subheader("Instruction")
st.write("With this application - QSAR page, it is possible to obtain information about the activity class of ligands towards the AhR receptor based on EC50 values. All you need is a SMILES of your molecule. An alternative way, you can draw a molecule.")
st.write('')
st.write('---')
st.subheader("References")
st.write(":sparkle: AhR database was curated by MSc Paulina Wojtyło.")
st.write(":sparkle: For model development Mljar tool was used. Plonska, A.; Plonski, P. MLJAR: State-of-the-Art Automated Machine Learning Framework for Tabular Data, version 0.10.3. https://github.com/mljar/mljar-supervised.")
st.write(":sparkle: Molecular representation used in the research was Mordred descriptors. Moriwaki H, Tian Y-S, Kawashita N, Takagi T (2018) Mordred: a molecular descriptor calculator. Journal of Cheminformatics 10:4. doi: 10.1186/s13321-018-0258-y")
st.write(":sparkle: The background image was download from Uniprot.")

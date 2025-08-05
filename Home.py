# Natalia Lapinska 

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.


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

# #background of webpage
# page_bg_img = f"""
# <style>
# [data-testid="stAppViewContainer"] > .main {{
# background-image: url("https:/raw.githubusercontent.com/nczub/AhR_Streamlit/AhR_background.svg");
# background-size: cover;
# background-position: top;
# background-repeat: repeat;
# background-attachment: local;
# backgroud-color: #45745c;
# }}
# [data-testid="stHeader"] {{
# background: rgba(0,0,0,0);

# }}
# </style>
# """
# st.markdown(page_bg_img, unsafe_allow_html=True)

page_bg_img = f"""
<style>
[data-testid="stAppViewContainer"] > .main {{
    background-image: url("https://raw.githubusercontent.com/nczub/AhR_Streamlit/main/AhR_background.svg");
    background-size: cover;
    background-position: top;
    background-repeat: repeat;
    background-attachment: local;
    # background-color: #45745c;
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
st.write("With this application - QSAR page, it is possible to obtain information about the activity class of ligands towards the AhR receptor based on EC50 and IC50 values. All you need is a SMILES of your molecule. An alternative way, you can draw a molecule. A threshold of active/inactive molecules was set at 1000 nM.")
st.write('')
st.write('---')
st.subheader("Publication")
st.write("A reference to an article that provides a detailed description of the database and the model used within it:")
st.write("*Initial Development of Automated Machine Learning-Assisted Prediction Tools for Aryl Hydrocarbon Receptor Activators*, Wojtyło Paulina, Natalia Łapińska, Lucia Bellagamba, Emidio Camaioni, Aleksander Mendyk, and Stefano Giovagnoli. 2024 Pharmaceutics, https://doi.org/10.3390/pharmaceutics16111456")

st.subheader("References")
st.write(":sparkle: AhR databases were curated by MSc Paulina Wojtyło.")
st.write(":sparkle: AhR QSAR models were developed by MSc Paulina Wojtyło with the support of Natalia Łapińska.")
st.write(":sparkle: For model development Mljar tool was used. Plonska, A.; Plonski, P. MLJAR: State-of-the-Art Automated Machine Learning Framework for Tabular Data, version 0.10.3. https://github.com/mljar/mljar-supervised.")
st.write(":sparkle: Molecular representation used in the research was *Mordred descriptors* and *Morgan Fingerprints*.")
st.write("Moriwaki H, Tian Y-S, Kawashita N, Takagi T (2018) Mordred: a molecular descriptor calculator. Journal of Cheminformatics 10:4. doi: 10.1186/s13321-018-0258-y")
st.write("RDKit blog Greg Landrum https://greglandrum.github.io/rdkit-blog/posts/2023-01-18-fingerprint-generator-tutorial.html")
st.write(":sparkle: The background image was download from Uniprot.")
st.subheader('License')
st.write('GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007')
st.write('Copyright (C) 2024 Natalia Łapińska')
st.write('This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. https://www.gnu.org/licenses/gpl-3.0.html')
st.write('DISCLAIMER OF LIABILITY')
st.write("The author of this software shall not be liable for any special, incidental, consequential, or indirect damages resulting from the use, misuse, or inability to use this software, including but not limited to, damages for loss of profits, business interruption, or loss of data. The software is provided 'as is' and the author make no warranties, either express or implied, regarding the software's fitness for a particular purpose or its accuracy and reliability.")


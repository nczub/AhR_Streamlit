import streamlit as st

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
    text-align: right;
    padding: 10px;
    font-size: 20px;
"""
st.markdown(
    f'<div style="{footer_style}">Copyright (C) Natalia Łapińska 2024</div>',
    unsafe_allow_html=True
)


st.title('Contact')
st.write("My name is **Natalia Łapińska** and I am the author of ***AhR app*** and also ***SerotoninAI*** (https://serotoninai.streamlit.app/).")
st.write('A few words about me, I am a Polish woman 🇵🇱 who fell in love with the world of artificial intelligence and the positive opportunities it gives us to improve all areas of research.')
st.write("Wanting to define myself, one term is not enough. Here's a list of the ones that best fit scientist :microscope:, pharmacist 💊, PhD 👩‍🎓, Python programmer :snake:, cheminformatics :female-technologist:, data scientist 📊 and AI specialist🔥.")
st.write("If you want to develop AI-based drug discovery together with me, I invite you to contact me!")
st.markdown('<i class="fab fa-linkedin"></i>  www.linkedin.com/in/natalia-czub', unsafe_allow_html=True)
st.markdown('<i class="fab fa-github-square"></i>  https://github.com/nczub', unsafe_allow_html=True)
st.markdown('<i class="fa-solid fa-message"></i>  lapinska.natalia@inbox.eu', unsafe_allow_html=True)

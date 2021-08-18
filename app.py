import streamlit as st
from multiapp import MultiApp
import code_KinomeTerrain
import process_new_data # import your app modules here

app = MultiApp()
st.set_page_config(
    page_title="KinomeTerrain App",
    page_icon=":shark:",
    layout="wide")
# Add all your application here
app.add_app("KinomeTerrain", code_KinomeTerrain.app)
app.add_app("Analysis your data !", process_new_data.app)
# The main app
app.run()


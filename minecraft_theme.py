"""
Minecraft-style theme for Streamlit.
Pure CSS - no images needed. Call apply_minecraft_theme() once, near the top of app.py.
To tweak colours, edit the values in the CSS below.
"""

import streamlit as st

MINECRAFT_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Press+Start+2P&family=VT323&display=swap');

/* ---------- Page background: dirt block pattern ---------- */
.stApp {
    background-color: #5a3d22;
    background-image:
        linear-gradient(45deg, rgba(0,0,0,.18) 25%, transparent 25%, transparent 75%, rgba(0,0,0,.18) 75%),
        linear-gradient(45deg, rgba(0,0,0,.18) 25%, transparent 25%, transparent 75%, rgba(0,0,0,.18) 75%);
    background-size: 48px 48px;
    background-position: 0 0, 24px 24px;
}

/* ---------- Grass top bar ---------- */
[data-testid="stHeader"] {
    background: #5d9b3a;
    border-bottom: 6px solid #3b6b22;
    box-shadow: 0 6px 0 #7a4f2a;
}

/* ---------- Main content = dark stone panel ---------- */
.block-container {
    background: rgba(32, 32, 32, 0.93);
    border: 4px solid #000;
    box-shadow: inset 4px 4px 0 #5c5c5c, inset -4px -4px 0 #1a1a1a;
    padding: 2.5rem 2.5rem 3rem 2.5rem !important;
    margin-top: 1.2rem;
    margin-bottom: 1.2rem;
}

/* ---------- Fonts ---------- */
.stApp p, .stApp label, .stApp li, .stApp input, .stApp textarea,
.stApp [data-testid="stMarkdownContainer"], .stApp [data-testid="stCaptionContainer"] {
    font-family: 'VT323', monospace !important;
    font-size: 1.3rem;
    color: #f2f2f2;
}
h1, h2, h3, h4, .main-title, .section-title {
    font-family: 'Press Start 2P', monospace !important;
    line-height: 1.6 !important;
    letter-spacing: 0;
}
h1 { font-size: 1.5rem !important; color: #fcfc54 !important; text-shadow: 3px 3px 0 #3f3f15; }
h2 { font-size: 1.15rem !important; color: #55ff55 !important; text-shadow: 3px 3px 0 #153f15; }
h3 { font-size: 0.95rem !important; color: #55ffff !important; text-shadow: 2px 2px 0 #153f3f; }
h4 { font-size: 0.85rem !important; color: #ffaa00 !important; text-shadow: 2px 2px 0 #3f2a00; }
.main-title { text-align: center; font-size: 1.6rem; color: #fcfc54; text-shadow: 4px 4px 0 #3f3f15; margin: 10px 0 18px 0; }
.mc-subtitle { text-align: center; font-family: 'VT323', monospace; font-size: 1.7rem; color: #ffaa00; text-shadow: 2px 2px 0 #3f2a00; }

/* ---------- Sidebar = stone ---------- */
[data-testid="stSidebar"] {
    background-color: #6b6b6b;
    background-image:
        linear-gradient(rgba(0,0,0,.12) 2px, transparent 2px),
        linear-gradient(90deg, rgba(0,0,0,.12) 2px, transparent 2px);
    background-size: 32px 32px;
    border-right: 6px solid #1a1a1a;
}
[data-testid="stSidebar"] * { text-shadow: 2px 2px 0 #3f3f3f; }
[data-testid="stSidebar"] h1 { font-size: 1.1rem !important; }

/* ---------- Buttons = Minecraft buttons ---------- */
.stButton > button, .stDownloadButton > button, [data-testid="stFormSubmitButton"] > button {
    background: #7d7d7d;
    color: #fff;
    border: 3px solid;
    border-color: #dcdcdc #373737 #373737 #dcdcdc;
    border-radius: 0;
    font-family: 'Press Start 2P', monospace;
    font-size: 0.7rem;
    padding: 0.8rem 1rem;
    text-shadow: 2px 2px 0 #3f3f3f;
    box-shadow: 0 0 0 3px #000;
    transition: none;
}
.stButton > button:hover, .stDownloadButton > button:hover {
    background: #7e8ad0;
    color: #ffffa0;
    border-color: #e6e9ff #3a3f73 #3a3f73 #e6e9ff;
}
.stButton > button:active, .stDownloadButton > button:active {
    border-color: #373737 #dcdcdc #dcdcdc #373737;
}
.stButton > button p, .stDownloadButton > button p {
    font-family: 'Press Start 2P', monospace !important;
    font-size: 0.7rem !important;
    color: #fff !important;
}

/* ---------- Inputs ---------- */
.stTextInput input, .stTextArea textarea, .stNumberInput input {
    background: #000 !important;
    color: #fff !important;
    border: 3px solid #a0a0a0 !important;
    border-radius: 0 !important;
}
[data-baseweb="select"] > div {
    background: #000 !important;
    border: 3px solid #a0a0a0 !important;
    border-radius: 0 !important;
}

/* ---------- Metrics = inventory slots ---------- */
[data-testid="stMetric"] {
    background: #8b8b8b;
    border: 4px solid;
    border-color: #373737 #ffffff #ffffff #373737;
    padding: 12px 16px;
}
[data-testid="stMetric"] * { color: #1a1a1a !important; text-shadow: none !important; }
[data-testid="stMetricValue"] { font-family: 'Press Start 2P', monospace !important; font-size: 1.1rem !important; }

/* ---------- Alerts = chat / achievement boxes ---------- */
[data-testid="stAlert"] {
    border-radius: 0 !important;
    border: 3px solid #000 !important;
    box-shadow: inset 3px 3px 0 rgba(255,255,255,.25), inset -3px -3px 0 rgba(0,0,0,.35);
}

/* ---------- Images, tables, misc ---------- */
[data-testid="stImage"] img { border: 4px solid #000; border-radius: 0; }
[data-testid="stDataFrame"] { border: 4px solid #000; }
hr { border: 0; border-top: 4px solid #000; box-shadow: 0 4px 0 #5c5c5c; }
[data-testid="stExpander"] { border: 3px solid #000; border-radius: 0; background: rgba(60,60,60,.6); }
[data-testid="stProgress"] > div > div > div > div { background-color: #55ff55; border-radius: 0; }
[data-baseweb="tab"] { font-family: 'Press Start 2P', monospace !important; font-size: 0.65rem; }
</style>
"""


def apply_minecraft_theme():
    """Inject the Minecraft CSS into the page."""
    st.markdown(MINECRAFT_CSS, unsafe_allow_html=True)

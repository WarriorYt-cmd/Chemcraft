"""
ChemBot - the AI helper inside ChemCraft Lab.

The API key is NEVER written in code.  It is read from Streamlit "Secrets"
(Streamlit Cloud -> your app -> Settings -> Secrets) or from environment
variables.  See secrets.toml.example for the exact text to paste.

Two providers are supported:

 1) Anthropic (Claude)       ANTHROPIC_API_KEY  (+ optional ANTHROPIC_MODEL)
 2) Any "OpenAI-compatible"  AI_API_KEY, AI_BASE_URL, AI_MODEL
    API, e.g. Google Gemini, Groq, OpenRouter - several have free tiers.
"""

import os
import requests
import streamlit as st

# A small, fast, low-cost model. Change it here if you prefer another one.
DEFAULT_ANTHROPIC_MODEL = "claude-haiku-5-5"

MAX_HISTORY_MESSAGES = 8      # how many past messages are sent each time
MAX_USER_CHARS = 1200         # longest question accepted
MAX_OUTPUT_TOKENS = 700       # keeps answers short and cheap
MAX_QUESTIONS_PER_SESSION = 30

SYSTEM_PROMPT = """You are ChemBot, the friendly assistant inside ChemCraft Lab, a Minecraft-themed \
cheminformatics learning app for undergraduate chemistry students.

You help with: SMILES notation, 2D/3D molecular structures, molecular descriptors (molecular weight, \
LogP, TPSA, hydrogen-bond donors and acceptors, rotatable bonds, rings, Fsp3, molar refractivity), \
Lipinski's Rule of Five, structure-property relationships, melting and boiling points, intermolecular \
forces, functional groups, and how to use the app's pages.

Rules:
- Keep answers clear and fairly short (under about 200 words unless the student asks for more).
- Explain ideas simply, with a small example when helpful. A light Minecraft-style touch is fine, but do not overdo it.
- If a molecule's data is supplied below, use it. Do not invent numbers. If you are unsure of an \
experimental value, say so and suggest PubChem.
- If asked something unrelated to chemistry, cheminformatics or this app, politely steer back.
- Do not give instructions for making dangerous substances (weapons, explosives, illegal drugs)."""


def _secret(name, default=None):
    """Read a value from Streamlit secrets, else from environment variables."""
    try:
        if name in st.secrets:
            return str(st.secrets[name])
    except Exception:
        pass
    return os.environ.get(name, default)


def get_ai_config():
    """Return a config dict, or None if no API key has been set up."""
    provider = (_secret("AI_PROVIDER") or "").strip().lower()

    if not provider:
        if _secret("ANTHROPIC_API_KEY"):
            provider = "anthropic"
        elif _secret("AI_API_KEY"):
            provider = "openai_compatible"

    if provider == "anthropic" and _secret("ANTHROPIC_API_KEY"):
        return {
            "provider": "anthropic",
            "key": _secret("ANTHROPIC_API_KEY"),
            "model": _secret("ANTHROPIC_MODEL", DEFAULT_ANTHROPIC_MODEL),
        }

    if provider in ("openai_compatible", "openai") and _secret("AI_API_KEY"):
        base = _secret("AI_BASE_URL")
        model = _secret("AI_MODEL")
        if base and model:
            return {
                "provider": "openai_compatible",
                "key": _secret("AI_API_KEY"),
                "base_url": base,
                "model": model,
            }

    return None


def _trim_history(messages):
    """Keep the last few messages and make sure the first one is from the user."""
    trimmed = messages[-MAX_HISTORY_MESSAGES:]
    while trimmed and trimmed[0]["role"] != "user":
        trimmed = trimmed[1:]
    return trimmed


def _friendly_error(status_code):
    if status_code in (401, 403):
        return "The AI key was rejected. Please check the API key in your app's Secrets."
    if status_code == 404:
        return "The AI model name was not found. Please check the model name in your Secrets."
    if status_code == 429:
        return "The AI is busy or the usage limit was reached. Please wait a moment and try again."
    return f"The AI service returned an error (code {status_code}). Please try again later."


def ask_ai(messages, molecule_context=""):
    """
    Send the chat to the AI.

    messages: list of {"role": "user"|"assistant", "content": str}
    Returns (ok, text).  When ok is False, text is a friendly error message.
    """
    config = get_ai_config()
    if config is None:
        return False, "ChemBot is not set up yet. Add an API key in Secrets (see the setup box)."

    history = _trim_history(messages)
    if not history:
        return False, "Please type a question first."

    system = SYSTEM_PROMPT
    if molecule_context:
        system += "\n\nMolecule the student is currently looking at:\n" + molecule_context

    try:
        if config["provider"] == "anthropic":
            response = requests.post(
                "https://api.anthropic.com/v1/messages",
                headers={
                    "x-api-key": config["key"],
                    "anthropic-version": "2023-06-01",
                    "content-type": "application/json",
                },
                json={
                    "model": config["model"],
                    "max_tokens": MAX_OUTPUT_TOKENS,
                    "system": system,
                    "messages": history,
                },
                timeout=60,
            )
            if response.status_code != 200:
                return False, _friendly_error(response.status_code)
            blocks = response.json().get("content", [])
            text = "".join(b.get("text", "") for b in blocks if b.get("type") == "text")

        else:
            response = requests.post(
                config["base_url"].rstrip("/") + "/chat/completions",
                headers={
                    "Authorization": f"Bearer {config['key']}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": config["model"],
                    "max_tokens": MAX_OUTPUT_TOKENS,
                    "messages": [{"role": "system", "content": system}] + history,
                },
                timeout=60,
            )
            if response.status_code != 200:
                return False, _friendly_error(response.status_code)
            text = response.json()["choices"][0]["message"]["content"]

        text = (text or "").strip()
        if not text:
            return False, "The AI sent an empty answer. Please try asking again."
        return True, text

    except requests.exceptions.Timeout:
        return False, "The AI took too long to answer. Please try again."
    except Exception:
        return False, "Could not reach the AI service. Please check your internet connection and try again."

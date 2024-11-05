import rdkit
from rdkit.Chem import Draw  # https://github.com/rdkit/rdkit/issues/2458
import base64
from io import BytesIO
import streamlit as st


@st.cache_data
def mol_to_image(smi: str) -> str:
    """Returns molecular image as data URI."""
    mol = rdkit.Chem.MolFromSmiles(smi)
    pil_image = Draw.MolToImage(mol, size=(200, 200))

    with BytesIO() as buffer:
        pil_image.save(buffer, "png")
        data = base64.encodebytes(buffer.getvalue()).decode("utf-8")

    return f"data:image/png;base64,{data}"

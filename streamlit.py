import streamlit as st
import pandas as pd
import numpy as np
from rdkit.Chem import DataStructs
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from fingerprints import (
    RDKIT_MACCS,
    PubChem,
    ECFP6,
    KlekotaRothFingerprint,
    CDK_PUBCHEM,
)

## Fingerprint classes to use
fingerprint_classes = {
    "RDKIT_MACCS": RDKIT_MACCS(),
    "RDKIT_PubChem": PubChem(),
    "ECFP6": ECFP6(),
    "Klekota_and_Roth": KlekotaRothFingerprint(),
    "CDK_PUBCHEM": CDK_PUBCHEM(),
}


# Set page title
st.set_page_config(page_title="Tanimoto Similarity")

st.title("Tanimoto Similarity Clustering of Chemical Compounds")

# Sidebar for uploading file
uploaded_file = st.sidebar.file_uploader("Choose a CSV file", type="csv")

# Placeholder for data
df = None

if uploaded_file:
    # Load the data
    df = pd.read_csv(uploaded_file)

    # Display data
    st.write("Uploaded Data")
    st.write(df.head())

    # Dynamically select SMILES and ID columns
    smiles_col = st.sidebar.selectbox("Select SMILES column", df.columns)
    id_col = st.sidebar.selectbox("Select ID column", df.columns)

    fp_algorithm = st.sidebar.radio(
        "Select fingerprint algorithm",
        options=[
            "RDKIT_PubChem",
            "RDKIT_MACCS",
            "ECFP6",
            "Klekota_and_Roth",
            "CDK_PUBCHEM",
        ],
    )

    # Similarity threshold slider
    threshold = st.sidebar.slider(
        "Select threshold dissimilarity value", 0.0, 1.0, 0.25
    )

    # Select the fingerprint calculation method based on user choice
    fingerprints = []
    fingerprint_instance = fingerprint_classes.get(fp_algorithm)

    for smiles in df[smiles_col]:
        try:
            fingerprint = fingerprint_instance.calculate_fingerprint(smiles)
            fingerprints.append(fingerprint)
        except Exception as e:
            st.error(f"Error calculating fingerprint for {smiles}: {e}")

    # Calculate Tanimoto similarity matrix
    num_mols = len(fingerprints)
    similarity_matrix = np.zeros((num_mols, num_mols))

    for i in range(num_mols):
        for j in range(num_mols):
            similarity_matrix[i, j] = DataStructs.TanimotoSimilarity(
                fingerprints[i], fingerprints[j]
            )

    # Convert similarity to distance matrix
    distance_matrix = 1 - similarity_matrix

    # Create linkage matrix for hierarchical clustering
    linkage_matrix = linkage(squareform(distance_matrix), method="average")

    # Plot dendrogram
    fig, ax = plt.subplots(figsize=(10, 6))
    dendrogram(linkage_matrix, labels=df[id_col].tolist(), ax=ax)

    # Mark the similarity threshold on the dendrogram
    threshold_value = 1 - threshold  # Convert similarity threshold to dissimilarity
    ax.hlines(threshold_value, *ax.get_xlim(), colors="red", linestyles="dashed")
    ax.text(
        0,
        threshold_value + 0.02,
        f"Threshold: {threshold:.2f}",
        color="red",
        fontsize=12,
    )

    plt.title("Hierarchical Clustering of Chemical Compounds")
    plt.xlabel("Compound ID")
    plt.ylabel("Dissimilarity Score")
    st.pyplot(fig)

    # Show distance matrix
    st.write("Distance Matrix")
    st.dataframe(pd.DataFrame(distance_matrix, columns=df[id_col], index=df[id_col]))

    # Threshold line on dendrogram
    st.write(f"Dissimilarity Threshold: {threshold}")

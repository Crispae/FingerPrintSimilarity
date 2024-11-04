# Chemical similarity using fingerprint

## Overview

This  application allows users to upload a CSV file containing chemical compound data and calculates Tanimoto similarity clustering among the compounds. The app utilizes various fingerprinting algorithms to compute similarities and visualizes the results in a hierarchical dendrogram.

## Features

- **Upload CSV File:** Users can upload a CSV file containing SMILES representations of chemical compounds.
- **Fingerprint Algorithms:** Choose from multiple fingerprinting algorithms (RDKIT MACCS, RDKIT PubChem, ECFP6, Klekota and Roth, CDK PubChem).
- **Dissimilarity Threshold:** Set a threshold value to visualize similarity on the dendrogram.
- **Dendrogram Visualization:** View hierarchical clustering results with a threshold line for easy interpretation.

## Requirements

To run this app

You can install the necessary Python libraries using pip. Here's an example of the requirements:

```bash
git clone https://github.com/Crispae/FingerPrintSimilarity
cd FingerPrintSimilarity
pip install -r requirements.txt
```

## Running with Docker

If you prefer to run the app in a Docker container, follow these steps:

Build the Docker image:

```bash
docker build -t fingerprint_similarity_app .
```
Run the Docker container:

```bash
docker run -d -p 8501:8501 fingerprint_similarity_app
```

Access the app: Open your web browser and go to http://localhost:8501 to view the app.

## Acknowledgments
- This app is developed based on https://gitlab.opencode.de/BfR/ShinyTanimotoGrouping.git
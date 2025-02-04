import os
import pandas as pd
import scanpy as sc

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.dirname(__file__))  # Gets py_files path
BASE_DIR = os.path.dirname(BASE_DIR)  # Moves one level up to RNASEQ/

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load gene expression matrix
def load_expression_matrix():
    print("Loading partial gene expression data...")

    # Construct file paths
    matrix_path = os.path.join(DATA_DIR, "matrix.csv")
    metadata_path = os.path.join(DATA_DIR, "metadata.csv")
    
    # Ensure files exist before proceeding
    if not os.path.exists(matrix_path):
        raise FileNotFoundError(f"Matrix file not found: {matrix_path}")
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    # Load data
    expression_matrix = pd.read_csv(matrix_path, index_col=0)
    metadata = pd.read_csv(metadata_path)

    # Convert to AnnData format
    print("Converting to AnnData format...")
    adata = sc.AnnData(expression_matrix)
    adata.obs = metadata.set_index("sample_name")

    # Save raw data
    raw_data_path = os.path.join(OUTPUT_DIR, "raw_data.h5ad")
    print(f"Saving raw data to {raw_data_path}...")
    adata.write_h5ad(raw_data_path)

    print("Data loading complete. Proceeding to Quality Control...")

    return adata

if __name__ == "__main__":
    load_expression_matrix()


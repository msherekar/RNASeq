import scanpy as sc
import os

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "outputs")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Perform clustering
def cluster_cells():
    print("Clustering cells...")

    # Construct input file path
    reduced_data_path = os.path.join(DATA_DIR, "reduced_data.h5ad")
    if not os.path.exists(reduced_data_path):
        raise FileNotFoundError(f"Reduced data file not found: {reduced_data_path}")

    # Load dimensionality-reduced data
    adata = sc.read_h5ad(reduced_data_path)

    # Compute clustering using Leiden algorithm
    print("Computing Leiden clustering...")
    sc.tl.leiden(adata, resolution=0.5)

    # Save clustering results
    clustered_data_path = os.path.join(OUTPUT_DIR, "clustered_data.h5ad")
    print(f"Saving clustered data to {clustered_data_path}...")
    adata.write_h5ad(clustered_data_path)

    print("Clustering complete. Proceeding to Differential Expression Analysis...")

if __name__ == "__main__":
    cluster_cells()

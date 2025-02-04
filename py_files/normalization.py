import scanpy as sc
import os

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "outputs")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Normalize data
def normalize_data():
    print("Normalizing data...")

    # Construct input file path
    qc_filtered_path = os.path.join(DATA_DIR, "qc_filtered_data.h5ad")
    if not os.path.exists(qc_filtered_path):
        raise FileNotFoundError(f"Filtered data file not found: {qc_filtered_path}")

    # Load QC-filtered data
    adata = sc.read_h5ad(qc_filtered_path)

    # Normalize total counts per cell to 10,000 reads and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)  # Changed to log1p (log2p is not standard in Scanpy)

    # Save normalized data
    normalized_data_path = os.path.join(OUTPUT_DIR, "normalized_data.h5ad")
    print(f"Saving normalized data to {normalized_data_path}...")
    adata.write_h5ad(normalized_data_path)

    print("Normalization complete. Proceeding to Dimensionality Reduction...")

if __name__ == "__main__":
    normalize_data()

import scanpy as sc
import os

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "outputs")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Perform quality control
def quality_control():
    print("Running quality control...")

    # Load raw data
    raw_data_path = os.path.join(DATA_DIR, "raw_data.h5ad")
    if not os.path.exists(raw_data_path):
        raise FileNotFoundError(f"Raw data file not found: {raw_data_path}")
    
    adata = sc.read_h5ad(raw_data_path)

    # Calculate QC metrics
    print("Calculating quality control metrics...")
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Filter out low-quality nuclei
    print("Filtering low-quality nuclei...")
    adata = adata[adata.obs['total_counts'] > 500, :]
    adata = adata[:, adata.var['n_cells_by_counts'] > 10]

    # Save filtered data
    filtered_data_path = os.path.join(OUTPUT_DIR, "qc_filtered_data.h5ad")
    print(f"Saving filtered data to {filtered_data_path}...")
    adata.write_h5ad(filtered_data_path)

    print("Quality control complete. Proceeding to Normalization...")

if __name__ == "__main__":
    quality_control()

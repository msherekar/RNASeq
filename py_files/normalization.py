# RNA-Seq Pipeline: Normalization (normalization.py)

import scanpy as sc
import os

# Define file paths
data_dir = "outputs/"
output_dir = "outputs/"
os.makedirs(output_dir, exist_ok=True)

# Normalize data
def normalize_data():
    print("Normalizing data...")
    adata = sc.read_h5ad(os.path.join(data_dir, "qc_filtered_data.h5ad"))
    
    # Normalize total counts per cell to 10,000 reads and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Save normalized data
    adata.write_h5ad(os.path.join(output_dir, "normalized_data.h5ad"))
    print("Normalization complete. Proceeding to Dimensionality Reduction...")

if __name__ == "__main__":
    normalize_data()

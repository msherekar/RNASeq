# RNA-Seq Pipeline: Quality Control (quality_control.py)

import scanpy as sc
import os

# Define file paths
data_dir = "outputs/"
output_dir = "outputs/"
os.makedirs(output_dir, exist_ok=True)

# Perform quality control
def quality_control():
    print("Running quality control...")
    adata = sc.read_h5ad(os.path.join(data_dir, "raw_data.h5ad"))
    
    # Calculate QC metrics
    print("Calculating quality control metrics...")
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    # Filter out low-quality nuclei
    print("Filtering low-quality nuclei...")
    adata = adata[adata.obs['total_counts'] > 500, :]
    adata = adata[:, adata.var['n_cells_by_counts'] > 10]
    
    # Save filtered data
    print("Saving filtered data...")
    adata.write_h5ad(os.path.join(output_dir, "qc_filtered_data.h5ad"))
    print("Quality control complete. Proceeding to Normalization...")

if __name__ == "__main__":
    quality_control()

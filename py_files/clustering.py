# RNA-Seq Pipeline: Clustering (clustering.py)

import scanpy as sc
import os

# Define file paths
data_dir = "outputs/"
output_dir = "outputs/"
os.makedirs(output_dir, exist_ok=True)

# Perform clustering
def cluster_cells():
    print("Clustering cells...")
    adata = sc.read_h5ad(os.path.join(data_dir, "reduced_data.h5ad"))
    
    # Compute clustering using Leiden algorithm
    sc.tl.leiden(adata, resolution=0.5)
    
    # Save clustering results
    adata.write_h5ad(os.path.join(output_dir, "clustered_data.h5ad"))
    print("Clustering complete. Proceeding to Differential Expression Analysis...")

if __name__ == "__main__":
    cluster_cells()

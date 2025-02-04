# RNA-Seq Pipeline: Dimensionality Reduction (dimensionality_reduction.py)

import scanpy as sc
import os

# Define file paths
data_dir = "outputs/"
output_dir = "outputs/"
os.makedirs(output_dir, exist_ok=True)

# Perform dimensionality reduction
def reduce_dimensions():
    print("Performing PCA and UMAP...")
    adata = sc.read_h5ad(os.path.join(data_dir, "normalized_data.h5ad"))
    
    # Compute PCA
    sc.pp.pca(adata, n_comps=50)
    
    # Compute UMAP and t-SNE
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    
    # Save results
    adata.write_h5ad(os.path.join(output_dir, "reduced_data.h5ad"))
    print("Dimensionality Reduction complete. Proceeding to Clustering...")

if __name__ == "__main__":
    reduce_dimensions()

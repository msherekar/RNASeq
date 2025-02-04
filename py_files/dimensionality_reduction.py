import scanpy as sc
import os

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "outputs")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Perform dimensionality reduction
def reduce_dimensions():
    print("Performing PCA and UMAP...")

    # Construct input file path
    normalized_data_path = os.path.join(DATA_DIR, "normalized_data.h5ad")
    if not os.path.exists(normalized_data_path):
        raise FileNotFoundError(f"Normalized data file not found: {normalized_data_path}")

    # Load normalized data
    adata = sc.read_h5ad(normalized_data_path)

    # Compute PCA
    print("Computing PCA...")
    sc.pp.pca(adata, n_comps=50)

    # Compute UMAP and t-SNE
    print("Computing nearest neighbors for UMAP and t-SNE...")
    sc.pp.neighbors(adata)

    print("Computing UMAP...")
    sc.tl.umap(adata)

    print("Computing t-SNE...")
    sc.tl.tsne(adata)

    # Save results
    reduced_data_path = os.path.join(OUTPUT_DIR, "reduced_data.h5ad")
    print(f"Saving reduced data to {reduced_data_path}...")
    adata.write_h5ad(reduced_data_path)

    print("Dimensionality Reduction complete. Proceeding to Clustering...")

if __name__ == "__main__":
    reduce_dimensions()

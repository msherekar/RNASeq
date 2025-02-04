import scanpy as sc
import os

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "outputs")

# Ensure the necessary input file exists
clustered_data_path = os.path.join(DATA_DIR, "clustered_data.h5ad")
if not os.path.exists(clustered_data_path):
    raise FileNotFoundError(f"❌ Clustered data file not found: {clustered_data_path}")

# Perform differential expression analysis
def differential_expression():
    print("🔹 Performing differential expression analysis...")

    # Load clustered data
    adata = sc.read_h5ad(clustered_data_path)

    # Check if Leiden clustering exists
    if "leiden" not in adata.obs.columns:
        raise ValueError("❌ Leiden clustering is missing. Ensure `clustering.py` ran successfully.")

    # Identify marker genes for each cluster
    print("📊 Identifying marker genes using Wilcoxon rank-sum test...")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    # Verify if results exist
    if "rank_genes_groups" not in adata.uns:
        raise ValueError("❌ `rank_genes_groups` not found in AnnData object. Differential expression analysis may have failed.")

    # Save results
    diff_exp_data_path = os.path.join(DATA_DIR, "differential_expression.h5ad")
    print(f"📁 Saving differential expression results to {diff_exp_data_path}...")
    adata.write_h5ad(diff_exp_data_path)

    print("✅ Differential expression analysis complete. Proceeding to Pathway Analysis...")

if __name__ == "__main__":
    differential_expression()

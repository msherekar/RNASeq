import scanpy as sc
import os
import matplotlib.pyplot as plt

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "outputs")

# Ensure the necessary input file exists
diff_exp_data_path = os.path.join(DATA_DIR, "differential_expression.h5ad")
if not os.path.exists(diff_exp_data_path):
    raise FileNotFoundError(f"❌ Differential expression data file not found: {diff_exp_data_path}")

def generate_plots():
    print("🔹 Generating visualization plots...")

    # Load differential expression data
    adata = sc.read_h5ad(diff_exp_data_path)

    # Check if `rank_genes_groups` exists
    if "rank_genes_groups" not in adata.uns:
        raise ValueError("❌ `rank_genes_groups` not found in AnnData object. Ensure `differential_expression.py` ran successfully.")

    # Generate heatmap for top marker genes
    print("📊 Generating heatmap for top marker genes...")
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, show=False)
    heatmap_path = os.path.join(DATA_DIR, "marker_genes_heatmap.png")
    plt.savefig(heatmap_path)
    plt.close()  # Prevents execution from stopping
    print(f"📁 Heatmap saved to {heatmap_path}")

    # Generate dot plot
    print("📊 Generating dot plot for marker genes...")
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=10, show=False)
    dotplot_path = os.path.join(DATA_DIR, "marker_genes_dotplot.png")
    plt.savefig(dotplot_path)
    plt.close()
    print(f"📁 Dot plot saved to {dotplot_path}")

    # Generate violin plot
    print("📊 Generating violin plot for marker genes...")
    sc.pl.rank_genes_groups_violin(adata, n_genes=10, show=False)
    violin_path = os.path.join(DATA_DIR, "marker_genes_violin.png")
    plt.savefig(violin_path)
    plt.close()
    print(f"📁 Violin plot saved to {violin_path}")

    print("✅ Visualization plots generated successfully.")

if __name__ == "__main__":
    generate_plots()

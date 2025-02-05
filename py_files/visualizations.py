import scanpy as sc
import os
import matplotlib.pyplot as plt

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "outputs")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

# Ensure results directory exists
os.makedirs(RESULTS_DIR, exist_ok=True)

# Ensure the necessary input files exist
clustered_data_path = os.path.join(DATA_DIR, "clustered_data.h5ad")
if not os.path.exists(clustered_data_path):
    raise FileNotFoundError(f"❌ Clustered data file not found: {clustered_data_path}")

diff_exp_data_path = os.path.join(DATA_DIR, "differential_expression.h5ad")
if not os.path.exists(diff_exp_data_path):
    raise FileNotFoundError(f"❌ Differential expression data file not found: {diff_exp_data_path}")

def generate_plots():
    print("🔹 Generating visualization plots...")

    # Load clustered data
    cdata = sc.read_h5ad(clustered_data_path)

    # Check if Leiden clustering exists
    if "leiden" not in cdata.obs.columns:
        raise ValueError("❌ Leiden clustering is missing. Ensure `clustering.py` ran successfully.")

    # UMAP Plot
    print("📊 Generating UMAP plot for clustered cells...")
    sc.pl.umap(cdata, color=['leiden'], show=False)
    umap_path = os.path.join(RESULTS_DIR, "umap_clusters.png")
    plt.savefig(umap_path)
    plt.close()
    print(f"📁 UMAP plot saved to {umap_path}")

    # t-SNE Plot
    print("📊 Generating t-SNE plot for clustered cells...")
    sc.pl.tsne(cdata, color=['leiden'], show=False)
    tsne_path = os.path.join(RESULTS_DIR, "tsne_clusters.png")
    plt.savefig(tsne_path)
    plt.close()
    print(f"📁 t-SNE plot saved to {tsne_path}")

    # PCA Plot
    print("📊 Generating PCA plot for clustered cells...")
    sc.pl.pca(cdata, color=['leiden'], show=False)
    pca_path = os.path.join(RESULTS_DIR, "pca_clusters.png")
    plt.savefig(pca_path)
    plt.close()
    print(f"📁 PCA plot saved to {pca_path}")

    # Violin Plots for QC metrics
    print("📊 Generating violin plots for QC metrics...")
    sc.pl.violin(cdata, ['total_counts', 'n_genes_by_counts'], groupby='leiden', show=False)
    violin_qc_path = os.path.join(RESULTS_DIR, "violin_qc.png")
    plt.savefig(violin_qc_path)
    plt.close()
    print(f"📁 Violin QC plot saved to {violin_qc_path}")

    
    #Load differential expression data
    adata = sc.read_h5ad(diff_exp_data_path)

    # Check if `rank_genes_groups` exists
    if "rank_genes_groups" not in adata.uns:
        raise ValueError("❌ `rank_genes_groups` not found in AnnData object. Ensure `differential_expression.py` ran successfully.")

    # Generate heatmap for top marker genes
    print("📊 Generating heatmap for top marker genes...")
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, show=False)
    heatmap_path = os.path.join(RESULTS_DIR, "marker_genes_heatmap.png")
    plt.savefig(heatmap_path)
    plt.close()  # Prevents execution from stopping
    print(f"📁 Heatmap saved to {heatmap_path}")

    # Generate dot plot
    print("📊 Generating dot plot for marker genes...")
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=10, show=False)
    dotplot_path = os.path.join(RESULTS_DIR, "marker_genes_dotplot.png")
    plt.savefig(dotplot_path)
    plt.close()
    print(f"📁 Dot plot saved to {dotplot_path}")

    # Generate violin plot
    print("📊 Generating violin plot for marker genes...")
    sc.pl.rank_genes_groups_violin(adata, n_genes=10, show=False)
    violin_path = os.path.join(RESULTS_DIR, "marker_genes_violin.png")
    plt.savefig(violin_path)
    plt.close()
    print(f"📁 Violin plot saved to {violin_path}")

    print("✅ Visualization plots generated successfully.")

if __name__ == "__main__":
    generate_plots()

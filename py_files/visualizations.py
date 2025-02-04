# RNA-Seq Pipeline: Visualizations (visualizations.py)

import scanpy as sc
import os
import matplotlib.pyplot as plt

# Define file paths
data_dir = "outputs/"
output_dir = "outputs/"
plot_dir = "outputs/plots/"
os.makedirs(plot_dir, exist_ok=True)

# Generate visualizations
def generate_plots():
    print("Generating visualizations...")
    adata = sc.read_h5ad(os.path.join(data_dir, "clustered_data.h5ad"))
    
    # UMAP Plot
    sc.pl.umap(adata, color=['leiden'], save="_clusters.png")
    
    # t-SNE Plot
    sc.pl.tsne(adata, color=['leiden'], save="_tsne_clusters.png")
    
    # Violin Plots for QC metrics
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts'], groupby='leiden', save="_qc_metrics.png")
    
    # Heatmap of top marker genes
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, save="_marker_genes.png")
    
    print("Plots saved in outputs/plots/")

if __name__ == "__main__":
    generate_plots()

# RNA-Seq Pipeline: Differential Expression Analysis (differential_expression.py)

import scanpy as sc
import os

# Define file paths
data_dir = "outputs/"
output_dir = "outputs/"
os.makedirs(output_dir, exist_ok=True)

# Perform differential expression analysis
def differential_expression():
    print("Performing differential expression analysis...")
    adata = sc.read_h5ad(os.path.join(data_dir, "clustered_data.h5ad"))
    
    # Identify marker genes for each cluster
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
    
    # Save results
    adata.write_h5ad(os.path.join(output_dir, "differential_expression.h5ad"))
    print("Differential expression analysis complete. Proceeding to Pathway Analysis...")

if __name__ == "__main__":
    differential_expression()


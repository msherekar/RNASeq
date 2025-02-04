# RNA-Seq Pipeline: Pathway Analysis (pathway_analysis.py)

import gseapy as gp
import scanpy as sc
import os
import pandas as pd

# Define file paths
data_dir = "outputs/"
output_dir = "outputs/"
os.makedirs(output_dir, exist_ok=True)

# Perform pathway enrichment analysis
def pathway_analysis():
    print("Performing pathway enrichment analysis...")
    adata = sc.read_h5ad(os.path.join(data_dir, "differential_expression.h5ad"))
    
    # Extract differentially expressed genes
    result = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(100)  # Take top 100 genes
    genes = result.values.flatten()
    
    # Run Gene Set Enrichment Analysis (GSEA) using KEGG pathways
    gsea_results = gp.enrichr(gene_list=genes.tolist(), 
                              gene_sets='KEGG_2019_Human', 
                              organism='Human', 
                              outdir=output_dir)
    
    # Save results
    gsea_results.res2d.to_csv(os.path.join(output_dir, "pathway_analysis_results.csv"))
    print("Pathway analysis complete. Pipeline execution finished!")

if __name__ == "__main__":
    pathway_analysis()

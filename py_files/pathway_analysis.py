import gseapy as gp
import scanpy as sc
import os
import pandas as pd

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "outputs")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Perform pathway enrichment analysis
def pathway_analysis():
    print("Performing pathway enrichment analysis...")

    # Construct input file path
    diff_exp_data_path = os.path.join(DATA_DIR, "differential_expression.h5ad")
    if not os.path.exists(diff_exp_data_path):
        raise FileNotFoundError(f"Differential expression data file not found: {diff_exp_data_path}")

    # Load differential expression data
    adata = sc.read_h5ad(diff_exp_data_path)

    # Extract differentially expressed genes
    print("Extracting top 100 differentially expressed genes...")
    if 'rank_genes_groups' not in adata.uns:
        raise ValueError("rank_genes_groups not found in AnnData object. Ensure differential expression analysis was run successfully.")

    result = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(100)  # Take top 100 genes
    genes = result.values.flatten()

    if genes.size == 0:
        raise ValueError("No genes extracted for pathway analysis. Check differential expression results.")

    # Run Gene Set Enrichment Analysis (GSEA) using KEGG pathways
    print("Running Gene Set Enrichment Analysis (GSEA) using KEGG 2019 Human...")
    gsea_results = gp.enrichr(
        gene_list=genes.tolist(),
        gene_sets="KEGG_2019_Human",
        organism="Human",
        outdir=OUTPUT_DIR
    )

    # Save results
    pathway_results_path = os.path.join(OUTPUT_DIR, "pathway_analysis_results.csv")
    print(f"Saving pathway analysis results to {pathway_results_path}...")
    gsea_results.res2d.to_csv(pathway_results_path)

    print("Pathway analysis complete. Pipeline execution finished!")

if __name__ == "__main__":
    pathway_analysis()

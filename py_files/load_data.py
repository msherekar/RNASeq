# RNA-Seq Pipeline: Load Data (load_data.py)

import os
import scanpy as sc
import pandas as pd
import numpy as np
import dask.dataframe as dd


# Define file paths (to be set in config.yaml)
data_dir = "/Users/mukulsherekar/pythonProject/RNASeq/"
output_dir = "/Users/mukulsherekar/pythonProject/RNASeq/outputs/"
os.makedirs(output_dir, exist_ok=True)

# Load gene expression matrix
def load_expression_matrix():
    print("Loading partial gene expression data using ...")
    
    # Load data with Dask (parallelized for large files)
    # expression_matrix = dd.read_csv(os.path.join(data_dir, "matrix.csv"), assume_missing=True, sample=1000000).compute()
    expression_matrix = pd.read_csv(os.path.join(data_dir, "matrix.csv"),index_col=0)
    
    # Load metadata using Pandas
    print("Loading metadata...")
    metadata = pd.read_csv(os.path.join(data_dir, "metadata.csv"))
    
    
    # Convert to AnnData format
    print("Converting to AnnData format...")
    adata = sc.AnnData(expression_matrix)
    adata.obs = metadata.set_index("sample_name")
    
    # Save raw data
    print("Saving raw data...")
    adata.write_h5ad(os.path.join(output_dir, "raw_data.h5ad"))
    print("Data loading complete. Proceeding to Quality Control...")
    
    return adata

if __name__ == "__main__":
    load_expression_matrix()

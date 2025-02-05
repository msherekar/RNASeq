* A toy pipeline to process RNA Seq data (see readme_data.txt for details) and run a toy machine learning pipellines.
* In the data file (matrix.csv), one row is for every cell and one column for every gene sequenced. 
* The values of the matrix represent total reads (introns + exons) for that gene (column) for that cell (row).
* Processing: Load->QC->Normalization->Dimensionality Reduction->Clustering->Differential Expression Analysis -> Pathway Analysis - > Visualizations
* Clustering was done using Leiden Algorithm. Each clustern represent groups of cells (or nuclei) with similar gene expression profiles, identified using graph-based clustering.
* The goal of ML models is to predict cluster label so when new data comes in, models can predict cluster type and hence type of cells.

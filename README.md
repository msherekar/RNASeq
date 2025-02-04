* A toy pipeline to process RNA Seq data (see readme_data.txt for details) and run a toy ML classification model.
* In the data file (matrix.csv), one row is for every cell and one column for every gene sequenced. 
* The values of the matrix represent total reads (introns + exons) for that gene (column) for that cell (row).
* Processing=> Load->QC->Normalization->Dimensionality Reduction->Clustering->Differential Expression Analysis -> Pathway Analysis - > Visualizations -> ML

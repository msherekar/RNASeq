import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report, confusion_matrix

# Step 1: Load normalized RNA-seq data
adata = sc.read_h5ad("/Users/mukulsherekar/pythonProject/RNASeq/outputs/normalized_data.h5ad")

# Inspect the dataset: print the shape and a glimpse of gene (variable) and sample (observation) names.
print("Data shape (samples x genes):", adata.shape)
print("First 5 gene names:", adata.var_names[:5])
print("First 5 sample names:", adata.obs_names[:5])

# Step 2: Feature Selection using Highly Variable Genes

# Assuming 'adata' has been loaded from normalized_data.h5ad in Step 1

# Compute highly variable genes using Scanpy's built-in function.
# The 'seurat' flavor is a popular choice; n_top_genes selects the top 2000 variable genes.
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)

# Inspect the number of genes identified as highly variable
num_hv_genes = adata.var['highly_variable'].sum()
print("Number of highly variable genes identified:", num_hv_genes)

# Subset the AnnData object to include only the highly variable genes
adata_selected = adata[:, adata.var['highly_variable']].copy()

print("Shape of the data after selecting highly variable genes (samples x genes):", adata_selected.shape)


# Step 3: Data Splitting, Scaling, and Model Training using cell_type_designation_label

# Extract the feature matrix (X) from the AnnData object
X = adata_selected.X
# If the feature matrix is in a sparse format, convert it to a dense array:
if hasattr(X, "toarray"):
    X = X.toarray()

# Inspect available columns to confirm the label column name
print("Available columns for labels:", adata_selected.obs.columns.tolist())

# Extract the labels using 'cell_type_designation_label'
y = adata_selected.obs['cell_type_designation_label'].values

# Print unique labels to understand the distribution
unique_labels = np.unique(y)
print("Unique labels in cell_type_designation_label:", unique_labels)

# Split the data into training and testing sets (80% train, 20% test), preserving label distribution
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

print("Training set shape:", X_train.shape)
print("Testing set shape:", X_test.shape)

# Standardize the features using StandardScaler
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train a logistic regression classifier as an example ML model
clf = LogisticRegression(max_iter=1000, random_state=42)
clf.fit(X_train_scaled, y_train)

# Predict on the test set and evaluate performance
y_pred = clf.predict(X_test_scaled)
accuracy = accuracy_score(y_test, y_pred)
print("Test set accuracy:", accuracy)

# Step 4: Model Evaluation and Interpretation

# Generate and print a classification report with precision, recall, and f1-score for each cell type
print("Classification Report:")
report = classification_report(y_test, y_pred)
print(report)

# Compute the confusion matrix to see the detailed breakdown of predictions vs. true labels
cm = confusion_matrix(y_test, y_pred)

# Because there may be many classes, the confusion matrix might be crowded.
# Optionally, you might consider plotting a confusion matrix for the top N most frequent classes.
# For now, we plot the full confusion matrix.
plt.figure(figsize=(12, 10))
sns.heatmap(cm, cmap="viridis", cbar=True,
            xticklabels=np.unique(y), yticklabels=np.unique(y))
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Confusion Matrix")
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()

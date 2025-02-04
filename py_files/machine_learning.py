import scanpy as sc
import numpy as np
import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "outputs")

def classification_analysis():
    print("🔹 Running machine learning classification...")

    # Ensure the necessary input file exists
    normalized_data_path = os.path.join(DATA_DIR, "normalized_data.h5ad")
    if not os.path.exists(normalized_data_path):
        raise FileNotFoundError(f"❌ Normalized data file not found: {normalized_data_path}")

    # Step 1: Load normalized RNA-seq data
    print("📥 Loading normalized RNA-seq data...")
    adata = sc.read_h5ad(normalized_data_path)

    # Step 2: Feature Selection using Highly Variable Genes
    print("📊 Selecting highly variable genes...")
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

    # Check if 'cell_type_designation_label' is available in metadata
    if "cell_type_designation_label" not in adata.obs.columns:
        raise ValueError("❌ Column 'cell_type_designation_label' not found in AnnData object. Check metadata.")

    # Extract highly variable genes
    adata_selected = adata[:, adata.var["highly_variable"]].copy()
    print(f"✅ Selected {adata_selected.shape[1]} highly variable genes.")

    # Step 3: Prepare Features and Labels
    print("🧬 Extracting feature matrix and labels...")
    X = adata_selected.X
    if hasattr(X, "toarray"):  # Convert sparse to dense if needed
        X = X.toarray()

    y = adata_selected.obs["cell_type_designation_label"].values

    # Print unique labels to understand distribution
    unique_labels = np.unique(y)
    print(f"🔬 Unique labels: {unique_labels}")

    print("🔍 Class distribution before train-test split:", Counter(y))

    # Replace rare classes with "Other"
    min_class_size = 2  # Classes must have at least 2 samples
    class_counts = Counter(y)

    # Modify labels before splitting
    y_merged = np.array(["Other" if class_counts[label] < min_class_size else label for label in y])

    print("🔍 Class distribution after merging rare classes:", Counter(y_merged))

    # Split the data into training and testing sets (80% train, 20% test)
    print("📌 Splitting data into train and test sets...")
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_merged, test_size=0.2, random_state=42, stratify=y_merged
    )

    # Step 4: Scale Features
    print("📏 Standardizing features...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Step 5: Train a Machine Learning Model (Logistic Regression)
    print("🤖 Training Logistic Regression Model...")
    clf = LogisticRegression(max_iter=1000, random_state=42)
    clf.fit(X_train_scaled, y_train)

    # Step 6: Evaluate Model Performance
    print("📊 Evaluating model performance...")
    y_pred = clf.predict(X_test_scaled)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"✅ Test set accuracy: {accuracy:.4f}")

    # Generate and print classification report
    report = classification_report(y_test, y_pred)
    print("\n🔹 Classification Report:\n", report)

    # Step 7: Compute Confusion Matrix
    print("📌 Computing confusion matrix...")
    cm = confusion_matrix(y_test, y_pred)

    # Save classification results to a text file
    results_file = os.path.join(DATA_DIR, "classification_results.txt")
    with open(results_file, "w") as f:
        f.write("📊 Machine Learning Classification Results\n")
        f.write("=" * 50 + "\n")
        f.write(f"✅ Test set accuracy: {accuracy:.4f}\n\n")
        f.write("🔹 Classification Report:\n")
        f.write(report + "\n")
        f.write("🔹 Confusion Matrix:\n")
        f.write(np.array2string(cm, separator=', ') + "\n")

    print(f"📁 Classification results saved to {results_file}")
    # Updated labels for confusion matrix
    merged_unique_labels = np.unique(y_merged)

    # Plot confusion matrix
    plt.figure(figsize=(12, 10))
    sns.heatmap(cm, cmap="viridis", cbar=True,
                xticklabels=merged_unique_labels, yticklabels=merged_unique_labels, annot=True, fmt="d")
    plt.xlabel("Predicted Label")
    plt.ylabel("True Label")
    plt.title("Confusion Matrix")
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()

    # Save the confusion matrix plot
    conf_matrix_path = os.path.join(DATA_DIR, "confusion_matrix.png")
    plt.savefig(conf_matrix_path)
    print(f"📁 Confusion matrix saved to {conf_matrix_path}")

    plt.close()
    print("🎯 Machine learning classification completed!")

if __name__ == "__main__":
    classification_analysis()

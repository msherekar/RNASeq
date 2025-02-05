import scanpy as sc
import pandas as pd
import numpy as np
import os
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Define paths relative to RNASEQ/
DATA_DIR = os.path.join(BASE_DIR, "outputs")

def supervised_learning():
    print("🔹 Performing supervised learning...")

    # Ensure the necessary input file exists
    clustered_data_path = os.path.join(DATA_DIR, "clustered_data.h5ad")
    if not os.path.exists(clustered_data_path):
        raise FileNotFoundError(f"❌ Clustered data file not found: {clustered_data_path}")
    
    # Step 1: Load clustered RNA-seq data
    print("📥 Loading clustered RNA-seq data...")
    adata = sc.read_h5ad(clustered_data_path)

    # Step 2: Prepare Features and Labels
    print("🧬 Extracting feature matrix and labels...")
    X = adata.X  # Gene expression matrix
    y = adata.obs["leiden"].astype(str)  # Ensure clusters are categorical

    # Step 3: Train-Test Split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
    # Encode labels
    label_encoder = LabelEncoder()
    y_train = label_encoder.fit_transform(y_train)
    y_test = label_encoder.transform(y_test)
    
    # Step 4: Train Random Forest Classifier
    print("🌲 Training Random Forest Classifier...")
    rf_model = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1, bootstrap=False)
    rf_model.fit(X_train, y_train)

    # Step 5: Train XGBoost Classifier
    print("🚀 Training XGBoost Classifier...")
    xgb_model = XGBClassifier(n_estimators=100, learning_rate=0.1, use_label_encoder=False, eval_metric="mlogloss", random_state=42, tree_method="hist", device="cuda")
    xgb_model.fit(X_train, y_train)

    # Step 6: Evaluate Model Performance
    print("📊 Evaluating model performance...")

    results_file = os.path.join(DATA_DIR, "classification_results.txt")

    with open(results_file, "w") as f:
        f.write("📊 Machine Learning Classification Results\n")
        f.write("=" * 50 + "\n")

        for model_name, model in [("Random Forest", rf_model), ("XGBoost", xgb_model)]:
            print(f"🔹 Evaluating {model_name}...")
            y_pred = model.predict(X_test)
            accuracy = accuracy_score(y_test, y_pred)
            report = classification_report(y_test, y_pred)

            print(f"✅ {model_name} Test Set Accuracy: {accuracy:.4f}\n")
            print("\n🔹 Classification Report:\n", report)

            # Save results to file
            f.write(f"\n🔹 {model_name} Results:\n")
            f.write(f"✅ Test set accuracy: {accuracy:.4f}\n\n")
            f.write("🔹 Classification Report:\n")
            f.write(report + "\n")

            # Step 7: Compute Confusion Matrix
            print(f"📌 Computing confusion matrix for {model_name}...")
            cm = confusion_matrix(y_test, y_pred)

            f.write("🔹 Confusion Matrix:\n")
            f.write(np.array2string(cm, separator=', ') + "\n")

            # Save confusion matrix plot
            plt.figure(figsize=(12, 10))
            sns.heatmap(cm, cmap="viridis", cbar=True, 
                        xticklabels=np.unique(y_test), 
                        yticklabels=np.unique(y_test), 
                        annot=True, fmt="d")
            plt.xlabel("Predicted Label")
            plt.ylabel("True Label")
            plt.title(f"{model_name} Confusion Matrix")
            plt.xticks(rotation=90)
            plt.yticks(rotation=0)
            plt.tight_layout()

            # Save the confusion matrix plot
            conf_matrix_path = os.path.join(DATA_DIR, f"confusion_matrix_{model_name.replace(' ', '_')}.png")
            plt.savefig(conf_matrix_path)
            print(f"📁 Confusion matrix saved to {conf_matrix_path}")

            plt.close()

    print(f"📁 Classification results saved to {results_file}")
    print("🎯 Machine learning classification completed!")

if __name__ == "__main__":
    supervised_learning()

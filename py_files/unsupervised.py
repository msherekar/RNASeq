import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import os
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.mixture import GaussianMixture
import hdbscan
from sklearn.preprocessing import StandardScaler

# Get the absolute path of the RNASEQ main directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE_DIR, "outputs")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

def unsupervised_learning():
    print("🔹 Performing unsupervised learning for novel subpopulation discovery...")

    # Ensure the necessary input file exists
    clustered_data_path = os.path.join(DATA_DIR, "clustered_data.h5ad")
    if not os.path.exists(clustered_data_path):
        raise FileNotFoundError(f"❌ Clustered data file not found: {clustered_data_path}")
    
    # Step 1: Load clustered RNA-seq data
    print("📥 Loading clustered RNA-seq data...")
    adata = sc.read_h5ad(clustered_data_path)

    # Reduce number of genes (select only highly variable genes)
    print("🔹 Selecting highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable]

    # Convert sparse matrix safely
    if isinstance(adata.X, np.ndarray):
        X = adata.X
    else:
        X = adata.X.toarray()

    # Standardize the data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    device = torch.device("mps") if torch.backends.mps.is_available() else torch.device("cpu")
    # Define batch generator to prevent memory issues
    def batch_loader(X, batch_size=512):
        for i in range(0, X.shape[0], batch_size):
            yield torch.tensor(X[i:i+batch_size], dtype=torch.float32).to(device)

    # Define VAE model
    class VAE(nn.Module):
        def __init__(self, input_dim, latent_dim=5):  # Reduce latent space to 5 dimensions
            super(VAE, self).__init__()
            self.encoder = nn.Sequential(
                nn.Linear(input_dim, 256),
                nn.ReLU(),
                nn.Linear(256, 128),
                nn.ReLU(),
                nn.Linear(128, latent_dim * 2)  # Mean and log variance
            )
            self.decoder = nn.Sequential(
                nn.Linear(latent_dim, 128),
                nn.ReLU(),
                nn.Linear(128, 256),
                nn.ReLU(),
                nn.Linear(256, input_dim),
                nn.Sigmoid()
            )

        def reparameterize(self, mu, log_var):
            std = torch.exp(0.5 * log_var)
            eps = torch.randn_like(std)
            return mu + eps * std

        def forward(self, x):
            x = self.encoder(x)
            mu, log_var = x.chunk(2, dim=-1)  # Split mean and log variance
            z = self.reparameterize(mu, log_var)
            x_recon = self.decoder(z)
            return x_recon, mu, log_var

    # Train VAE
    def train_vae(X, latent_dim=5, epochs=50, batch_size=512, learning_rate=1e-3):
        #device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        device = torch.device("mps") if torch.backends.mps.is_available() else torch.device("cpu")
        input_dim = X.shape[1]
        model = VAE(input_dim, latent_dim).to(device)
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)

        loss_fn = nn.MSELoss()

        for epoch in range(epochs):
            model.train()
            epoch_loss = 0
            for batch in batch_loader(X, batch_size):
                optimizer.zero_grad()
                x_recon, mu, log_var = model(batch)
                recon_loss = loss_fn(x_recon, batch)

                # KL Divergence loss
                kl_loss = -0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())
                loss = recon_loss + kl_loss

                loss.backward()
                optimizer.step()
                epoch_loss += loss.item()

            print(f"Epoch {epoch+1}/{epochs}, Loss: {epoch_loss:.4f}")

        model.eval()
        latent_repr = []
        with torch.no_grad():
            for batch in batch_loader(X, batch_size):
                _, mu, _ = model(batch)
                latent_repr.append(mu.cpu().numpy())
        
        return np.vstack(latent_repr)

    # Generate latent space representation
    print("🔍 Training VAE for latent representation...")
    latent_space = train_vae(X_scaled, latent_dim=5)

    # Apply HDBSCAN clustering
    print("📌 Applying HDBSCAN on latent space...")
    hdbscan_clusterer = hdbscan.HDBSCAN(min_cluster_size=10, min_samples=5)
    hdbscan_labels = hdbscan_clusterer.fit_predict(latent_space)

    # Apply Gaussian Mixture Model clustering
    print("📌 Applying Gaussian Mixture Model on latent space...")
    gmm = GaussianMixture(n_components=10, covariance_type='full', random_state=42)
    gmm_labels = gmm.fit_predict(latent_space)

    # Save new clusters to AnnData object
    adata.obs["HDBSCAN_Cluster"] = hdbscan_labels.astype(str)
    adata.obs["GMM_Cluster"] = gmm_labels.astype(str)
    updated_data_path = os.path.join(DATA_DIR, "clustered_data_with_new_clusters.h5ad")
    adata.write(updated_data_path)
    print(f"📁 Updated clustered data saved to {updated_data_path}")

    # Compute UMAP if not already computed
    if "X_umap" not in adata.obsm:
        print("📌 Computing UMAP...")
        sc.pp.neighbors(adata, use_rep="X_pca")  # Build nearest-neighbor graph
        sc.tl.umap(adata)  # Compute UMAP

    # Plot UMAP colored by HDBSCAN Clusters
    sc.pl.umap(adata, color="HDBSCAN_Cluster", title="UMAP of HDBSCAN Clusters", show=False)
    hbd_path = os.path.join(RESULTS_DIR, "hbd_clusters.png")
    plt.savefig(hbd_path)
    plt.close()
    
    # Plot UMAP colored by GMM Clusters
    sc.pl.umap(adata, color="GMM_Cluster", title="UMAP of GMM Clusters", show=False)
    gmm_path = os.path.join(RESULTS_DIR, "gmm_clusters.png")
    plt.savefig(gmm_path)
    plt.close()

    print("🎯 Unsupervised learning and clustering completed!")

if __name__ == "__main__":
    unsupervised_learning()

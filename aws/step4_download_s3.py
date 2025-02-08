# step4_download_s3.py
import boto3
import config

def download_results():
    s3 = boto3.client("s3")
    print("📥 Downloading results from S3...")
    s3.download_file(config.S3_BUCKET, "results/gene_counts_matrix.csv", "gene_counts_matrix.csv")
    print("✅ Download complete.")

if __name__ == "__main__":
    download_results()

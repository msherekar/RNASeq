# step4_download_s3.py
import boto3
import config
import logging

def download_results():
    s3 = boto3.client("s3")
    logging.info("📥 Downloading results from S3...")
    try:
        s3.download_file(config.S3_BUCKET, "results/gene_counts_matrix.csv", "gene_counts_matrix.csv")
        logging.info("✅ Download complete.")
    except Exception as e:
        logging.error("❌ Error downloading results from S3: %s", e)

if __name__ == "__main__":
    download_results()

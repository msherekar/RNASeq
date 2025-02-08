# step1_upload_s3.py
import boto3
import config

def upload_fastq_to_s3():
    s3 = boto3.client("s3")
    print("📤 Uploading FASTQ files to S3...")
    for file in config.FASTQ_FILES:
        s3.upload_file(file, config.S3_BUCKET, f"input/{file}")
    print("✅ FASTQ files uploaded.")

if __name__ == "__main__":
    upload_fastq_to_s3()
    print("🚀 Uploading FASTQ files to S3...")
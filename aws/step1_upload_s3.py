# step1_upload_gz_to_s3.py
import os
import boto3
import config

def upload_to_s3():
    """
    Uploads all .gz FASTQ files specified in config.FASTQ_FILES to S3.
    
    Files are uploaded as-is to S3 to save on bandwidth,
    with the intention of extracting them later in the cloud.
    """
    s3 = boto3.client("s3")
    
    # Iterate over each file defined in the config FASTQ_FILES list
    for fastq_file in config.FASTQ_FILES:
        # Construct the local file path. Assuming this script is run from RNASEQ/aws,
        # and the data folder is in RNASEQ/data.
        file_path = os.path.join("..", "data", fastq_file)
        
        # Define the S3 key. Here we place the file in a "data" folder in the S3 bucket.
        s3_key = os.path.join("data", fastq_file)
        
        print(f"📤 Uploading {file_path} as {s3_key} to S3 bucket {config.S3_BUCKET}...")
        try:
            s3.upload_file(file_path, config.S3_BUCKET, s3_key)
            print(f"✅ {fastq_file} uploaded successfully.")
        except Exception as e:
            print(f"❌ Failed to upload {fastq_file}: {e}")

if __name__ == "__main__":
    upload_to_s3()

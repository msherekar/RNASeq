# main.py
import logging
import step1_upload_s3
import step2_launch_ec2
import step3_run_pipeline
import step4_download_s3
import step5_terminate_ec2

def setup_logging():
    # Set up the root logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # File handler logs to pipeline.log
    fh = logging.FileHandler("pipeline.log")
    fh.setLevel(logging.INFO)
    
    # Console handler logs to stdout
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    # Formatter for both handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    # Add handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

if __name__ == "__main__":
    setup_logging()

    logging.info("🚀 Starting AWS RNA-Seq Pipeline")
    print("🚀 Starting AWS RNA-Seq Pipeline")

    # Step 1: Upload FASTQ files to S3
    logging.info("🔼 Step 1: Uploading FASTQ files to S3...")
    print("🔼 Step 1: Uploading FASTQ files to S3...")
    step1_upload_s3.upload_to_s3()
    logging.info("✅ Step 1 complete")
    print("✅ Step 1 complete")

    # Step 2: Launch EC2 instance
    logging.info("🌐 Step 2: Launching EC2 instance...")
    print("🌐 Step 2: Launching EC2 instance...")
    instance_id, public_ip = step2_launch_ec2.launch_ec2("First_Pipeline")
    logging.info("✅ Step 2 complete - Instance ID: %s, Public IP: %s", instance_id, public_ip)
    print("✅ Step 2 complete")

    # Step 3: Run RNA-Seq pipeline on EC2
    logging.info("⚙️ Step 3: Running RNA-Seq pipeline on EC2...")
    print("⚙️ Step 3: Running RNA-Seq pipeline on EC2...")
    step3_run_pipeline.run_rna_seq_on_ec2(public_ip)
    logging.info("✅ Step 3 complete")
    print("✅ Step 3 complete")

    # Step 4: Download results from S3
    logging.info("📥 Step 4: Downloading results from S3...")
    print("📥 Step 4: Downloading results from S3...")
    step4_download_s3.download_results()
    logging.info("✅ Step 4 complete")
    print("✅ Step 4 complete")

    # Step 5: Terminate EC2 instance
    logging.info("💀 Step 5: Terminating EC2 instance...")
    print("💀 Step 5: Terminating EC2 instance...")
    step5_terminate_ec2.terminate_ec2()
    logging.info("✅ Step 5 complete")
    print("✅ Step 5 complete")

    logging.info("🎉 RNA-Seq pipeline completed successfully!")
    print("🎉 RNA-Seq pipeline completed successfully!")

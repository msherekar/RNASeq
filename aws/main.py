# main.py
import step1_upload_s3
import step2_launch_ec2
import step3_run_pipeline
import step4_download_s3
import step5_terminate_ec2

if __name__ == "__main__":
    print("🚀 Starting AWS RNA-Seq Pipeline")

    # Step 1: Upload FASTQ files to S3
    step1_upload_s3.upload_fastq_to_s3()

    # Step 2: Launch EC2 instance
    instance_id, public_ip = step2_launch_ec2.launch_ec2()

    # Step 3: Run RNA-Seq pipeline on EC2
    step3_run_pipeline.run_rna_seq_on_ec2(public_ip)

    # Step 4: Download results from S3
    step4_download_s3.download_results()

    # Step 5: Terminate EC2 instance
    step5_terminate_ec2.terminate_ec2()

    print("🎉 RNA-Seq pipeline completed successfully!")

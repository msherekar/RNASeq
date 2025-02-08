# step3_run_pipeline.py
import paramiko
import config

def run_rna_seq_on_ec2(public_ip):
    print("🔑 Connecting to EC2 instance...")
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(public_ip, username="ubuntu", key_filename=config.SSH_KEY_PATH)

    print("⚙️ Running RNA-Seq pipeline...")
    commands = f"""
    sudo apt update && sudo apt install -y fastqc fastp star subread awscli;
    aws s3 cp s3://{config.S3_BUCKET}/input/sample_R1.fastq.gz .;
    aws s3 cp s3://{config.S3_BUCKET}/input/sample_R2.fastq.gz .;
    fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o qc_reports/;
    fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz;
    STAR --genomeDir {config.GENOME_INDEX} --readFilesIn trimmed_R1.fastq.gz trimmed_R2.fastq.gz --runThreadN 8 --outSAMtype BAM SortedByCoordinate;
    featureCounts -T 8 -a {config.ANNOTATION_FILE} -o gene_counts.txt sample_Aligned.sortedByCoord.out.bam;
    aws s3 cp gene_counts.txt s3://{config.S3_BUCKET}/results/gene_counts_matrix.csv;
    """
    stdin, stdout, stderr = ssh.exec_command(commands)
    print(stdout.read().decode(), stderr.read().decode())
    print("✅ RNA-Seq pipeline completed!")

if __name__ == "__main__":
    with open("instance_info.txt") as f:
        _, public_ip = f.read().splitlines()
    run_rna_seq_on_ec2(public_ip)

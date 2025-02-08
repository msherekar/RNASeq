# config.py
AWS_REGION = "us-east-1"
S3_BUCKET = "my-rna-seq-bucket"
INSTANCE_TYPE = "r5.4xlarge"
KEY_NAME = "my-key"
SECURITY_GROUP = "my-security-group"
IAM_ROLE = "MyEC2RoleForRNASeq"
AMI_ID = "ami-0abcdef1234567890"
SSH_KEY_PATH = "my-key.pem"

# RNA-Seq pipeline files
FASTQ_FILES = ["sample_R1.fastq.gz", "sample_R2.fastq.gz"]
GENOME_INDEX = "genome_index/"
ANNOTATION_FILE = "annotation.gtf"

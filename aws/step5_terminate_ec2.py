# step5_terminate_ec2.py
import boto3

def terminate_ec2():
    with open("instance_info.txt") as f:
        instance_id, _ = f.read().splitlines()
    
    print("🛑 Terminating EC2 instance...")
    ec2 = boto3.client("ec2")
    ec2.terminate_instances(InstanceIds=[instance_id])
    print("✅ EC2 instance terminated.")

if __name__ == "__main__":
    terminate_ec2()

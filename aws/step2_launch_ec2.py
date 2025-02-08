# step2_launch_ec2.py
import boto3
import config

def launch_ec2():
    ec2 = boto3.resource("ec2", region_name=config.AWS_REGION)
    print("🚀 Launching EC2 instance...")
    instance = ec2.create_instances(
        ImageId=config.AMI_ID,
        InstanceType=config.INSTANCE_TYPE,
        KeyName=config.KEY_NAME,
        SecurityGroups=[config.SECURITY_GROUP],
        IamInstanceProfile={'Name': config.IAM_ROLE},
        MinCount=1,
        MaxCount=1
    )[0]

    instance.wait_until_running()
    instance.load()
    print(f"✅ EC2 instance launched: {instance.public_ip_address}")
    return instance.instance_id, instance.public_ip_address

if __name__ == "__main__":
    instance_id, public_ip = launch_ec2()
    with open("instance_info.txt", "w") as f:
        f.write(f"{instance_id}\n{public_ip}")

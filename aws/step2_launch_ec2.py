import boto3
import config

def launch_ec2(instance_name):
    """
    Launches an EC2 instance with a specified name.
    
    Args:
        instance_name (str): The name to assign to the EC2 instance.
    
    Returns:
        tuple: (instance_id, public_ip)
    """
    ec2 = boto3.resource("ec2", region_name=config.AWS_REGION)
    print(f"🚀 Launching EC2 instance with name: {instance_name}...")

    instance = ec2.create_instances(
        ImageId=config.AMI_ID,
        InstanceType=config.INSTANCE_TYPE,
        KeyName=config.KEY_NAME,
        SecurityGroups=[config.SECURITY_GROUP],
        IamInstanceProfile={'Name': config.IAM_ROLE},
        MinCount=1,
        MaxCount=1,
        TagSpecifications=[
            {
                "ResourceType": "instance",
                "Tags": [
                    {"Key": "Name", "Value": instance_name}
                ]
            }
        ]
    )[0]

    instance.wait_until_running()
    instance.load()

    print(f"✅ EC2 instance launched: {instance_name}")
    print(f"📌 Instance ID: {instance.instance_id}")
    print(f"🌍 Public IP: {instance.public_ip_address}")

    return instance.instance_id, instance.public_ip_address, instance_name

if __name__ == "__main__":
    instance_name = "My-RNASeq-Instance"  # Change or pass dynamically
    instance_id, public_ip, instance_name = launch_ec2(instance_name)
    
    with open("instance_info.txt", "w") as f:
        f.write(f"Instance Name: {instance_name}\nInstance ID: {instance_id}\nPublic IP: {public_ip}")

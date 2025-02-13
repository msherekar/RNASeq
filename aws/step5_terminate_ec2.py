import boto3

def get_instance_id_by_name(instance_name):
    """
    Retrieves the instance ID of an EC2 instance by its name tag.

    Args:
        instance_name (str): The name of the EC2 instance to find.

    Returns:
        str: The instance ID if found, otherwise None.
    """
    ec2 = boto3.client("ec2")

    response = ec2.describe_instances(
        Filters=[
            {"Name": "tag:Name", "Values": [instance_name]},
            {"Name": "instance-state-name", "Values": ["running", "pending", "stopping", "stopped"]}
        ]
    )

    instances = [
        res["Instances"][0]["InstanceId"]
        for res in response["Reservations"] if res["Instances"]
    ]

    if not instances:
        print(f"❌ No running/stopped instances found with name: {instance_name}")
        return None

    return instances[0]  # Return the first matching instance ID

def terminate_ec2():
    """
    Reads instance details from `instance_info.txt` and terminates the instance.
    """
    try:
        with open("instance_info.txt", "r") as f:
            lines = f.read().splitlines()
            instance_name = lines[0].split(": ")[1]  # Extract name from "Instance Name: <name>"
            instance_id = lines[1].split(": ")[1]  # Extract ID from "Instance ID: <id>"

        print(f"🛑 Terminating EC2 instance: {instance_name} (ID: {instance_id})...")
        ec2 = boto3.client("ec2")
        ec2.terminate_instances(InstanceIds=[instance_id])
        print(f"✅ EC2 instance {instance_name} terminated.")

    except FileNotFoundError:
        print("❌ instance_info.txt not found. Cannot terminate instance.")
    except IndexError:
        print("❌ instance_info.txt format is incorrect. Ensure it contains instance details.")

if __name__ == "__main__":
    terminate_ec2()

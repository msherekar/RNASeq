# step5_terminate_ec2.py
import boto3
import logging

def get_instance_id_by_name(instance_name):
    """
    Retrieves the instance ID of an EC2 instance by its name tag.

    Args:
        instance_name (str): The name of the EC2 instance to find.

    Returns:
        str: The instance ID if found, otherwise None.
    """
    ec2 = boto3.client("ec2")
    try:
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
            logging.error("❌ No running/stopped instances found with name: %s", instance_name)
            return None

        return instances[0]  # Return the first matching instance ID

    except Exception as e:
        logging.error("❌ An error occurred while retrieving instance ID for %s: %s", instance_name, e)
        return None

def terminate_ec2():
    """
    Reads instance details from `instance_info.txt` and terminates the instance.
    Assumes instance_info.txt contains two lines:
      Line 1: instance_id
      Line 2: public_ip
    """
    try:
        with open("instance_info.txt", "r") as f:
            lines = f.read().splitlines()
            if len(lines) < 1:
                logging.error("instance_info.txt does not contain the required instance details.")
                return
            instance_id = lines[0].strip()

        logging.info("🛑 Terminating EC2 instance with ID: %s...", instance_id)
        ec2 = boto3.client("ec2")
        response = ec2.terminate_instances(InstanceIds=[instance_id])
        logging.info("✅ Termination initiated for instance ID: %s. Response: %s", instance_id, response)

    except FileNotFoundError:
        logging.error("❌ instance_info.txt not found. Cannot terminate instance.")
    except IndexError:
        logging.error("❌ instance_info.txt format is incorrect. Ensure it contains instance details.")
    except Exception as e:
        logging.error("❌ An error occurred while terminating the EC2 instance: %s", e)

if __name__ == "__main__":
    # Configure logging if not already configured by the main pipeline.
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    terminate_ec2()

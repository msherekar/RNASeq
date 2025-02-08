# utils.py
import requests
import os
import config

def download_file(url, output_path):
    """Download a file from a given URL in chunks."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(config.DOWNLOAD_CHUNK_SIZE):
                file.write(chunk)
    
    print(f"✅ Download completed: {output_path}")

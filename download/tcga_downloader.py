# tcga_downloader.py
import requests
import os
import config
import utils

def fetch_tcga_files(case_id, data_type="RNA-Seq"):
    """Fetch available RNA-Seq files for a given TCGA case ID."""
    url = f"{config.TCGA_BASE_URL}/files"
    params = {
        "filters": '{"op":"and","content":[{"op":"in","content":{"field":"cases.case_id","value":["' + case_id + '"]}},{"op":"in","content":{"field":"data_type","value":["' + data_type + '"]}}]}',
        "format": "JSON",
        "size": "100"
    }
    response = requests.get(url, params=params)
    response.raise_for_status()
    files = response.json().get("data", {}).get("hits", [])
    return files

def download_tcga_file(file_id):
    """Download a TCGA file using its file ID."""
    url = f"{config.TCGA_BASE_URL}/data/{file_id}"
    output_path = os.path.join(config.TCGA_DOWNLOAD_DIR, f"{file_id}.gz")
    
    utils.download_file(url, output_path)
    print(f"✅ TCGA file downloaded: {output_path}")

if __name__ == "__main__":
    case_id = "TCGA-02-0047"  # Example case ID
    files = fetch_tcga_files(case_id)
    
    if files:
        for file in files[:2]:  # Download first two files as an example
            download_tcga_file(file["file_id"])
    else:
        print("No TCGA RNA-Seq files found for this case.")

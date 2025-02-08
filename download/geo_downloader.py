# geo_downloader.py
import requests
import os
import config
import utils

def fetch_geo_series(accession_id):
    """Fetch metadata for a GEO series using its accession ID."""
    url = f"{config.GEO_BASE_URL}?acc={accession_id}&format=json"
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

def download_geo_data(accession_id):
    """Download RNA-Seq data from GEO using an accession ID."""
    output_path = os.path.join(config.GEO_DOWNLOAD_DIR, f"{accession_id}.tar")
    
    geo_metadata = fetch_geo_series(accession_id)
    if "supplementary_files" in geo_metadata:
        file_url = geo_metadata["supplementary_files"][0]["url"]
        utils.download_file(file_url, output_path)
        print(f"✅ GEO dataset downloaded: {output_path}")
    else:
        print("No supplementary files available for this GEO dataset.")

if __name__ == "__main__":
    accession_id = "GSE134131"  # Example GEO dataset
    download_geo_data(accession_id)

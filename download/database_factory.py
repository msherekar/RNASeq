# database_factory.py
from tcga_downloader import download_tcga_file
from geo_downloader import download_geo_data

def download_rna_seq_data(database, identifier):
    """Dispatch function to download RNA-Seq data from the correct database."""
    if database.lower() == "tcga":
        download_tcga_file(identifier)
    elif database.lower() == "geo":
        download_geo_data(identifier)
    else:
        print(f"❌ Unsupported database: {database}")


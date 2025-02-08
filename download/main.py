# main.py
import database_factory

if __name__ == "__main__":
    print("🚀 RNA-Seq Data Downloader")
    
    # Example: Fetch data from TCGA
    database = "TCGA"  # Change to "GEO" for GEO datasets
    identifier = "TCGA-02-0047"  # Change to accession ID for GEO
    
    database_factory.download_rna_seq_data(database, identifier)

    print("🎉 RNA-Seq data download completed!")

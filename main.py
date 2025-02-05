# RNA-Seq Pipeline: Main Script (main.py)

import os
import logging
import time
import traceback

# Import modules from py_files directory
from py_files import (
    load_data,
    quality_control,
    normalization,
    dimensionality_reduction,
    clustering,
    differential_expression,
    pathway_analysis,
    visualizations,
    supervised, 
    unsupervised,  
)

# Configure logging
logging.basicConfig(
    filename="pipeline.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

def run_pipeline():
    start_time = time.time()
    logging.info("Starting RNA-Seq analysis pipeline...")
    print("\n🚀 Starting RNA-Seq Analysis Pipeline...\n")

    try:
        # Step 1: Load Data
        print("🔹 Step 1: Loading Data...")
        step_start = time.time()
        load_data.load_expression_matrix()
        logging.info(f"Step 1 completed in {time.time() - step_start:.2f} seconds")

        # Step 2: Quality Control
        print("🔹 Step 2: Performing Quality Control...")
        step_start = time.time()
        quality_control.quality_control()
        logging.info(f"Step 2 completed in {time.time() - step_start:.2f} seconds")

        # Step 3: Normalization
        print("🔹 Step 3: Normalizing Data...")
        step_start = time.time()
        normalization.normalize_data()
        logging.info(f"Step 3 completed in {time.time() - step_start:.2f} seconds")

        # Step 4: Dimensionality Reduction
        print("🔹 Step 4: Performing Dimensionality Reduction...")
        step_start = time.time()
        dimensionality_reduction.reduce_dimensions()
        logging.info(f"Step 4 completed in {time.time() - step_start:.2f} seconds")

        # Step 5: Clustering
        print("🔹 Step 5: Clustering Cells...")
        step_start = time.time()
        clustering.cluster_cells()
        logging.info(f"Step 5 completed in {time.time() - step_start:.2f} seconds")

        # Step 6: Differential Expression Analysis
        print("🔹 Step 6: Performing Differential Expression Analysis...")
        step_start = time.time()
        differential_expression.differential_expression()
        logging.info(f"Step 6 completed in {time.time() - step_start:.2f} seconds")

        # Step 7: Pathway Analysis
        print("🔹 Step 7: Running Pathway Analysis...")
        step_start = time.time()
        pathway_analysis.pathway_analysis()
        logging.info(f"Step 7 completed in {time.time() - step_start:.2f} seconds")

        # Step 8: Generate Visualizations
        print("🔹 Step 8: Generating Visualizations...")
        step_start = time.time()
        visualizations.generate_plots()
        logging.info(f"Step 8 completed in {time.time() - step_start:.2f} seconds")

        # Step 9: Supervised Learning
        print("🔹 Step 9: Running Supervised Learning...")
        step_start = time.time()
        supervised.supervised_learning()  
        logging.info(f"Step 9 completed in {time.time() - step_start:.2f} seconds")

        # Step: 10: Unsupervised Learning
        print("🔹 Step 10: Running Unsupervised Learning...")
        step_start = time.time()
        unsupervised.unsupervised_learning()
        logging.info(f"Step 10 completed in {time.time() - step_start:.2f} seconds")

        # Pipeline completed
        total_time = time.time() - start_time
        logging.info(f"RNA-Seq pipeline completed successfully in {total_time:.2f} seconds!")
        print(f"\n✅ Pipeline completed successfully in {total_time:.2f} seconds!\n")

    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        print("\n❌ An error occurred during pipeline execution.")
        print(traceback.format_exc())  # Print full traceback for debugging

if __name__ == "__main__":
    run_pipeline()

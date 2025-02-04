# RNA-Seq Pipeline: Main Script (main.py)

import os
import logging
import time
import load_data
import quality_control
import normalization
import dimensionality_reduction
import clustering
import differential_expression
import pathway_analysis
import visualizations

# Configure logging
logging.basicConfig(
    filename="outputs/pipeline.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

def run_pipeline():
    start_time = time.time()
    logging.info("Starting RNA-Seq analysis pipeline...")
    try:
        # Step 1: Load Data
        step_start = time.time()
        logging.info("Loading data...")
        load_data.load_expression_matrix()
        logging.info(f"Step 1 completed in {time.time() - step_start:.2f} seconds")
        
        # Step 2: Quality Control
        step_start = time.time()
        logging.info("Running quality control...")
        quality_control.quality_control()
        logging.info(f"Step 2 completed in {time.time() - step_start:.2f} seconds")
        
        # Step 3: Normalization
        step_start = time.time()
        logging.info("Normalizing data...")
        normalization.normalize_data()
        logging.info(f"Step 3 completed in {time.time() - step_start:.2f} seconds")
        
        # Step 4: Dimensionality Reduction
        step_start = time.time()
        logging.info("Performing dimensionality reduction...")
        dimensionality_reduction.reduce_dimensions()
        logging.info(f"Step 4 completed in {time.time() - step_start:.2f} seconds")
        
        # Step 5: Clustering
        step_start = time.time()
        logging.info("Clustering cells...")
        clustering.cluster_cells()
        logging.info(f"Step 5 completed in {time.time() - step_start:.2f} seconds")
        
        # Step 6: Differential Expression Analysis
        step_start = time.time()
        logging.info("Performing differential expression analysis...")
        differential_expression.differential_expression()
        logging.info(f"Step 6 completed in {time.time() - step_start:.2f} seconds")
        
        # Step 7: Pathway Analysis
        step_start = time.time()
        logging.info("Running pathway analysis...")
        pathway_analysis.pathway_analysis()
        logging.info(f"Step 7 completed in {time.time() - step_start:.2f} seconds")
        
        # Step 8: Generate Visualizations
        step_start = time.time()
        logging.info("Generating visualizations...")
        visualizations.generate_plots()
        logging.info(f"Step 8 completed in {time.time() - step_start:.2f} seconds")
        
        total_time = time.time() - start_time
        logging.info(f"RNA-Seq pipeline completed successfully in {total_time:.2f} seconds!")
        print(f"Pipeline completed in {total_time:.2f} seconds")
    
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    run_pipeline()

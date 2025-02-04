=================================
Human Cortex
=================================

This data set includes single-nucleus transcriptomes from 76,533 total nuclei derived from 2 post-mortem human brain specimens, to survey 
cell type diversity in the primary motor cortex (M1C or M1). The data were generated as part of a BICCN collaboration to characterize cell type 
diversity in M1 across species and data modalities. In total, 127 transcriptomic cell types were identified. 

Samples were processed using the 10x Chromium Single Cell 3â€™ Reagent Kit v3. 10x chip loading and sample processing was done according to the 
manufacturerâ€™s protocol. Gene expression was quantified using the default 10x Cell Ranger v3 pipeline except substituting the curated genome 
annotation used for SMART-seq v4 quantification. Introns were annotated as â€œmRNA,â€ and intronic reads were included in expression quantification.


Gene expression data matrices (matrix.csv)
    This file csv contains one row for every cell in the dataset and one column for every gene sequenced. The values of the matrix represent total reads (introns + exons) 
	for that gene (column) for that cell (row).

		
Trimmed means (trimmed-means.csv) 
	A table of trimmed means for each gene (rows) in each cluster (columns).  
	Trimmed means are calculated by first normalizing gene expression as follows: norm_data = log2(CPM(exons+introns)), 
	and then calculating the average expression of the middle 50% of the data 
	(e.g., after excluding the 25% highest and 25% lowest expression values) independently for each gene and each cluster.
	The first row lists the cluster name (cluster_label), which matches the cell type alias shown in the Transcriptomic Explorer.
	The first column lists the unique gene identifier (gene), which in most cases is the gene symbol.


Cell metadata (metadata.csv)
* Each item of this table (except "sample_name") has three columns:
	[item]_label
		Name of the item (e.g., "V1C" would be an example under "brain_region_label")
	[item]_order
		Order that the item will be displayed on the Transcriptomics Explorer 
	[item]_color
		Color that the item will be displayed on the Transcriptomics Explorer 

* Items in the sample information table:
	sample_name
		Unique sample identifier
	cluster
		Cell type cluster name
	cell_type_accession
		Cell type accession ID (see https://portal.brain-map.org/explore/classes/nomenclature for details)
	cell_type_alias
		Cell type alias (see https://portal.brain-map.org/explore/classes/nomenclature for details).  This is the same as "cluster".
	cell_type_alt_alias
		Cell type alternative alias, if any (see https://portal.brain-map.org/explore/classes/nomenclature for details)
	cell_type_designation
		Cell type label (see https://portal.brain-map.org/explore/classes/nomenclature for details)
	class
		Broad cell class (for example, "GABAergic", "Non-neuronal", and "Glutamatergic")
	subclass
		Cell type subclass (for example, "SST", "L6 CT", and "Astrocyte")
	external_donor_name
		Unique (de-identified) identifier for each human donor
	donor_sex
		Biological sex of the donor
	cortical_layer
		Cortical layer targeted for sampling
	region
		Brain region targeted for sampling
	full_genotype
		(These slots are left blank for human)
	
	
TSNE coordinates (2d-coordinates.zip)
t-Distributed Stochastic Neighbor Embedding (t-SNE) coordinates for each sample shown on the Transcriptomics Explorer.  t-SNE is a method for dimensionality reduction of gene expression that is  well suited for data visualization (as of 1 October 2019, a comprehensive t-SNE resource is available here: https://lvdmaaten.github.io/tsne/)
	sample_name
		Unique sample identifier
	tsne_1
		First t-SNE coordinate
	tsne_2
		Second t-SNE coordinate

		
Taxonomy of clusters (dend.json)
	Serialized cluster hierarchy with all node information embedded in json format.
	The dendrogram shown at the top of the Transcriptomics Explorer, including the underlying cluster order, is derived from this file.
	

Taxonomy information (taxonomy.txt)
	Tracking taxonomy meta-data is critical for reproducibility.  This file is a draft of taxonomy meta-data to be stored.  See the "Tracking taxonomies" section at https://portal.brain-map.org/explore/classes/nomenclature for details of each descriptor.
	

Gene information (**STORED ELSEWHERE**)
* To access this file, please use the following link: http://celltypes.brain-map.org/api/v2/well_known_file_download/694416044
* Within that zip file, the gene information is located in "human_MTG_2018-06-14_genes-rows.csv".  All other files can be ignored.
	gene
		Gene symbol
	chromosome
		Chromosome location of gene
	entrez_id
		NCBI Entrez ID
	gene_name
		Gene name
	mouse_homologenes
		Mouse ortholog

		
Gene ".gtf" file (**STORED ELSEWHERE**)
* To access this file, please use the following link: http://celltypes.brain-map.org/api/v2/well_known_file_download/502175284
.gtf is a standard format for localizing various aspects of transcripts within a specific genome and information about this format is plentiful.
As of 1 October 2019, one active link describing this format is here: https://www.gencodegenes.org/pages/data_format.html
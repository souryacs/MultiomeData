RNA velocity analysis of multi-omic data (multivelo)
==================================================

multivelo (https://github.com/welch-lab/MultiVelo) is a package to perform RNA velocity analysis from multi-omic single cell (scRNA-seq + scATC-seq) data.

Assumes that the multi-omic data is processed by Signac.

	- *1_Extract_Seurat_WNN_Out.R*: R Script to extract the clustering information from the input Signac object. 

		- *Note*: User needs to check the code and edit the clustering resolution (and corresponding fields). The current code assumes the clustering resolution = 0.3. For any other resolution, user needs to edit the values and respective fields (check lines 32-42 of the script).

		- User needs to provide the input Signac object file and the output directory name to store the output files

	- *Run_velocyto_script.sh*: Script to run Velocyto to generate the .loom files from the Cellranger-arc output .bam file. This .loom file is an input to multielo.

		- User needs to check individual parameters of this script and edit accordingly.

	- *2_MultiVelo_Script.py*: Script to run multivelo

		- User needs to check the parameter section, and edit. 
		- They need to provide the input files according to the outputs of the previous two scripts.
		- Same as scvelo, user needs to provide one target gene list file (first column containing this list) whose trajectory need to be examined. 

		- Output plots will be stored under the specied output directory (check the parameter section)

		- Check the multivelo manual for the details of the output files.

Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org


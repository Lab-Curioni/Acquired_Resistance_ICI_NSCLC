# manuscript_resistance_to_immunotherapy
Code related to the manuscript resistance to immunotherapy NSCLC by hiltbrunner et al.



## DNA analysis
Snakemake was used to run the pipeline described in the following paragraph.
The preprocessing of the dna analysis makes use of the public pipeline NGS-pipe. For this project the git commit c738bfc has been used. Please refer to the public wiki here https://github.com/cbg-ethz/NGS-pipe for more information on this part.

The inhouse pipeline pipeline_dna (commit 00cb440e ) has been built on top of this to do the postprocessing. 

Both commits have been copied into this git into the correct folder structure [pipeline_dna_commit_00cb440e and subfolder NGS-pipe] that the pipeline needs in order to run. The versions can be found in the versions.md file in this folder.

### WES analysis

Project specific files used: 
- pipeline_dna.wes.SNAKE
- pipeline_dna.wes_snake.py
- pipeline_dna.config_wes.json

### WGS analysis

Project specific files used:
- pipeline_dna.wgs.SNAKE
- pipeline_dna.wgs_snake.py
- pipeline_dna.config_wgs.json



## RNA analysis
Snakemake was used to run the pipeline described in the following paragraph.
For the rna analysis an inhouse pipeline (commit 4c35df23) has been used to perform the analysis. 
The commit has been copied into this git (pipeline_rna_commit_4c35df23).
The software versions used can be found in the versions.md file in this folder. 

### Gene set variation
R script used to determine gene set enrichment scores and generate the heatmaps can be found here: gene_set_variation_analysis .

### log2FC
R script used to determine log2 fold changes in gene expression between tumor pairs can be found here: RNAseq_log2FC .



## Circos plots
The scripts that have been used to create the plots can be found in the folder circos_plots.

The circos plots were generated using circos:

'''
Circos: An information aesthetic for comparative genomics
Krzywinski et al., Genome Res. Published in Advance June 18, 2009, doi:10.1101/gr.092759.109
''' 

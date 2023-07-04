/path/to/utilities/sharedPrograms/snakemake/snakemake_4.8.0/bin/snakemake --notemp -s /path/to//nexusGit/pipeline_rna/config_snake_files/tcga_rna.snake  --configfile /path/to//nexusGit/pipeline_rna/config_snake_files/tcga_config.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 120 -p -k -n
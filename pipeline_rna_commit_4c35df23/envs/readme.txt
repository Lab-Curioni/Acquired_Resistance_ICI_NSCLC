###   storing the environments as textfile   ###
conda env export | grep -v prefix > icgc_pipeline_conda_env.yml
conda env export | grep -v prefix > tcga_cohort_conda_export.yml


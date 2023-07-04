# Running the RNA Pipeline

# Set Up

The Index for Star has to be pregenerated. 
For generating an index the build_STAR_index.sh skript in the scripts folder can be used.
Only generate a new one if the index for your reference vs. readlength combination is not available yet in this location.

# On LeoMed

## Load Environment
Using `source activate rna_pipeline` the standard environment with predefined software versions is used to run the pipeline.

If non-standard software versions are used:
- push the rna_environment.yaml file from the pipeline git into the project git.
- change the environment name inside the file and adapt the software versions or add additional software.
- Using an ssh tunnel (see wiki_page about leomed), create the environment using `conda env create -f rna_environment.yaml`
- Load the corresponding environment using `source activate ENVIRONMENT_NAME`


### Software not available in Conda
If additional software needs to be installed that is a non conda package they should go into the location:



### Scripts

- Currently we provide absolute path to the individual git if custom scripts are used. This needs to be changed as soon as possible. 
- All future Python Scripts should be written with Version 3+. 

## Dry Run

With a preloaded environment a dry run can be performed using:
```
snakemake -s example.snake --configfile config.json -n -p
```

## Run Pipeline
```
bsub -o **SAMPLE_NAME**.out -e **SAMPLE_NAME**.err -n 1 -W 3000 -R "rusage[mem=250]" -J MTBZ snakemake --notemp -s **example.snake** --configfile **config.json** --cluster "bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}" -j 100 -p -k
```

# On Euler

See information about "Running the RNA Pipeline" on LeoMed.

## R
Currently the R packages are installed by each user and stored in user specific directories. Necessary packages have therefore to be installed by each user separately.
R packages necessary for the ICGC pipeline:
-QoRTs



# Running the tcga cohort comparison

## On Euler

## Load environment
Using `source activate tcga_cohort` the standard conda environment with predefined software versions is used to run the pipeline.
In the config file it is then not necessary to provide an absolute path to the tool as all software installed in the conda environment is added to the PATH.
If the environment tcga_cohort is not available but conda is installed the environment can be generated using the text file tcga_cohort_conda_env.txt by executing 'conda create --name tcga_cohort --file tcga_cohort_conda_env.txt'


# Running the icgc pipeline for the Tumor Profiler project

## On Euler

## Load environment
Using `source activate icgc_pipeline` the standard conda environment with predefined software versions is used to run the pipeline. This includes also htseq-count (version 0.9.1) that is not included in the rna_pipeline conda environment.

For running the ICGC pipeline on one single sample (as it is done in the TP project) it is required to generate a local project copy of the icgc_snake.py file that is in the snakefiles/ directory. The "shell:" output of the "rule rnaseq" must be replaced with:

'python2.7 {params.generate_readStats} {params.indir} {params.samples}; date > {output}'

The rest of the code compares different samples of the same experiment to check for outliers or suboptimal samples. This is not applicable for an experiment using one sample.

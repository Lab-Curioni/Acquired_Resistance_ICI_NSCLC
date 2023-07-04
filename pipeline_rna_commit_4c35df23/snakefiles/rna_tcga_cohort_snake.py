import os, glob, sys

include: "./misc_snake_neu.py"

# This is the default order in which the programs are executed
# If the user specified a different order the user specified version is chosen.
FASTQIN = FASTQDIR
STARALIGNOUT = OUTDIR + 'star/'
RSEMCALCEXPRIN = STARALIGNOUT
RSEMCALCEXPROUT = OUTDIR + 'rsem/'
STRIPTABSIN = RSEMCALCEXPROUT
STRIPTABSOUT = RSEMCALCEXPROUT
PRUNEISOIN = RSEMCALCEXPROUT
PRUNEISOOUT = RSEMCALCEXPROUT
NORMQUANTIN = RSEMCALCEXPROUT
NORMQUANTOUT = RSEMCALCEXPROUT
CHANGEHEADIN = RSEMCALCEXPROUT
CHANGEHEADOUT = RSEMCALCEXPROUT
PARSECOHORTIN = RSEMCALCEXPROUT
PARSECOHORTOUT = OUTDIR + 'output_tcga/'
CLEANTABLEIN = PARSECOHORTOUT
CLEANTABLEOUT = OUTDIR + 'cleaned/'
ADDDGIDBIN = CLEANTABLEOUT
ADDDGIDBOUT = OUTDIR + 'final_table/'
FILTERGENESIN = ADDDGIDBOUT
FILTERGENESOUT = ADDDGIDBOUT
BOXPLOTIN = CHANGEHEADOUT
BOXPLOTOUT = OUTDIR + 'plots/'
DGIDBIN = CLEANTABLEOUT
DGIDBOUT = OUTDIR + 'dgidb/'

# Check if the user specified the proper input and output directories
if not 'FASTQDIR' in globals():
    print('You have to specify the root directory of the fastq files!')
    sys.exit(1)
if not 'OUTDIR' in globals():
    print('You have to specify the root directory where the results will be generated!')
    sys.exit(1)
if not 'PAIREDEND' in globals():
    print('You have to specify if you have Pairedend-data or not!')
    sys.exit(1)

# Include the rules
include: "./align.py"
include: "./rna_tcga_cohort_rules.py"
#include: "../snakefiles/stats.py"
#include: "../snakefiles/warnings.py"

localrules: tcga_rnaseq
rule tcga_rnaseq:
    input:
       expand(STARALIGNOUT + '{sample}_Aligned.toTranscriptome.out.bam',sample = getSampleNames()),
       expand(RSEMCALCEXPROUT + '{sample}.isoforms.results', sample = getSampleNames()),
       expand(STRIPTABSOUT + '{sample}.isoforms.results.bak', sample = getSampleNames()),
       expand(PRUNEISOOUT + '{sample}.genes.results.pruned.txt', sample = getSampleNames()),
       expand(NORMQUANTOUT + '{sample}.isoforms.results.normalized.txt', sample = getSampleNames()),
       expand(NORMQUANTOUT + '{sample}.genes.results.pruned.normalized.txt', sample = getSampleNames()),
       expand(CHANGEHEADOUT + '{sample}.genes.results.pruned.normalized.header.txt', sample = getSampleNames()),
       #expand(PARSECOHORTOUT + '{sample}.final_rna_table.txt', sample = getSampleNames()),
       expand(CLEANTABLEOUT + '{sample}.final_rna_table.txt', sample = getSampleNames()),
       expand(BOXPLOTOUT + '{sample}_boxplots_success.txt', sample = getSampleNames()),
       expand(DGIDBOUT + '{sample}.dgidb.txt', sample = getSampleNames()),
       expand(ADDDGIDBOUT + '{sample}.final_rna_table_dgidb_results.txt', sample = getSampleNames()),
       expand(ADDDGIDBOUT + '{sample}.final_rna_table_dgidb_results_filtered_clinical_genes.txt', sample = getSampleNames())
    output:
        OUTDIR + 'rna_tcga_complete.txt'
    params:
        lsfoutfile = OUTDIR + 'complete.lsfout.log',
        lsferrfile = OUTDIR + 'complete.lsferr.log'
    shell:
        'touch {output}'

import os, glob, sys

# This function adapts the config object to include full path information
include: "../common/misc/misc_snake_germline.py"
postProcessConfigMap()

# Check if the uses specified the proper input and output directories
if not 'FASTQDIR' in globals():
    print('You have to specify the root directory of the fastq files!')
    sys.exit(1)
if not 'OUTDIR' in globals():
    print('You have to specify the root directory where the results will be generated!')
    sys.exit(1)
if not 'TMPDIR' in globals():
    print('You have to specify the root directory where temporary files will be stored!')
    sys.exit(1)
# This is the default order in which the programs are executed
# If the user specified a different order the user specified version is chose
if not 'SEQPURGEIN' in globals():
    SEQPURGEIN = FASTQDIR
if not 'CLIPPEDOUT' in globals():
    CLIPPEDOUT = OUTDIR + 'clipped/'
if not 'BWAIN' in globals():
    BWAIN = CLIPPEDOUT
if not 'BWAOUT' in globals():
    BWAOUT = OUTDIR + 'bwa/'
if not 'FIXMATEANDSORTIN' in globals():
    FIXMATEANDSORTIN = BWAOUT
if not 'FIXMATEANDSORTOUT' in globals():
    FIXMATEANDSORTOUT = OUTDIR + 'fixmatepair/'
if not 'MERGEBAMSIN' in globals():
    MERGEBAMSIN = FIXMATEANDSORTOUT
if not 'MERGEBAMSOUT' in globals():
    MERGEBAMSOUT = OUTDIR + 'mergebams/'
if not 'NOSECONDARYALNIN' in globals():
    NOSECONDARYALNIN = MERGEBAMSOUT
if not 'NOSECONDARYALNOUT' in globals():
    NOSECONDARYALNOUT = OUTDIR + 'noSecondaryAln/'
if not 'MARKPCRDUBLICATESIN' in globals():
    MARKPCRDUBLICATESIN = NOSECONDARYALNOUT
if not 'MARKPCRDUBLICATESOUT' in globals():
    MARKPCRDUBLICATESOUT = OUTDIR + 'markDuplicates/'
if not 'REMOVEDUBLICATESIN' in globals():
    REMOVEPCRDUBLICATESIN = MARKPCRDUBLICATESOUT
if not 'REMOVEPCRDUBLICATESOUT' in globals():
    REMOVEPCRDUBLICATESOUT = OUTDIR + 'removedDuplicates/'
if not 'NAMESORTEDIN' in globals():
    NAMESORTEDIN = REMOVEPCRDUBLICATESOUT
if not 'NAMESORTEDOUT' in globals():
    NAMESORTEDOUT = OUTDIR + 'namesorted/'
if not 'CLIPOVERLAPIN' in globals():
    CLIPOVERLAPIN = NAMESORTEDOUT
if not 'CLIPOVERLAPOUT' in globals():
    CLIPOVERLAPOUT = OUTDIR + 'clippedOverlaps/'
if not 'SORTEDIN' in globals():
    SORTEDIN = CLIPOVERLAPOUT
if not 'SORTEDOUT' in globals():
    SORTEDOUT = OUTDIR + 'sorted/'
if not 'BASERECALIBRATIONIN' in globals():
    BASERECALIBRATIONIN = SORTEDOUT
if not 'BASERECALIBRATIONOUT' in globals():
    BASERECALIBRATIONOUT = OUTDIR + 'recalibratedBases/'
if not 'HAPLOTYPECALLERIN'  in globals():
    HAPLOTYPECALLERIN = BASERECALIBRATIONOUT
if not 'HAPLOTYPECALLEROUT' in globals():
    HAPLOTYPECALLEROUT = OUTDIR + 'haplotypeCaller/'
if not 'GENOTYPERIN' in globals():
    GENOTYPERIN = HAPLOTYPECALLEROUT
if not '' in globals():
    GENOTYPEROUT = OUTDIR + 'variants/'
# Definition of some constantly used lists
SINGLEFASTQFILES = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq.gz')]
if not SINGLEFASTQFILES :
    SINGLEFASTQFILES = [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq')]
PAIREDFASTQFILES = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*.fastq.gz')]
if not PAIREDFASTQFILES:
    PAIREDFASTQFILES = [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*.fastq')]
PAIREDFASTQFILESWITHOUTR = [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
if not PAIREDFASTQFILESWITHOUTR:
    PAIREDFASTQFILESWITHOUTR = [file.replace(FASTQDIR, '').replace('_R1.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq')]
PAIREDUNPAIRED = ['_PAIRED','_UNPAIRED']
MATEPAIR = ['_R1', '_R2']
SAMPLENAMES=getSampleNames()
EXPERIMENTS=getExperiments()
PAIREDFASTQFILESWITHDIR = [FASTQDIR + x for x in PAIREDFASTQFILES]
# Create all necessary directories
# The creation is necesseccary because the otherwise the log files cannot be written on some cluster systems
def createDirs():
    for d in CLIPPEDOUT, BWAOUT:
        for s in SAMPLENAMES:
            if not os.path.exists(d + '/' + s + '/PAIREDEND/'):
                os.makedirs(d + '/' + s + '/PAIREDEND/')


# Include the rules
include: "../common/seqpurge/seq_purge_snake.py"
include: "../common/align/align_snake.py"
include: "../common/bam_mod/bam_mod.py"
include: "../common/variants/call_variants.py"
include: "../common/variants/variant_mod.py"
include: "../common/stats/stats_snake.py"
# This rule defines which files should be created
rule wes:
	input:
		expand(MERGEBAMSOUT + '{sample}.bam', sample = SAMPLENAMES),
		expand(MERGEBAMSOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
		expand(MERGEBAMSOUT + '{sample}.hsmetrics.txt', sample = SAMPLENAMES),
		expand(NOSECONDARYALNOUT + '{sample}.bam', sample = SAMPLENAMES),
		expand(MARKPCRDUBLICATESOUT + '{sample}.bam', sample = SAMPLENAMES), 
		expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam', sample = SAMPLENAMES),
		expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
		expand(REMOVEPCRDUBLICATESOUT + '{sample}.hsmetrics.txt', sample = SAMPLENAMES),
		expand(REMOVEPCRDUBLICATESOUT + '{sample}.bai', sample = SAMPLENAMES),
		expand(NAMESORTEDOUT + '{sample}.bam', sample = SAMPLENAMES),
		expand(CLIPOVERLAPOUT + '{sample}.bam', sample = SAMPLENAMES),
		expand(CLIPOVERLAPOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
                #expand(CLIPOVERLAPOUT + '{sample}.hsmetrics.txt', sample = SAMPLENAMES),
		expand(SORTEDOUT + '{sample}.bam', sample = SAMPLENAMES),
		expand(BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table', sample = SAMPLENAMES),
		expand(BASERECALIBRATIONOUT + '{sample}.bam', sample = SAMPLENAMES),
		expand(HAPLOTYPECALLEROUT + '{sample}.g.vcf', sample = SAMPLENAMES),
		expand(GENOTYPEROUT + '{experiment}.vcf', experiment = EXPERIMENTS),
		expand(GENOTYPEROUT + '{experiment}.snps.recal', experiment = EXPERIMENTS),
		expand(GENOTYPEROUT + '{experiment}.indels.recal', experiment = EXPERIMENTS),
		expand(GENOTYPEROUT + '{experiment}.snps.recal.vcf', experiment = EXPERIMENTS),
		expand(GENOTYPEROUT + '{experiment}.recal.vcf', experiment = EXPERIMENTS),
		expand(GENOTYPEROUT + '{experiment}.final.vcf', experiment = EXPERIMENTS)
	output:
        	OUTDIR + 'complete.txt'
	params:
	        lsfoutfile = OUTDIR + 'complete.lsfout.log',
       		lsferrfile = OUTDIR + 'complete.lsferr.log'
	shell:
        	'date > {output}'

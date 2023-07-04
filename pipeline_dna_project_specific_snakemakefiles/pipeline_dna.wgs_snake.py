import os, glob, sys

# This function adapts the config object to include full path information
include: "../common/misc/misc_snake.py"
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
if not 'WHICHCLIP' in globals():
    print('You have to specify if you want to use Cliptrim or Seqpurge!')
    sys.exit(1)

print(WHICHCLIP)
# This is the default order in which the programs are executed
# If the user specified a different order the user specified version is chosen.
if WHICHCLIP is True:
        if not 'SEQPURGEIN' in globals():
                SEQPURGEIN = FASTQDIR
else:
        if not 'CLIPTRIMIN' in globals():
                CLIPTRIMIN = FASTQDIR
if not 'CLIPPEDOUT' in globals():
    CLIPPEDOUT = OUTDIR + 'clipped/'
if not 'BWAIN' in globals():
    BWAIN = CLIPPEDOUT
if not 'BWAOUT' in globals():
    BWAOUT = OUTDIR + 'bwa/'
if not 'FIXMATEANDSORTIN' in globals():
    FIXMATEANDSORTIN = BWAOUT
if not 'FIXMATEANDSORTOUT' in globals():
    FIXMATEANDSORTOUT = OUTDIR + 'fix_sorted/'
if not 'MERGEBAMSIN' in globals():
    MERGEBAMSIN = FIXMATEANDSORTOUT
if not 'MERGEBAMSOUT' in globals():
    MERGEBAMSOUT = OUTDIR + 'merged/'
if not 'NOSECONDARYALNIN' in globals():
    NOSECONDARYALNIN = MERGEBAMSOUT
if not 'NOSECONDARYALNOUT' in globals():
    NOSECONDARYALNOUT = OUTDIR + 'noSecondaryAln/'
if not 'MARKPCRDUBLICATESIN' in globals():
    MARKPCRDUBLICATESIN = NOSECONDARYALNOUT
if not 'MARKPCRDUBLICATESOUT' in globals():
    MARKPCRDUBLICATESOUT = OUTDIR + 'markedDuplicates/'
if not 'REMOVEPCRDUBLICATESIN' in globals():
    REMOVEPCRDUBLICATESIN = MARKPCRDUBLICATESOUT
if not 'REMOVEPCRDUBLICATESOUT' in globals():
    REMOVEPCRDUBLICATESOUT = OUTDIR + 'removedPcrDuplicates/'
if not 'REALIGNINDELSIN' in globals():
    REALIGNINDELSIN = REMOVEPCRDUBLICATESOUT
if not 'REALIGNINDELSOUT' in globals():
    REALIGNINDELSOUT = OUTDIR + 'realignedIndels/'
if not 'BASERECALIBRATIONIN' in globals():
    BASERECALIBRATIONIN = REALIGNINDELSOUT
if not 'BASERECALIBRATIONOUT' in globals():
    BASERECALIBRATIONOUT = OUTDIR + 'recalibratedBases/'
if not 'BICSEQ2IN' in globals():
    BICSEQ2IN = REMOVEPCRDUBLICATESOUT
if not 'BICSEQ2OUT' in globals():
    BICSEQ2OUT = OUTDIR + 'bicseq2/'
if not 'MPILEUPIN' in globals():
    MPILEUPIN = REMOVEPCRDUBLICATESOUT
if not 'MPILEUPOUT' in globals():
    MPILEUPOUT = OUTDIR + 'mpileup/'
if not 'VARSCANCNVIN' in globals():
    VARSCANCNVIN = MPILEUPOUT
if not 'VARSCANCNVOUT' in globals():
    VARSCANCNVOUT = OUTDIR + 'varscan_cnv/'
if not 'EXCAVATORIN' in globals():
    EXCAVATORIN = BASERECALIBRATIONOUT
if not 'EXCAVATOROUT' in globals():
    EXCAVATOROUT = OUTDIR + 'copynumber/excavator/'
if not 'DATABASEQUERY' in globals():
    DATABASEQUERY = OUTDIR + 'databaseQuery/'
if not 'DOWNLOADCLINICALTRIALSOUT' in globals():
    DOWNLOADCLINICALTRIALSOUT = OUTDIR + 'clinicalTrials/'

# Definition of some constantly used lists
# TODO(jochen): Add lists for BAM files (maybe even for each level - align, sort, merge, ...)
SINGLEFASTQFILES = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq.gz')]
PAIREDFASTQFILES = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*.fastq.gz')]
PAIREDFASTQFILESORPAHNS = [file.replace(FASTQDIR, '').replace('PAIREDEND','PAIREDEND/ORPHAN').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*.fastq.gz')]
PAIREDFASTQFILESWITHOUTR = [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
PAIREDUNPAIRED = ['_PAIRED','_UNPAIRED']
MATEPAIR = ['_R1', '_R2']
SAMPLENAMES=getSampleNames()
print(PAIREDFASTQFILESORPAHNS)

# Create all necessary directories
# The creation is necesseccary because the otherwise the log files cannot be written on some cluster systems
def createDirs():
    for d in CLIPPEDOUT, BWAOUT, FIXMATEANDSORTOUT:
        for s in SAMPLENAMES:
            if not os.path.exists(d + '/' + s + '/PAIREDEND/'):
                os.makedirs(d + '/' + s + '/PAIREDEND/')
    for d in MERGEBAMSOUT, NOSECONDARYALNOUT, MARKPCRDUBLICATESOUT, REMOVEPCRDUBLICATESOUT, REALIGNINDELSOUT, BASERECALIBRATIONOUT, FREEBAYESOUT, HAPLOTYPECALLEROUT:
        if not os.path.exists(d):
            os.makedirs(d)
#createDirs()

# Here are some directories for rules which are currently not used!
if not 'BWAALNIN' in globals():
    BWAALNIN = CLIPPEDOUT
if not 'BWAALNOUT' in globals():
    BWAALNOUT = OUTDIR + 'bwa_aln/'
if not 'BOWTIEIN' in globals():
    BOWTIEIN = CLIPPEDOUT
if not 'BOWTIEOUT' in globals():
    BOWTIEOUT = OUTDIR + 'bowtie2_out/'
if not 'YARAIN' in globals():
    YARAIN = OUTDIR + '.yara_in'
if not 'YARAOUT' in globals():
    YARAOUT = OUTDIR + '.yara_out'
if not 'REASSIGNONEMAPPINGQUALIN' in globals():
    REASSIGNONEMAPPINGQUALIN = OUTDIR + '.reassing_one_mapping_quality_in'
if not 'REASSIGNONEMAPPINGQUALOUT' in globals():
    REASSIGNONEMAPPINGQUALOUT = OUTDIR + '.reassing_one_mapping_quality_out'
if not 'SOAPIN' in globals():
    SOAPIN = OUTDIR + '.yara_in'
if not 'SOAPOUT' in globals():
    SOAPOUT = OUTDIR + '.soap_out'
if not 'SL_DIR' in globals():
    SL_DIR = OUTDIR + 'spezialitaetenliste_database/'

# Include the rules
if WHICHCLIP is True:
        include: "../common/seqpurge/seq_purge_snake.py"
else:
        include: "../common/clip_trim/clip_trim_snake.py"
include: "../common/align/align_snake.py"
include: "../common/bam_mod/wgs_bam_mod_snake.py"
include: "../common/stats/stats_snake.py"
include: "../common/copy_number/copy_number_snake.py"
include: "../common/database_query/database_query_snake.py"
#include: "../common/variants/variants_snake.py"
#include: "../common/variants/variant_mod_snake.py"

#print(expand(VARSCANCOMPLETEOUT + '{tumorNormalMatching}.vcf', tumorNormalMatching = getNormalTumorFiles()))
# This rule defines which files should be created
rule wgs:
    input:    
        #expand(MERGEBAMSOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(MERGEBAMSOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
	expand(MERGEBAMSOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
	#expand(NOSECONDARYALNOUT + '{sample}.bam', sample = SAMPLENAMES),
	expand(NOSECONDARYALNOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
	expand(NOSECONDARYALNOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
	#expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam', sample = SAMPLENAMES),
	expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
	expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
	expand(BICSEQ2OUT + '{sample}/{contigNames}.seq', sample = SAMPLENAMES, contigNames = getContigNames()),
	os.path.realpath(config['resources'][ORGANISM]['reference']) + '_contigs/chr1.fasta',
	expand(BICSEQ2OUT + '{sample}/configNorm.txt', sample = SAMPLENAMES),
	expand(BICSEQ2OUT + '{sample}/paramsEstimate.txt', sample = SAMPLENAMES),
	expand(BICSEQ2OUT + '{tumorNormalMatching}/configSeg.txt', tumorNormalMatching = getNormalTumorFiles()),
	expand(BICSEQ2OUT + '{tumorNormalMatching}.cnvsRaw.txt', tumorNormalMatching = getNormalTumorFiles()),
	expand(BICSEQ2OUT + '{tumorNormalMatching}.cnvsGenotype.txt', tumorNormalMatching = getNormalTumorFiles()),
	expand(BICSEQ2OUT + '{tumorNormalMatching}.filtered.txt', tumorNormalMatching = getNormalTumorFiles()),
	expand(DATABASEQUERY + '{tumorNormalMatching}.filtered.CNVannotated.overview.txt', tumorNormalMatching = getNormalTumorFiles()),
	expand(DATABASEQUERY + '{tumorNormalMatching}.filtered.CNVannotated.overview.dgidb.txt.CompleteTable.txt', tumorNormalMatching = getNormalTumorFiles()),
	expand(DATABASEQUERY + '{tumorNormalMatching}.filtered.CNVannotated.overview.dgidb.txt.CompleteTable.ClinicalTrials.txt', tumorNormalMatching = getNormalTumorFiles()),
	expand(DATABASEQUERY + '{tumorNormalMatching}.filtered.CNVannotated.overview.databaseQueriesCNV.txt', tumorNormalMatching = getNormalTumorFiles()),
	#expand(VARSCANCNVOUT + '{tumorNormalMatching}.copynumber', tumorNormalMatching = getNormalTumorFiles()),
    output:
        OUTDIR + 'complete.txt'
    params:
        lsfoutfile = OUTDIR + 'complete.lsfout.log',
        lsferrfile = OUTDIR + 'complete.lsferr.log',
        #scratch = config['tools']['wes']['scratch'],
        #mem = config['tools']['wes']['mem'],
        #time = config['tools']['wes']['time']
    shell:
        'date > {output}'


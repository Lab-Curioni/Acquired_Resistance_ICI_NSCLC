import os, glob, sys
from itertools import chain

# functionality like getSampleNames
include: "../NGS-pipe/snake/common/misc/misc_snake.py"

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
if not 'SEQPURGEOUT' in globals():
    SEQPURGEOUT = OUTDIR + 'seqpurge/'
if not 'CLIPTRIMIN' in globals():
    CLIPTRIMIN = FASTQDIR
if not 'CLIPPEDOUT' in globals():
    CLIPPEDOUT = OUTDIR + 'clipped/'
if not 'HLATYPEIN' in globals():
    HLATYPEIN = FASTQDIR
if not 'HLATYPEOUT' in globals():
    HLATYPEOUT = OUTDIR + 'hla_type/'

if not 'BWAIN' in globals():
    BWAIN = SEQPURGEOUT
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
# conpair is a hg19-only tool!
if not 'CONPAIRIN' in globals():
    CONPAIRIN = MERGEBAMSOUT
if not 'CONPAIROUT' in globals():
    CONPAIROUT = OUTDIR + 'conpair/'
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
if not 'NAMESORTEDIN' in globals():
    NAMESORTEDIN = REALIGNINDELSOUT
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
if not 'FREEBAYESIN' in globals():
    FREEBAYESIN = REALIGNINDELSOUT
if not 'FREEBAYESOUT' in globals():
    FREEBAYESOUT = OUTDIR + 'variants/freebayes/raw/'
if not 'MPILEUPIN' in globals():
    MPILEUPIN = BASERECALIBRATIONOUT
if not 'MPILEUPOUT' in globals():
    MPILEUPOUT = OUTDIR + 'mpileup/'
if not 'VARSCANSNPIN' in globals():
    VARSCANSNPIN = MPILEUPOUT
if not 'VARSCANSNPOUT' in globals():
    VARSCANSNPOUT = OUTDIR + 'variants/varscan_snp/raw/'
if not 'VARSCANSOMATICIN' in globals():
    VARSCANSOMATICIN = MPILEUPOUT
if not 'VARSCANSOMATICOUT' in globals():
    VARSCANSOMATICOUT = OUTDIR + 'variants/varscan_somatic/raw/'
if not 'VARSCANUPDATEHEADERIN' in globals():
    VARSCANUPDATEHEADERIN = VARSCANSOMATICOUT
if not 'VARSCANUPDATEHEADEROUT' in globals():
    VARSCANUPDATEHEADEROUT = OUTDIR + 'variants/varscan_somatic/complete_raw/'
if not 'VARSCANCOMPLETEIN' in globals():
    VARSCANCOMPLETEIN = VARSCANUPDATEHEADEROUT
if not 'VARSCANCOMPLETEOUT' in globals():
    VARSCANCOMPLETEOUT = OUTDIR + 'variants/varscan_somatic/combined_raw/'
if not 'VARSCANSOMATICFILTERIN' in globals():
    VARSCANSOMATICFILTERIN = VARSCANCOMPLETEOUT
if not 'VARSCANSOMATICFILTEROUT' in globals():
    VARSCANSOMATICFILTEROUT = OUTDIR + 'variants/varscan_somatic/filtered/'
if not 'BCFTOOLSIN' in globals():
    BCFTOOLSIN = MPILEUPOUT
if not 'BCFTOOLSOUT' in globals():
    BCFTOOLSOUT = OUTDIR + 'bcftools/raw/'
if not 'HAPLOTYPECALLERIN' in globals():
    HAPLOTYPECALLERIN = BASERECALIBRATIONOUT
if not 'HAPLOTYPECALLEROUT' in globals():
    HAPLOTYPECALLEROUT = OUTDIR + 'variants/GATK/raw/'
if not 'MUTECT2IN' in globals():
    MUTECT2IN = BASERECALIBRATIONOUT
if not 'MUTECT2OUT' in globals():
    MUTECT2OUT = OUTDIR + 'variants/mutect2/'
if not 'MUTECT1IN' in globals():
    MUTECT1IN = BASERECALIBRATIONOUT
if not 'MUTECT1OUT' in globals():
    MUTECT1OUT = OUTDIR + 'variants/mutect1/raw/'
if not 'MUTECT1FILTERIN' in globals():
    MUTECT1FILTERIN = MUTECT1OUT
if not 'MUTECT1FILTEROUT' in globals():
    MUTECT1FILTEROUT = OUTDIR + 'variants/mutect1/filtered/'
if not 'JOINTSNVMIX2_075_IN' in globals():
    JOINTSNVMIX2_075_IN = BASERECALIBRATIONOUT
if not 'JOINTSNVMIX2_075_OUT' in globals():
    JOINTSNVMIX2_075_OUT = OUTDIR + 'variants/jointSNVMix2_075/raw/'
if not 'JOINTSNVMIX2IN' in globals():
    JOINTSNVMIX2IN = BASERECALIBRATIONOUT
if not 'JOINTSNVMIX2OUT' in globals():
    JOINTSNVMIX2OUT = OUTDIR + 'variants/jointSNVMix2/raw/'
if not 'SOMATICSNIPERIN' in globals():
    SOMATICSNIPERIN = BASERECALIBRATIONOUT
if not 'SOMATICSNIPEROUT' in globals():
    SOMATICSNIPEROUT = OUTDIR + 'variants/somaticSniper/raw/'
if not 'VARDICTIN' in globals():
    VARDICTIN = BASERECALIBRATIONOUT
if not 'VARDICTOUT' in globals():
    VARDICTOUT = OUTDIR + 'variants/varDict/raw/'
if not 'SOMATICSEQOUT' in globals():
    SOMATICSEQOUT = OUTDIR + 'variants/somatic_seq/'
if not 'STRELKAIN' in globals():
    STRELKAIN = BASERECALIBRATIONOUT
if not 'STRELKAOUT' in globals():
    STRELKAOUT = OUTDIR + 'variants/strelka/'
if not 'STRELKAFILTERIN' in globals():
    STRELKAFILTERIN = STRELKAOUT
if not 'STRELKAFILTEROUT' in globals():
    STRELKAFILTEROUT = OUTDIR + 'variants/strelka/filtered/'
if not 'GATKVARIANTCOMBINEOUT' in globals():
    GATKVARIANTCOMBINEOUT = OUTDIR + 'variants/combined/'
if not 'PARSEANNOTATEDVCFIN' in globals():
    PARSEANNOTATEDVCFIN = GATKVARIANTCOMBINEOUT
if not 'PARSEANNOTATEDVCFOUT' in globals():
    PARSEANNOTATEDVCFOUT = OUTDIR + 'databaseQuery/'
#if not 'DATABASEQUERY' in globals():
#    DATABASEQUERY = OUTDIR + 'databaseQuery/'
if not 'EXCAVATORIN' in globals():
    EXCAVATORIN = BASERECALIBRATIONOUT
if not 'EXCAVATOROUT' in globals():
    EXCAVATOROUT = OUTDIR + 'copynumber/excavator/'
if not 'EXCAVATOR2IN' in globals():
    EXCAVATOR2IN = BASERECALIBRATIONOUT
if not 'EXCAVATOR2OUT' in globals():
    EXCAVATOR2OUT = OUTDIR + 'copynumber/excavator2/'
if not 'VARSCANCNVIN' in globals():
    VARSCANCNVIN = MPILEUPOUT
if not 'VARSCANCNVOUT' in globals():
    VARSCANCNVOUT = OUTDIR + 'varscan_cnv/'
if not 'BICSEQ2IN' in globals(): 
    BICSEQ2IN = REMOVEPCRDUBLICATESOUT
if not 'BICSEQ2OUT' in globals():
    BICSEQ2OUT = OUTDIR + 'bicseq2/'
if not 'DGIDB_IN' in globals():
    DGIDB_IN = PARSEANNOTATEDVCFOUT	
if not 'DGIDB_OUT' in globals():
    DGIDB_OUT = OUTDIR + 'dgidb/'		
if not 'CLINICALTRIALS_IN' in globals():
    CLINICALTRIALS_IN = DGIDB_OUT	
if not 'CLINICALTRIALS_OUT' in globals():
    CLINICALTRIALS_OUT = OUTDIR + 'clinicalTrials/'		
if not 'SL_DIR' in globals():
    SL_DIR = OUTDIR + 'spezialitaetenliste_database/'
if not 'CREATEREFERENCEHEADERIN' in globals():
    CREATEREFERENCEHEADERIN = REMOVEPCRDUBLICATESOUT
if not 'CREATEREFERENCEHEADEROUT' in globals():
    CREATEREFERENCEHEADEROUT = OUTDIR + 'variants/'
if not 'SOMATICSIGNATURES_IN' in globals():
    SOMATICSIGNATURES_IN = GATKVARIANTCOMBINEOUT
if not 'SOMATICSIGNATURES_OUT' in globals():
    SOMATICSIGNATURES_OUT = OUTDIR + 'variants/somaticSignatures/'
if not 'ANNOTATECLINICAL_IN' in globals():
    ANNOTATECLINICAL_IN = PARSEANNOTATEDVCFOUT	
if not 'ANNOTATECLINICAL_OUT' in globals():
    ANNOTATECLINICAL_OUT = OUTDIR + 'clinicalAnnotation/'		
if not 'MUTATIONALBURDEN_OUT' in globals():
    MUTATIONALBURDEN_OUT = OUTDIR + 'variants/mutationalBurden/'

# Definition of some constantly used lists

SAMPLENAMES=getSampleNames()
SINGLEFASTQFILES=getSingleFastqFiles(SAMPLENAMES)
PAIREDFASTQFILES=getPairedFastqFiles(SAMPLENAMES)
PAIREDFASTQFILESWITHOUTR=getPairedFastqFilesWithoutR(SAMPLENAMES)

print(SAMPLENAMES)

# Include the rules
include: "../NGS-pipe/snake/common/clip_trim/clip_trim_snake.py"
include: "../NGS-pipe/snake/common/align/align_snake.py"
include: "../NGS-pipe/snake/common/bam_mod/bam_mod_snake.py"
include: "../common/bam_mod/bam_mod_snake.py"
include: "../NGS-pipe/snake/common/stats/stats_snake.py"
include: "../NGS-pipe/snake/common/copy_number/copy_number_snake.py"
include: "../common/copy_number/copy_number_snake.py"
include: "../NGS-pipe/snake/common/variants/variants_snake.py"
#include: "../../NGS-pipe/snake/common/variants/somatic_seq_snake.py"
include: "../NGS-pipe/snake/common/variants/variant_mod_snake.py"
include: "../common/variants/variant_mod_snake.py"

# This rule defines which files should be created
rule wes:
    input:
        [file.replace('.fastq.gz', '_fastqc.html') for file in list(chain.from_iterable(glob.glob(FASTQDIR + SAMPLE + '/PAIREDEND/*.fastq.gz') for SAMPLE in SAMPLENAMES))],
        [file.replace('.fastq.gz', '_fastqc.html').replace(FASTQDIR, SEQPURGEOUT) for file in list(chain.from_iterable(glob.glob(FASTQDIR + SAMPLE + '/PAIREDEND/*.fastq.gz') for SAMPLE in SAMPLENAMES))],
        [file.replace('.fastq.gz', '.count') for file in list(chain.from_iterable(glob.glob(FASTQDIR + SAMPLE + '/PAIREDEND/*.fastq.gz') for SAMPLE in SAMPLENAMES))],
        [file.replace('.fastq.gz', '.count').replace(FASTQDIR, SEQPURGEOUT) for file in list(chain.from_iterable(glob.glob(FASTQDIR + SAMPLE + '/PAIREDEND/*.fastq.gz') for SAMPLE in SAMPLENAMES))],
        expand(MERGEBAMSOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(MERGEBAMSOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(MERGEBAMSOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(NOSECONDARYALNOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(NOSECONDARYALNOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(NOSECONDARYALNOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(REMOVEPCRDUBLICATESOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(REALIGNINDELSOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(REALIGNINDELSOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(REALIGNINDELSOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(NAMESORTEDOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(CLIPOVERLAPOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(SORTEDOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(BASERECALIBRATIONOUT + '{sample}.bam', sample = SAMPLENAMES),
        expand(BASERECALIBRATIONOUT + '{sample}.bam.flagstat', sample = SAMPLENAMES),
        expand(BASERECALIBRATIONOUT + '{sample}.bam_stats/report.pdf', sample = SAMPLENAMES),
        expand(BASERECALIBRATIONOUT + '{sample}_base_recalibration_report.pdf', sample = SAMPLENAMES),
        expand(VARSCANCOMPLETEOUT + '{tumorNormalMatching}.vcf', tumorNormalMatching = getNormalTumorFiles()),
        expand(CREATEREFERENCEHEADEROUT + '{tumorNormalMatching}.referenceNames_forVCFheaderUpdate.txt', tumorNormalMatching = getNormalTumorFiles()),
        expand(MUTECT1FILTEROUT + '{tumorNormalMatching}.correctNames.snpEff.dbSNP.pass.vcf', tumorNormalMatching = getNormalTumorFiles()),
        expand(STRELKA2FILTEROUT + '{tumorNormalMatching}.correctNames.snpEff.dbSNP.cosmic.pass.vcf', tumorNormalMatching = getNormalTumorFiles()),
        expand(VARSCANSOMATICFILTEROUT + '{tumorNormalMatching}.correctNames.snpEff.dbSNP.cosmic.pass.vcf', tumorNormalMatching = getNormalTumorFiles()),
        expand(FREEBAYESFILTEROUT + '{tumorNormalMatching}.correctNames.snpEff.dbSNP.cosmic.pass.vcf', tumorNormalMatching = getNormalTumorFiles()),
        expand(GATKVARIANTCOMBINEOUT + '{tumorNormalMatching}.correctNames.snpEff.dbSNP.cosmic.pass.combined.vcf', tumorNormalMatching = getNormalTumorFiles()),
        expand(SOMATICSIGNATURES_OUT + '{tumorNormalMatching}.correctNames.snpEff.dbSNP.cosmic.pass.combined.somaticSignatures.txt', tumorNormalMatching = getNormalTumorFiles()),
        #expand(ANNOTATECLINICAL_OUT + '{tumorNormalMatching}.correctNames.snpEff.dbSNP.cosmic.pass.combined.clinicalAnnotation.txt', tumorNormalMatching = getNormalTumorFiles()),
        #expand(ANNOTATECLINICAL_OUT + '{tumorNormalMatching}.correctNames.snpEff.dbSNP.cosmic.pass.combined.clinicalAnnotation.txt_dgidbIndependent.txt.damaging.txt', tumorNormalMatching = getNormalTumorFiles()),
       # expand(MUTATIONALBURDEN_OUT + '{tumorNormalMatching}.correctNames.snpEff.dbSNP.cosmic.pass.combined.clinicalAnnotation.txt_dgidbIndependent.txt.mutationalBurden.txt', tumorNormalMatching = getNormalTumorFiles())
        #expand(EXCAVATOR2OUT + '{tumorNormalMatching}/' + '{tumorNormalMatching}.copyNumSuccess.txt', tumorNormalMatching = getNormalTumorFiles()),
        #expand(EXCAVATOR2OUT + '{tumorNormalMatching}.excavator2CNV.txt', tumorNormalMatching = getNormalTumorFiles())
        expand(EXCAVATOROUT + '{tumorNormalMatching}/' + '{tumorNormalMatching}.copyNumSuccess.txt', tumorNormalMatching = getNormalTumorFiles()),
        expand(EXCAVATOROUT + '{tumorNormalMatching}.excavatorCNV.txt', tumorNormalMatching = getNormalTumorFiles()),
        expand(FREEBAYESOUT + 'all.vcf') 
    output:
        OUTDIR + 'complete.txt'
    params:
        lsfoutfile = OUTDIR + 'complete.lsfout.log',
        lsferrfile = OUTDIR + 'complete.lsferr.log',
        mem = '1000',
        scratch = '1000',
        time = '1'
    benchmark:
        OUTDIR + 'complete.txt.benchmark'
    shell:
        'date > {output}'

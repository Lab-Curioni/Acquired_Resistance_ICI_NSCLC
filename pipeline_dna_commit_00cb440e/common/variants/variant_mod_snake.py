'''
# This rule uses picard tools createVCFdictionary to fill the header of VarScan and freebayes vcf files
rule updateVCFSequenceDictionary:
    input:
        vcf = VARSCANUPDATEHEADERIN + '{sample}.vcf',
        reference = config['resources'][ORGANISM]['reference']
    output:
        vcf = VARSCANUPDATEHEADEROUT + '{sample}.vcf'
    params:
        lsfoutfile = VARSCANUPDATEHEADEROUT + '{sample}.vcf.lsfout.log',
        lsferrfile = VARSCANUPDATEHEADEROUT + '{sample}.vcf.lsferr.log',
        scratch = config['tools']['picard']['updateVCFSequenceDictionary']['scratch'],
        mem = config['tools']['picard']['updateVCFSequenceDictionary']['mem'],
        time = config['tools']['picard']['updateVCFSequenceDictionary']['time']
    threads:
        int(config['tools']['picard']['updateVCFSequenceDictionary']['threads'])
    benchmark:
        VARSCANUPDATEHEADEROUT + '{sample}.vcf.benchmark'
    shell:
        '{config[tools][picard][call]} UpdateVcfSequenceDictionary INPUT={input.vcf} SD={input.reference} OUTPUT={output.vcf}'
'''

if not 'VARSCANSNPFILTERIN' in globals():
    VARSCANSNPFILTERIN = VARSCANSNPOUT
if not 'VARSCANSNPFILTEROUT' in globals():
    VARSCANSNPFILTEROUT = OUTDIR + 'variants/varscan_snp/filtered/'
# This rule filters annotated vcf files produced by VarScan2 nonsomatic
rule filterVarScanNonsomatic:
    input:
        vcf = VARSCANSNPFILTERIN + '{sample}.vcf'
    output:
        vcf = VARSCANSNPFILTEROUT + '{sample}.SNPpass.vcf'
    params:
        lsfoutfile = VARSCANSNPFILTEROUT + '{sample}.SNPpass.vcf.lsfout.log',
        lsferrfile = VARSCANSNPFILTEROUT + '{sample}.SNPpass.vcf.lsferr.log',
        scratch = config['tools']['varscanNonSomaticFilter']['scratch'],
        mem = config['tools']['varscanNonSomaticFilter']['mem'],
        time = config['tools']['varscanNonSomaticFilter']['time'],
        minVarSupport = config['tools']['varscanNonSomaticFilter']['minVarSupport'],
        pvalue = config['tools']['varscanNonSomaticFilter']['pvalue'],
        minNucCoverage = config['tools']['varscanNonSomaticFilter']['minNucCoverage'],
        filterStrands = config['tools']['varscanNonSomaticFilter']['filterStrands'],
        tumorFreqThreshold = config['tools']['varscanNonSomaticFilter']['tumorFreqThreshold'],
        filterHomopolymer = config['tools']['varscanNonSomaticFilter']['filterHomopolymer'],
        filterSilent = config['tools']['varscanNonSomaticFilter']['filterSilent'],
        filterCommon = config['tools']['varscanNonSomaticFilter']['filterCommon']
    threads:
        config['tools']['varscanNonSomaticFilter']['threads']
    benchmark:
        VARSCANSNPFILTEROUT + '{sample}.SNPpass.vcf.benchmark'
    shell:
        ('{config[tools][varscanNonSomaticFilter][call]} {input.vcf} {output.vcf} ' +
        '{params.minVarSupport} ' +
        '{params.pvalue} ' +
        '{params.minNucCoverage} ' +
        '{params.filterStrands} ' +
        '{params.tumorFreqThreshold} ' +
        '{params.filterHomopolymer} ' +
        '{params.filterSilent} ' +
        '{params.filterCommon}')

		
if not 'PARSEANNOTATEDVCFIN' in globals():
    PARSEANNOTATEDVCFIN = GATKVARIANTCOMBINEOUT
if not 'PARSEANNOTATEDVCF_OUT' in globals():
    PARSEANNOTATEDVCF_OUT = OUTDIR + 'databaseQuery/'
# This rule calls a python script parsing an annotated vcf file in order to create output files necessary for postprocessing, including database queries
rule parseAnnotatedVCF_forPostprocessingSteps:
    input:
        vcf = PARSEANNOTATEDVCFIN + '{tumor}_vs_{normal}.{placeholder}.vcf',
	pathwayDB = {config['resources'][ORGANISM]['pathwayDB']}
    output:
        outList_all = PARSEANNOTATEDVCF_OUT + '{tumor}_vs_{normal}.{placeholder}.overview.txt',
	outList_genes = PARSEANNOTATEDVCF_OUT + '{tumor}_vs_{normal}.{placeholder}.overview.txt_distinctGeneNames.txt',
	outList_cBioportal = PARSEANNOTATEDVCF_OUT + '{tumor}_vs_{normal}.{placeholder}.overview.txt_listForCbioportalQuery.txt'
    params:
        tumor = '{tumor}',
        lsfoutfile = PARSEANNOTATEDVCF_OUT + '{tumor}_vs_{normal}.{placeholder}.overview.lsfout.log',
        lsferrfile = PARSEANNOTATEDVCF_OUT + '{tumor}_vs_{normal}.{placeholder}.overview.lsferr.log',
        scratch = config['tools']['parseAnnotatedGenesFromVCF']['scratch'],
        mem = config['tools']['parseAnnotatedGenesFromVCF']['mem'],
        time = config['tools']['parseAnnotatedGenesFromVCF']['time']
    threads:
        config['tools']['parseAnnotatedGenesFromVCF']['threads']
    benchmark:
        PARSEANNOTATEDVCF_OUT + '{tumor}_vs_{normal}.{placeholder}.overview.benchmark'
    shell:
        'tumorName=$(basename {params.tumor}) ; {config[tools][parseAnnotatedGenesFromVCF][call]} {input.vcf} {output.outList_all} {input.pathwayDB} $tumorName'

if not 'SOMATICSIGNATURES_IN' in globals():
    SOMATICSIGNATURES_IN = GATKVARIANTCOMBINEOUT
if not 'SOMATICSIGNATURES_OUT' in globals():
    SOMATICSIGNATURES_OUT = OUTDIR + 'variants/somaticSignatures/'

# Compute somatic signatures and compare to cosmic signatures30
rule somaticSignatures:
    input:
        vcf = SOMATICSIGNATURES_IN + '{tumor}_vs_{normal}.{placeholder}.vcf'
    output:
        outFile= SOMATICSIGNATURES_OUT + '{tumor}_vs_{normal}.{placeholder}.somaticSignatures.txt'
    params:
        genomeVersion = config['tools']['somaticSignatures']['genomeVersion'],
        signaturesDirectory = config['tools']['somaticSignatures']['signaturesDirectory'],
        lsfoutfile = SOMATICSIGNATURES_OUT + '{tumor}_vs_{normal}.{placeholder}.somaticSignatures.lsfout.log',
        lsferrfile = SOMATICSIGNATURES_OUT + '{tumor}_vs_{normal}.{placeholder}.somaticSignatures.lsferr.log',
        scratch = config['tools']['somaticSignatures']['scratch'],
        mem = config['tools']['somaticSignatures']['mem'],
        time = config['tools']['somaticSignatures']['time'],
        tumorName = '{tumor}'
    threads:
        config['tools']['somaticSignatures']['threads']
    benchmark:
        SOMATICSIGNATURES_OUT + '{tumor}_vs_{normal}.{placeholder}.somaticSignatures.benchmark'
    shell:
        '{config[tools][somaticSignatures][call]} {input.vcf} {output.outFile} {params.genomeVersion} {params.signaturesDirectory} {params.tumorName}'


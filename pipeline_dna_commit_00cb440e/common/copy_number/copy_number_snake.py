### Rules for Excavator2

if not 'EXCAVATOR2IN' in globals():
    EXCAVATOR2IN = BASERECALIBRATIONOUT
if not 'EXCAVATOR2OUT' in globals():
    EXCAVATOR2OUT = OUTDIR + 'copynumber/excavator2/'

rule createExcavator2InputFiles_readInput:
    output:
        outReadInput = EXCAVATOR2OUT + '{tumor}_vs_{normal}_ReadInputExcavator2.txt'
    params:
        lsfoutfile = EXCAVATOR2OUT + '{tumor}_vs_{normal}.createExcavator2ReadInputFile.lsfout.log',
        lsferrfile = EXCAVATOR2OUT + '{tumor}_vs_{normal}.createExcavator2ReadInputFile.lsferr.log',
        scratch = config['tools']['excavator2']['createExcavator2ReadInputFile']['scratch'],
        mem = config['tools']['excavator2']['createExcavator2ReadInputFile']['mem'],
        time = config['tools']['excavator2']['createExcavator2ReadInputFile']['time'],
        out = EXCAVATOR2OUT,
        tumorBam = EXCAVATOR2IN + '{tumor}.bam',
        normalBam = EXCAVATOR2IN + '{normal}.bam'
    threads:
        config['tools']['excavator2']['createExcavator2ReadInputFile']['threads']
    benchmark:
        EXCAVATOR2OUT + '{tumor}_vs_{normal}.createExcavator2ReadInputFile.benchmark'
    log:
        EXCAVATOR2OUT + '{tumor}_vs_{normal}.createExcavator2ReadInputFile.log'
    shell:
        ('{config[tools][excavator2][createExcavator2ReadInputFile][call]} ' +
        '{output.outReadInput} ' +
        '{params.out} ' +
        '{params.tumorBam} ' +
        '{params.normalBam}')


# create the necessary input files
rule createExcavator2InputFiles_targetFile:
    input:
        wiggle = config['resources'][ORGANISM]['excavator2Wiggle'],
        reference = config['resources'][ORGANISM]['reference']
    output:
        outTarget =  EXCAVATOR2OUT + config['tools']['excavator2']['targetName'] + '_Excavator2.txt'
    params:
        lsfoutfile = EXCAVATOR2OUT + 'createExcavator2TargetFile.lsfout.log',
        lsferrfile = EXCAVATOR2OUT + 'createExcavator2TargetFile.lsferr.log',
        scratch = config['tools']['excavator2']['createExcavator2TargetFile']['scratch'],
        mem = config['tools']['excavator2']['createExcavator2TargetFile']['mem'],
        time = config['tools']['excavator2']['createExcavator2TargetFile']['time']
    threads:
        config['tools']['excavator2']['createExcavator2TargetFile']['threads']
    benchmark:
        EXCAVATOR2OUT + 'createExcavator2TargetFile.benchmark'
    log:
        EXCAVATOR2OUT + 'createExcavator2TargetFile.log'
    shell:
        ('{config[tools][excavator2][createExcavator2TargetFile][call]} ' +
        '{output.outTarget} ' +
        '{input.wiggle} ' +
        '{input.reference}')


# call excavator2 for exome sequencing copy number analysis
# first Target creation
rule excavator2TargetCreation:
    input:
        targetTXT = EXCAVATOR2OUT + config['tools']['excavator2']['targetName'] + '_Excavator2.txt',
        bedfile = config['resources'][ORGANISM]['regions']
    output:
        out = EXCAVATOR2OUT + config['tools']['excavator2']['targetName'] + '.targetSuccess.txt'
    params:
        lsfoutfile = EXCAVATOR2OUT + 'targetCreation.lsfout.log',
        lsferrfile = EXCAVATOR2OUT + 'targetCreation.lsferr.log',
        scratch = config['tools']['excavator2']['targetCreation']['scratch'],
        mem = config['tools']['excavator2']['targetCreation']['mem'],
        time = config['tools']['excavator2']['targetCreation']['time'],
        assembly = config['tools']['excavator2']['assembly'],
        window = config['tools']['excavator2']['window'],
        targetName = config['tools']['excavator2']['targetName']
    threads:
        config['tools']['excavator2']['targetCreation']['threads']
    benchmark:
        EXCAVATOR2OUT + 'targetCreation.benchmark'
    log:
        EXCAVATOR2OUT + 'targetCreation.log'
    shell:
        ('{config[tools][excavator2][targetCreation][call]} ' +
        '{input.targetTXT} ' +
        '{input.bedfile} ' +
        '{params.targetName} ' +
        '{params.window} '
        '{params.assembly} && ' +
        'touch {output.out}')


# call excavator2 for exome sequencing copy number analysis
# second copy call
rule excavator2CopyNum:
    input:
        tumor = EXCAVATOR2IN + '{tumor}.bam',
        normal = EXCAVATOR2IN + '{normal}.bam',
        targetCheck = EXCAVATOR2OUT + config['tools']['excavator2']['targetName'] + '.targetSuccess.txt',
        readInputTXT = EXCAVATOR2OUT + '{tumor}_vs_{normal}_ReadInputExcavator2.txt'
    output:
        out = EXCAVATOR2OUT + '{tumor}_vs_{normal}/' + '{tumor}_vs_{normal}.copyNumSuccess.txt',
        outResult = EXCAVATOR2OUT + '{tumor}_vs_{normal}/Results/{tumor}/FastCallResults_{tumor}.txt'
    params:
        lsfoutfile = EXCAVATOR2OUT + '{tumor}_vs_{normal}.copyNum.lsfout.log',
        lsferrfile = EXCAVATOR2OUT + '{tumor}_vs_{normal}.copyNum.lsferr.log',
        scratch = config['tools']['excavator2']['readInput']['scratch'],
        mem = config['tools']['excavator2']['readInput']['mem'],
        time = config['tools']['excavator2']['readInput']['time'],
        mode = config['tools']['excavator2']['readInput']['mode'],
        outputFolder = EXCAVATOR2OUT + '{tumor}_vs_{normal}',
        targetTXT = config['tools']['excavator2']['targetName'],
        assembly = config['tools']['excavator2']['assembly']
    threads:
        config['tools']['excavator2']['readInput']['threads']
    benchmark:
        EXCAVATOR2OUT + '{tumor}_vs_{normal}/' + '{tumor}_vs_{normal}.copyNum.benchmark'
    log:
	    EXCAVATOR2OUT + '{tumor}_vs_{normal}/' + '{tumor}_vs_{normal}.copyNum.log'
    shell:
        ('{config[tools][excavator2][readInput][call]} ' +
        '{input.readInputTXT} ' +
        '--processors {threads} ' +
        '--target {params.targetTXT} ' +
        '--assembly {params.assembly}')

# This rule reformats the excavator2 FastCall result into bed format
rule reformatExcavator2Result:
    input:
        inRes = EXCAVATOR2OUT + '{tumor}_vs_{normal}/Results/{tumor}/FastCallResults_{tumor}.txt'
    output:
        out = EXCAVATOR2OUT + '{tumor}_vs_{normal}.excavator2CNV.txt'
    params:
        lsfoutfile = EXCAVATOR2OUT + '{tumor}_vs_{normal}.excavator2CNV.lsfout.log',
        lsferrfile = EXCAVATOR2OUT + '{tumor}_vs_{normal}.excavator2CNV.lsferr.log',
        scratch = config['tools']['excavator2']['reformat']['scratch'],
        mem = config['tools']['excavator2']['reformat']['mem'],
        time = config['tools']['excavator2']['reformat']['time']
    threads:
        config['tools']['excavator2']['reformat']['threads']
    benchmark:
        EXCAVATOR2OUT + '{tumor}_vs_{normal}.excavator2CNV.benchmark'
    shell:
        'cp {input.inRes} {output.out} ; sed -i "s/Chromosome/#Chromosome/g" {output.out}'

### Rules for Excavator
# call excavator for exome sequencing copy number analysis
# create the necessary input files
rule createExcavatorInputFiles_readInput:
    output:
        outReadInput = EXCAVATOROUT + '{tumor}_vs_{normal}_ReadInputExcavator.txt'
    params:
        lsfoutfile = EXCAVATOROUT + '{tumor}_vs_{normal}.createExcavatorReadInputFile.lsfout.log',
        lsferrfile = EXCAVATOROUT + '{tumor}_vs_{normal}.createExcavatorReadInputFile.lsferr.log',
        scratch = config['tools']['excavator']['createExcavatorReadInputFile']['scratch'],
        mem = config['tools']['excavator']['createExcavatorReadInputFile']['mem'],
        time = config['tools']['excavator']['createExcavatorReadInputFile']['time'],
        assembly = config['tools']['excavator']['assembly'],
        experimentTarget = config['tools']['excavator']['targetName'],
        tumorBam = EXCAVATORIN + '{tumor}.bam',
        normalBam = EXCAVATORIN + '{normal}.bam'
    threads:
        config['tools']['excavator']['createExcavatorReadInputFile']['threads']
    benchmark:
        EXCAVATOROUT + '{tumor}_vs_{normal}.createExcavatorReadInputFile.benchmark'
    log:
        EXCAVATOROUT + '{tumor}_vs_{normal}.createExcavatorReadInputFile.log'
    shell:
        ('{config[tools][excavator][createExcavatorReadInputFile][call]} ' +
        '{output.outReadInput} ' +
        '{params.experimentTarget} ' +
        '{params.assembly} ' +
        '{params.tumorBam} ' +
        '{params.normalBam}')
    
# create the necessary input files
rule createExcavatorInputFiles_targetFile:
    input:
        wiggle = config['resources'][ORGANISM]['excavatorWiggle'],
        reference = config['resources'][ORGANISM]['reference']
    output:
        outTarget =  EXCAVATOROUT + config['tools']['excavator']['targetName'] + '_Excavator.txt'
    params:
        lsfoutfile = EXCAVATOROUT + 'createExcavatorTargetFile.lsfout.log',
        lsferrfile = EXCAVATOROUT + 'createExcavatorTargetFile.lsferr.log',
        scratch = config['tools']['excavator']['createExcavatorTargetFile']['scratch'],
        mem = config['tools']['excavator']['createExcavatorTargetFile']['mem'],
        time = config['tools']['excavator']['createExcavatorTargetFile']['time']
    threads:
        config['tools']['excavator']['createExcavatorTargetFile']['threads']
    benchmark:
        EXCAVATOROUT + 'createExcavatorTargetFile.benchmark'
    log:
        EXCAVATOROUT + 'createExcavatorTargetFile.log'
    shell:
        ('{config[tools][excavator][createExcavatorTargetFile][call]} ' +
        '{output.outTarget} ' +
        '{input.wiggle} ' +
        '{input.reference}')

# call excavator for exome sequencing copy number analysis
# first Target creation
rule excavatorTargetCreation:
    input:
        targetTXT = EXCAVATOROUT + config['tools']['excavator']['targetName'] + '_Excavator.txt',
        bedfile = config['resources'][ORGANISM]['regions']
    output:
        out = EXCAVATOROUT + config['tools']['excavator']['targetName'] + '.targetSuccess.txt'
    params:
        lsfoutfile = EXCAVATOROUT + 'targetCreation.lsfout.log',
        lsferrfile = EXCAVATOROUT + 'targetCreation.lsferr.log',
        scratch = config['tools']['excavator']['targetCreation']['scratch'],
        mem = config['tools']['excavator']['targetCreation']['mem'],
        time = config['tools']['excavator']['targetCreation']['time'],
        assembly = config['tools']['excavator']['assembly'],
	targetName = config['tools']['excavator']['targetName']
    threads:
        config['tools']['excavator']['targetCreation']['threads']
    benchmark:
        EXCAVATOROUT + 'targetCreation.benchmark'
    log:
        EXCAVATOROUT + 'targetCreation.log'
    shell:
        ('{config[tools][excavator][targetCreation][call]} ' +
        '{input.targetTXT} ' +
        '{input.bedfile} ' +
        '{params.targetName} ' +
        '--assembly {params.assembly} && ' +
        'touch {output.out}')
        
# call excavator for exome sequencing copy number analysis
# second copy call
rule excavatorCopyNum:
    input:
        tumor = EXCAVATORIN + '{tumor}.bam',
        normal = EXCAVATORIN + '{normal}.bam',
        targetCheck = EXCAVATOROUT + config['tools']['excavator']['targetName'] + '.targetSuccess.txt',
        readInputTXT = EXCAVATOROUT + '{tumor}_vs_{normal}_ReadInputExcavator.txt'
    output:
        out = EXCAVATOROUT + '{tumor}_vs_{normal}/' + '{tumor}_vs_{normal}.copyNumSuccess.txt',
        outResult = EXCAVATOROUT + '{tumor}_vs_{normal}/Results/{tumor}/FastCallResults_{tumor}.txt' 
    params:
        lsfoutfile = EXCAVATOROUT + '{tumor}_vs_{normal}.copyNum.lsfout.log',
        lsferrfile = EXCAVATOROUT + '{tumor}_vs_{normal}.copyNum.lsferr.log',
        scratch = config['tools']['excavator']['readInput']['scratch'],
        mem = config['tools']['excavator']['readInput']['mem'],
        time = config['tools']['excavator']['readInput']['time'],
        mode = config['tools']['excavator']['readInput']['mode'],
        outputFolder = EXCAVATOROUT + '{tumor}_vs_{normal}'
    threads:
        config['tools']['excavator']['readInput']['threads']
    benchmark:
        EXCAVATOROUT + '{tumor}_vs_{normal}/' + '{tumor}_vs_{normal}.copyNum.benchmark'
    log:
        EXCAVATOROUT + '{tumor}_vs_{normal}/' + '{tumor}_vs_{normal}.copyNum.log'
    shell:
        ('{config[tools][excavator][readInput][call]} ' +
        '{input.readInputTXT} ' +
        '{params.outputFolder} ' +
        '--mode {params.mode} && ' +
        'touch {output.out}')
		

# This rule reformats the excavator FastCall result into bed format
rule reformatExcavatorResult:
    input:
        inRes = EXCAVATOROUT + '{tumor}_vs_{normal}/Results/{tumor}/FastCallResults_{tumor}.txt'
    output:
        out = EXCAVATOROUT + '{tumor}_vs_{normal}.excavatorCNV.txt'
    params:
        lsfoutfile = EXCAVATOROUT + '{tumor}_vs_{normal}.excavatorCNV.lsfout.log',
        lsferrfile = EXCAVATOROUT + '{tumor}_vs_{normal}.excavatorCNV.lsferr.log',
        scratch = config['tools']['excavator']['reformat']['scratch'],
        mem = config['tools']['excavator']['reformat']['mem'],
        time = config['tools']['excavator']['reformat']['time']
    threads:
        config['tools']['excavator']['reformat']['threads']
    benchmark:
        EXCAVATOROUT + '{tumor}_vs_{normal}.excavatorCNV.benchmark'
    shell:
        'cp {input.inRes} {output.out} ; sed -i "s/Chromosome/#Chromosome/g" {output.out}'

if not 'FACETSIN' in globals():
    FACETSIN = BASERECALIBRATIONOUT
if not 'FACETSOUT' in globals():
    FACETSOUT = OUTDIR + 'copy_number/facets/'

rule facets_filter:
    input:
        cn = FACETSOUT + '{tumor}_vs_{normal}.cn'
    output:
        filteredCN = FACETSOUT + '{tumor}_vs_{normal}_filtered.cn'
    params:
        lsfoutfile = FACETSOUT + '{tumor}_vs_{normal}_filtered.cn.lsfout.log',
        lsferrfile = FACETSOUT + '{tumor}_vs_{normal}_filtered.cn.lsferr.log',
        scratch = config['tools']['facets']['facets_filter']['scratch'],
        mem = config['tools']['facets']['facets_filter']['mem'],
        time = config['tools']['facets']['facets_filter']['time'],
        colName_totalCopy = config['tools']['facets']['facets_filter']['colName_totalCopy']
    threads:
        config['tools']['facets']['facets_filter']['threads']
    benchmark:
        FACETSOUT + '{tumor}_vs_{normal}_filtered.cn.benchmark'
    shell:
        ('{config[tools][facets][facets_filter][call]} ' + 
        '--infile {input.cn} --outfile {output.filteredCN} ' +
        '--colName_totalCopy {params.colName_totalCopy}')

# calls bedtools intersect to annotate genes to the copy number regions
rule annotate_facets_CNVs_withBedtools:
    input:
        inRes = '{sample}.cn',
        inDB = config['resources'][ORGANISM]['geneAnnotationDB'] 
    output:
        out = '{sample}_annotated.txt'
    params:
        lsfoutfile = '{sample}.annotated.lsfout.log',
        lsferrfile = '{sample}.annotated.lsferr.log',
        scratch = config['tools']['bedtools']['intersect']['scratch'],
        mem = config['tools']['bedtools']['intersect']['mem'],
        time = config['tools']['bedtools']['intersect']['time']
    threads:
        config['tools']['bedtools']['intersect']['threads']
    benchmark:
        '{sample}.annotated.benchmark'
    shell:
        '{config[tools][bedtools][intersect][call]} intersect -a {input.inRes} -b {input.inDB} -wa -wb > {output.out}.temp1.txt ; ' +
        '{config[tools][bedtools][intersect][call]} intersect -a {input.inRes} -b {input.inDB} -v -wa > {output.out}.temp2.txt ; ' +
        'cat {output.out}.temp1.txt {output.out}.temp2.txt | sort -k1,1 -k2,2n > {output.out}'


if not 'DATABASEQUERY' in globals():
    DATABASEQUERY = OUTDIR + 'databaseQuery/'


# This rule calls a python script parsing and filtering an annotated facets cnv file in order to create output files necessary for postprocessing, 
# including database queries; script necessary for Facets WES copy number calls
rule parseAndFilter_Facets_Annotated:
    input:
        infile = FACETSOUT + '{sample}_annotated.txt'
    output:
        outList_all = DATABASEQUERY + '{sample}.overviewFiltered.txt',
	outList_genes = DATABASEQUERY + '{sample}.overviewFiltered.txt_distinctGeneNames.txt',
	outList_cBioportal = DATABASEQUERY + '{sample}.overviewFiltered.txt_listForCbioportalQuery.txt'
    params:
        lsfoutfile = DATABASEQUERY + '{sample}.parseAndFilter_Facets_Annotated.lsfout.log',
        lsferrfile = DATABASEQUERY + '{sample}.parseAndFilter_Facets_Annotated.lsferr.log',
        scratch = config['tools']['parseAndFilter_Facets_Annotated']['scratch'],
        mem = config['tools']['parseAndFilter_Facets_Annotated']['mem'],
        time = config['tools']['parseAndFilter_Facets_Annotated']['time'],
        probThreshold = config['tools']['parseAndFilter_Facets_Annotated']['probabilityThreshold']
    threads:
        config['tools']['parseAndFilter_Facets_Annotated']['threads']
    benchmark:
        DATABASEQUERY + '{sample}.parseAndFilter_Facets_Annotated.benchmark'
    shell:
        '{config[tools][parseAndFilter_Facets_Annotated][call]} {input.infile} {output.outList_all} {params.probThreshold}'


# This rule calls a python script parsing and filtering an annotated vcf file in order to create output files necessary for postprocessing, 
# including database queries; script necessary for Excavator WES copy number calls
rule parseAndFilterExcavatorAnnotated:
    input:
        infile = EXCAVATOROUT + '{sample}.txt'
    output:
        outList_all = DATABASEQUERY + '{sample}.overviewFiltered.txt',
	outList_genes = DATABASEQUERY + '{sample}.overviewFiltered.txt_distinctGeneNames.txt',
	outList_cBioportal = DATABASEQUERY + '{sample}.overviewFiltered.txt_listForCbioportalQuery.txt'
    params:
        lsfoutfile = DATABASEQUERY + '{sample}.parseAndFilter.lsfout.log',
        lsferrfile = DATABASEQUERY + '{sample}.parseAndFilter.lsferr.log',
        scratch = config['tools']['parseAndFilterExcavatorAnnotated']['scratch'],
        mem = config['tools']['parseAndFilterExcavatorAnnotated']['mem'],
        time = config['tools']['parseAndFilterExcavatorAnnotated']['time'],
        probThreshold = config['tools']['parseAndFilterExcavatorAnnotated']['probabilityThreshold']
    threads:
        config['tools']['parseAndFilterExcavatorAnnotated']['threads']
    benchmark:
        DATABASEQUERY + '{sample}.parseAndFilter.benchmark'
    shell:
        '{config[tools][parseAndFilterExcavatorAnnotated][call]} {input.infile} {output.outList_all} {params.probThreshold}'
		
if not 'BICSEQ2OUT' in globals():
    BICSEQ2OUT = OUTDIR + 'copynumber/bicseq2/'

# This rule calls a python script parsing an annotated vcf file in order to create output files necessary for postprocessing, including database queries
# Script necessary for BICSeq2 WGS copy number calls
rule processAnnotatedBicSeq2:
    input:
        infile = BICSEQ2OUT + '{sample}.txt'
    output:
        outList_all = DATABASEQUERY + '{sample}.overview.txt',
	outList_genes = DATABASEQUERY + '{sample}.overview.txt_distinctGeneNames.txt',
	outList_cBioportal = DATABASEQUERY + '{sample}.overview.txt_listForCbioportalQuery.txt'
    params:
        lsfoutfile = DATABASEQUERY + '{sample}.parseForOverview.lsfout.log',
        lsferrfile = DATABASEQUERY + '{sample}.parseForOverview.lsferr.log',
        scratch = config['tools']['processAnnotatedBicSeq2']['scratch'],
        mem = config['tools']['processAnnotatedBicSeq2']['mem'],
        time = config['tools']['processAnnotatedBicSeq2']['time']
    threads:
        int(config['tools']['processAnnotatedBicSeq2']['threads'])
    benchmark:
        DATABASEQUERY + '{sample}.parseForOverview.benchmark'
    shell:
        '{config[tools][processAnnotatedBicSeq2][call]} {input.infile} {output.outList_all}'


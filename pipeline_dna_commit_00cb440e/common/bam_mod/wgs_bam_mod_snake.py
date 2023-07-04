import ntpath

# BamClipOverlap - NGSbits
# Softclipping overlapping reads.
"""
rule bamClipOverlap:
	input:
		bam = CLIPOVERLAPIN + '{sample}.bam'
	output:
		bam_clipped = CLIPOVERLAPOUT + '{sample}.bam'
	params:
		lsfoutfile = CLIPOVERLAPOUT + '{sample}.bam.lsfout.log',
        	lsferrfile = CLIPOVERLAPOUT + '{sample}.bam.lsferr.log',
        	scratch = config['tools']['clipoverlap']['scratch'],
        	mem = config['tools']['clipoverlap']['mem'],
        	time = config['tools']['clipoverlap']['time'],
	threads:
		int(config['tools']['clipoverlap']['threads'])
	benchmark:
		CLIPOVERLAPOUT + '{sample}.clipoverlap.benchmark'
	shell:
		'{config[tools][clipoverlap][call]} -in {input.bam} -out {output.bam_clipped}'


# This rule sorts a bam file
rule NameSort:
	input:
        	bam=NAMESORTEDIN + '{sample}.bam',
		bai=NAMESORTEDIN + '{sample}.bam.bai'
	output:
        	bam=temp(NAMESORTEDOUT + '{sample}.bam')
	params:
		out = NAMESORTEDOUT,
		lsfoutfile = NAMESORTEDOUT + '{sample}.lsfout.log',
        	lsferrfile = NAMESORTEDOUT + '{sample}.lsferr.log',
        	scratch = config['tools']['picard']['SortSam']['scratch'],
        	mem = config['tools']['picard']['SortSam']['mem'],
        	time = config['tools']['picard']['SortSam']['time'],
        	max_records_in_ram = config['tools']['picard']['SortSam']['max_records_in_ram']
	benchmark:
		NAMESORTEDOUT + '{sample}.benchmark'
	threads:
        	int(config['tools']['picard']['SortSam']['threads'])
	shell:
        	('{config[tools][picard][call]} SortSam ' +
        	'INPUT={input.bam} ' +
        	'OUTPUT={output.bam} ' +
        	'SORT_ORDER= queryname ' +
        	'MAX_RECORDS_IN_RAM={params.max_records_in_ram} ' +
        	'TMP_DIR={TMPDIR} ; cp {input.bai} {params.out}' )

rule Sort:
    input:
        bam=SORTEDIN + '{sample}.bam'
    output:
        bam=temp(SORTEDOUT + '{sample}.bam')
    params:
        lsfoutfile = SORTEDOUT + '{sample}.lsfout.log',
        lsferrfile = SORTEDOUT + '{sample}.lsferr.log',
        scratch = config['tools']['picard']['SortSam']['scratch'],
        mem = config['tools']['picard']['SortSam']['mem'],
        time = config['tools']['picard']['SortSam']['time'],
        max_records_in_ram = config['tools']['picard']['SortSam']['max_records_in_ram']
    benchmark:
        NAMESORTEDOUT + '{sample}.benchmark'
    threads:
        int(config['tools']['picard']['SortSam']['threads'])
    shell:
        ('{config[tools][picard][call]} SortSam ' +
        'INPUT={input.bam} ' +
        'OUTPUT={output.bam} ' +
        'SORT_ORDER= coordinate ' +
	'VALIDATION_STRINGENCY=LENIENT ' +
        'MAX_RECORDS_IN_RAM={params.max_records_in_ram} ' +
        'TMP_DIR={TMPDIR} ')

"""

# This rule creates an indes of a BAM file
rule createIndex:
    input:
        bam = '{sample}.bam',
    output:
        idx = '{sample}.bai',
    params:
        lsfoutfile = '{sample}.bai.lsfout.log',
        lsferrfile = '{sample}.bai.lsferr.log',
        scratch = config['tools']['samtools']['index']['scratch'],
        mem = config['tools']['samtools']['index']['mem'],
        time = config['tools']['samtools']['index']['time']
    threads:
        int(config['tools']['samtools']['index']['threads'])
    benchmark:
        '{sample}.bai.benchmark'
    shell:
        '{config[tools][samtools][call]} index {input.bam} && mv {input.bam}.bai {output.idx}'

rule linkIndex:
    input:
        bai = '{sample}.bai',
    output:
        bai = '{sample}.bam.bai',
    params:
        lsfoutfile = '{sample}.bam.bai.lsfout.log',
        lsferrfile = '{sample}.bam.bai.lsferr.log',
        scratch = '1000', # config['tools']['samtools']['index']['scratch'],
        mem = '1000', # config['tools']['samtools']['index']['mem'],
        time = '1' #config['tools']['samtools']['index']['time']
    threads:
        1
    benchmark:
        '{sample}.bam.bai.benchmark'
    shell:
        'ln -s {input.bai} {output.bai}'

# This rule sorts a BAM file and fixes mate pair information if necessary
# This is a picard tool
rule fixMatePairAndSort:
    input:
        bam=FIXMATEANDSORTIN + '{sample}.bam'
    output:
        bam=temp(FIXMATEANDSORTOUT + '{sample}.bam')
    params:
        lsfoutfile = FIXMATEANDSORTOUT + '{sample}.lsfout.log',
        lsferrfile = FIXMATEANDSORTOUT + '{sample}.lsferr.log',
        scratch = config['tools']['picard']['fixMateInformation']['scratch'],
        mem = config['tools']['picard']['fixMateInformation']['mem'],
        time = config['tools']['picard']['fixMateInformation']['time'],
        sortOrder = config['tools']['picard']['fixMateInformation']['sortOrder'],
        assume_sorted = config['tools']['picard']['fixMateInformation']['assume_sorted'],
        max_records_in_ram = config['tools']['picard']['fixMateInformation']['max_records_in_ram']
    benchmark:
        FIXMATEANDSORTOUT + '{sample}.benchmark'
    threads:
        int(config['tools']['picard']['fixMateInformation']['threads'])
    log:
        FIXMATEANDSORTOUT + '{sample}.log'
    shell:
        ('{config[tools][picard][call]} FixMateInformation ' +
        'INPUT={input.bam} ' +
        'OUTPUT={output.bam} ' +
        'SORT_ORDER={params.sortOrder} ' +
        'ASSUME_SORTED={params.assume_sorted} ' +
        'MAX_RECORDS_IN_RAM={params.max_records_in_ram} ' +
        'TMP_DIR={TMPDIR} ' +
        '2> {log}')

# This functiom creates a list of BAM files created by the read mapper
def getAlignerBams():
    out = []
    if CLIPPEDOUT == FASTQDIR or WHICHCLIP is True:
        for f in PAIREDFASTQFILESWITHOUTR:
            out.append( f + '.bam')
    else:
        for f in PAIREDFASTQFILES:
            out.append(os.path.dirname(f) +'/ORPHAN/' + ntpath.basename(f) + '.bam')
        for f in PAIREDFASTQFILESWITHOUTR:
            out.append( f + '.bam')
    return out

# This function is a helper function to get the BAM file names which are then merged.
def getBamsToMerge(wildcards):
    out = []
    allBams = getAlignerBams()
    for bam in allBams:
        #print(bam.split("/")[0])
        if wildcards.sample == bam.split("/")[0]: #[0:len(wildcards.sample)]:
            out.append(MERGEBAMSIN + bam)
    if not out:
        return ['ERROR']
    return out

# This function prepends 'INPUT=' in front of every BAM file that is to be merged.
def prependBamsToMerge(wildcards):
    bamsToMerge = getBamsToMerge(wildcards)
    return ''.join(['INPUT='+bam+' ' for bam in bamsToMerge])

# This rule removes all alignments which are marked as secondary alignment
rule mergeBams:
    input:
        bams = getBamsToMerge
    output:
        bam=MERGEBAMSOUT + '{sample}.bam'
    params:
        lsfoutfile = MERGEBAMSOUT + '{sample}.lsfout.log',
        lsferrfile = MERGEBAMSOUT + '{sample}.lsferr.log',
        assume_sorted = config['tools']['picard']['mergeBams']['assume_sorted'],
        scratch = config['tools']['picard']['mergeBams']['scratch'],
        mem = config['tools']['picard']['mergeBams']['mem'],
        time = config['tools']['picard']['mergeBams']['time'],
        input = prependBamsToMerge
    threads:
        int(config['tools']['picard']['mergeBams']['threads'])
    benchmark:
        MERGEBAMSOUT + '{sample}.benchmark'
    log:
        log = MERGEBAMSOUT + '{sample}.log',
        metrics = MERGEBAMSOUT + '{sample}.metrics'
    shell:
        ('{config[tools][picard][call]} ' +
        'MergeSamFiles ' +
        '{params.input} ' +
        'OUTPUT={output.bam} ' +
        'ASSUME_SORTED={params.assume_sorted} ' +
        '2> {log.log}')

# This rule removes all alignments which are marked as secondary alignment
rule removeSecondaryAlignments:
    input:
        bam = NOSECONDARYALNIN+ '{sample}.bam',
    output:
        bam=temp(NOSECONDARYALNOUT + '{sample}.bam'),
    params:
        lsfoutfile = NOSECONDARYALNOUT + '{sample}.lsfout.log',
        lsferrfile = NOSECONDARYALNOUT + '{sample}.lsferr.log',
        scratch = config['tools']['samtools']['rmSecondary']['scratch'],
        mem = config['tools']['samtools']['rmSecondary']['mem'],
        time = config['tools']['samtools']['rmSecondary']['time']
    benchmark:
        NOSECONDARYALNOUT + '{sample}.benchmark'
    threads:
        int(config['tools']['samtools']['rmSecondary']['threads'])
    shell:
        '{config[tools][samtools][call]} view -bh -F 256 {input.bam} > {output.bam}'

# This rule markes PCR duplicates
# This is a picard tool
if not 'MARKPCRDUBLICATESIN' in globals():
    MARKPCRDUBLICATESIN = NOSECONDARYALNOUT
if not 'MARKPCRDUBLICATESOUT' in globals():
    MARKPCRDUBLICATESOUT = OUTDIR + 'markedDuplicates/'

rule markPCRDuplicates:
    input:
        bam=MARKPCRDUBLICATESIN + '{sample}.bam',
    output:
        bam=temp(MARKPCRDUBLICATESOUT + '{sample}.bam'),
    params:
        lsfoutfile = MARKPCRDUBLICATESOUT + '{sample}.lsfout.log',
        lsferrfile = MARKPCRDUBLICATESOUT + '{sample}.lsferr.log',
        removeDuplicates = config['tools']['picard']['markduplicates']['removeDuplicates'],
        createIndex = config['tools']['picard']['markduplicates']['createIndex'],
        scratch = config['tools']['picard']['markduplicates']['scratch'],
        mem = config['tools']['picard']['markduplicates']['mem'],
        time = config['tools']['picard']['markduplicates']['time'],
        assume_sorted = config['tools']['picard']['markduplicates']['assume_sorted'],
        max_records_in_ram = config['tools']['picard']['markduplicates']['max_records_in_ram'],
        max_file_handles_for_read_ends_map = config['tools']['picard']['markduplicates']['max_file_handles_for_read_ends_map']
    threads:
        int(config['tools']['picard']['markduplicates']['threads'])
    benchmark:
        MARKPCRDUBLICATESOUT + '{sample}.benchmark'
    log:
        log = MARKPCRDUBLICATESOUT + '{sample}.log',
        metrics = MARKPCRDUBLICATESOUT + '{sample}.metrics'
    shell:
        ('{config[tools][picard][call]} ' +
        'MarkDuplicates ' +
        'INPUT={input.bam} ' +
        'OUTPUT={output.bam}  ' +
        'ASSUME_SORTED={params.assume_sorted} ' +
        'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP={params.max_file_handles_for_read_ends_map} ' +
        'MAX_RECORDS_IN_RAM={params.max_records_in_ram} ' +
        'TMP_DIR={TMPDIR} ' +
        'METRICS_FILE={log.metrics} ' +
        'REMOVE_DUPLICATES={params.removeDuplicates} ' +
        'CREATE_INDEX={params.createIndex} ' +
        '2> {log.log}')

# This rule removed the previously marked PCR duplicates
if not 'REMOVEPCRDUBLICATESIN' in globals():
    REMOVEPCRDUBLICATESIN = MARKPCRDUBLICATESOUT
if not 'REMOVEPCRDUBLICATESOUT' in globals():
    REMOVEPCRDUBLICATESOUT = OUTDIR + 'removedPcrDuplicates/'

rule removePCRDuplicates:
    input:
        bam=REMOVEPCRDUBLICATESIN + '{sample}.bam',
    output:
        bam=REMOVEPCRDUBLICATESOUT + '{sample}.bam',
    params:
        lsfoutfile = REMOVEPCRDUBLICATESOUT + '{sample}.lsfout.log',
        lsferrfile = REMOVEPCRDUBLICATESOUT + '{sample}.lsferr.log',
        scratch = config['tools']['samtools']['rmDuplicates']['scratch'],
        mem = config['tools']['samtools']['rmDuplicates']['mem'],
        time = config['tools']['samtools']['rmDuplicates']['time']
    threads:
        int(config['tools']['samtools']['rmDuplicates']['threads'])
    benchmark:
        REMOVEPCRDUBLICATESOUT + '{sample}.benchmark'
    shell:
        '{config[tools][samtools][call]} view -bh -F 0x400 {input.bam} > {output.bam}'

# Rule to replace a certain mapping quality with a specified one
# This is a GATK tool
if not 'REASSIGNONEMAPPINGQUALIN' in globals():
    REASSIGNONEMAPPINGQUALIN = OUTDIR + '.reassing_one_mapping_quality_in'
if not 'REASSIGNONEMAPPINGQUALOUT' in globals():
    REASSIGNONEMAPPINGQUALOUT = OUTDIR + '.reassing_one_mapping_quality_out'

rule reassignOneMappingQualityFilter:
    input:
        bam=REASSIGNONEMAPPINGQUALIN + '{sample}.bam',
        bamIdx=REASSIGNONEMAPPINGQUALIN + '{sample}.bam.index',
    output:
        bam=REASSIGNONEMAPPINGQUALOUT + '{sample}.bam',
    params:
        scratch = config['tools']['GATK']['baseRecalibrator']['scratch'],
        mem = config['tools']['GATK']['baseRecalibrator']['mem'],
        time = config['tools']['GATK']['baseRecalibrator']['time'],
        lsfoutfile = REASSIGNONEMAPPINGQUALOUT + '{sample}.lsfout.log',
        lsferrfile = REASSIGNONEMAPPINGQUALOUT + '{sample}.lsferr.log'
    benchmark:
        REASSIGNONEMAPPINGQUALOUT + '{sample}.benchmark'
    threads:
        int(config['tools']['GATK']['baseRecalibrator']['threads'])
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T PrintReads ' +
        '-R {config[tools][reference]} ' +
        '-I {input.bam} ' +
        '-o {output.bam} ' +
        '-rf ReassignOneMappingQuality ' +
        '-RMQF {config[tools][GATK][reassignOneMappingQualityFilter][oldQual]} ' +
        '-RMQT {config[tools][GATK][reassignOneMappingQualityFilter][newQual]}')

#rule realignTargetCreation:
#    input:
#        bam = REALIGNINDELSIN + '{sample}.bam',
#        bai = REALIGNINDELSIN + '{sample}.bai'
#    output:
#        bed = temp(REALIGNINDELSOUT + '{sample}.bed'),
#    params:
#        lsfoutfile = REALIGNINDELSOUT + '{sample}_target_creation.lsfout.log',
#        lsferrfile = REALIGNINDELSOUT + '{sample}_target_creation.lsferr.log',
#        reference = {config['resources'][ORGANISM]['reference']},
#        scratch = config['tools']['GATK']['realignTargetCreator']['scratch'],
#        mem = config['tools']['GATK']['realignTargetCreator']['mem'],
#        time = config['tools']['GATK']['realignTargetCreator']['time']
#    benchmark:
#        REALIGNINDELSOUT + '{sample}_target_creation.benchmark'
#    threads:
#        int(config['tools']['GATK']['realignTargetCreator']['threads'])
#    shell:
#        ('{config[tools][GATK][call]} ' +
#        '-T RealignerTargetCreator ' +
#        '-R {params.reference} ' +
#        '-I {input.bam} ' +
#        '{config[tools][GATK][realignTargetCreator][known1]} ' +
#        '{config[tools][GATK][realignTargetCreator][known2]} ' +
#        '-o {output.bed} ' +
#        '-nt {config[tools][GATK][realignTargetCreator][threads]}')
#
#
## Rule to perform the indel realignment
## This is a GATK tool
#ruleorder: realignIndels > createIndex
#rule realignIndels:
#    input:
#        bam = REALIGNINDELSIN + '{sample}.bam',
#        bai = REALIGNINDELSIN + '{sample}.bai',
#        interval = REALIGNINDELSOUT + '{sample}.bed',
#    output:
#        bam = REALIGNINDELSOUT + '{sample}.bam',
#        bai = REALIGNINDELSOUT + '{sample}.bai'
#    params:
#        lsfoutfile = REALIGNINDELSOUT + '{sample}.lsfout.log',
#        lsferrfile = REALIGNINDELSOUT + '{sample}.lsferr.log',
#        reference = {config['resources'][ORGANISM]['reference']},
#        scratch = config['tools']['GATK']['realignTargetCreator']['scratch'],
#        mem = config['tools']['GATK']['realignTargetCreator']['mem'],
#        time = config['tools']['GATK']['realignTargetCreator']['time']
#    benchmark:
#        REALIGNINDELSOUT + '{sample}.benchmark'
#    threads:
#        int(config['tools']['GATK']['realignTargetCreator']['threads'])
#    shell:
#        ('{config[tools][GATK][call]} ' +
#        '-T IndelRealigner ' +
#        '-R {params.reference} ' +
#        '-I {input.bam} ' +
#        '{config[tools][GATK][realignTargetCreator][known1]} ' +
#        '{config[tools][GATK][realignTargetCreator][known2]} ' +
#        '{config[tools][GATK][realignTargetCreator][known3]} ' +
#        '-o {output.bam} ' +
#        '-targetIntervals {input.interval}')
#
def getSamplesFromExperimentId(wildcards):
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    expMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            lineSplit = line.strip().split()
            exp = lineSplit[0]
            sample = lineSplit[1]
            sampleType = lineSplit[2]
            tpoint = lineSplit[3]
            if exp not in expMap.keys():
                expMap[exp] = []
            expMap[exp].append(sample)
    if wildcards.experiment not in expMap.keys():
        #print("Unknown experiment")
        return "UnknownExperiment"
    return expMap[wildcards.experiment]

def getBamsFromExperimentId(wildcards):
    return expand('{sample}.bam', sample = getSamplesFromExperimentId(wildcards))

def getBaisFromExperimentId(wildcards):
    return expand('{sample}.bai', sample = getSamplesFromExperimentId(wildcards))

def getBamsToRealingFromExperimentId(wildcards):
    return expand(REALIGNINDELSIN + '{bam}', bam = getBamsFromExperimentId(wildcards))
def getBaisToRealingFromExperimentId(wildcards):
    return expand(REALIGNINDELSIN + '{bai}', bai = getBaisFromExperimentId(wildcards))
def prependBamsToRealign(wildcards):
    bamsToRealign = getBamsToRealingFromExperimentId(wildcards)
    return ''.join(['-I '+bam+' ' for bam in bamsToRealign])

# Rule to create a file containing the regions to perfrom the indel realignment on
# This is a GATK tool
rule realignTargetCreation:
    input:
        bam = getBamsToRealingFromExperimentId,
        bai = getBaisToRealingFromExperimentId
    output:
        bed = temp(REALIGNINDELSOUT + '{experiment}.bed'),
    params:
        lsfoutfile = REALIGNINDELSOUT + '{experiment}_target_creation.lsfout.log',
        lsferrfile = REALIGNINDELSOUT + '{experiment}_target_creation.lsferr.log',
        reference = config['resources'][ORGANISM]['reference'],
        scratch = config['tools']['GATK']['realignTargetCreator']['scratch'],
        mem = config['tools']['GATK']['realignTargetCreator']['mem'],
        time = config['tools']['GATK']['realignTargetCreator']['time'],
        input = prependBamsToRealign,
        known1 = config['tools']['GATK']['realignTargetCreator']['known1'],
        known2 = config['tools']['GATK']['realignTargetCreator']['known2']
    benchmark:
        REALIGNINDELSOUT + '{experiment}_target_creation.benchmark'
    threads:
        int(config['tools']['GATK']['realignTargetCreator']['threads'])
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T RealignerTargetCreator ' +
        '-R {params.reference} ' +
        '{params.input} ' +
        '{params.known1} ' +
        '{params.known2} ' +
        '-o {output.bed} ' +
        '-nt {threads}')

localrules: createRealingIndelsInOutMapping
rule createRealingIndelsInOutMapping:
    output:
        exp = REALIGNINDELSOUT + '{experiment}.map'
    run:
        import sys
        try:
            outFile = open(output.exp, 'x')
            bams = getBamsFromExperimentId(wildcards)
            for bam in bams:
                outFile.write(bam + "\t" + REALIGNINDELSOUT + "ORIGINAL_" + bam + "\n")
        except IOError:
            print("Could not open file: ", output.exp)

# Rule to perform the indel realignment
# This is a GATK tool
#ruleorder: realignIndels > createIndex
rule realignIndels:
    input:
        bam = getBamsToRealingFromExperimentId,
        bai = getBaisToRealingFromExperimentId,
        map = REALIGNINDELSOUT + '{experiment}.map',
        interval = REALIGNINDELSOUT + '{experiment}.bed'
    output:
        txt = REALIGNINDELSOUT + '{experiment}.realigned.txt',
    params:
        lsfoutfile = REALIGNINDELSOUT + '{experiment}.lsfout.log',
        lsferrfile = REALIGNINDELSOUT + '{experiment}.lsferr.log',
        reference = config['resources'][ORGANISM]['reference'],
        scratch = config['tools']['GATK']['realignIndels']['scratch'],
        mem = config['tools']['GATK']['realignIndels']['mem'],
        time = config['tools']['GATK']['realignIndels']['time'],
        input = prependBamsToRealign,
        known1 = config['tools']['GATK']['realignIndels']['known1'],
        known2 = config['tools']['GATK']['realignIndels']['known2']
    benchmark:
        REALIGNINDELSOUT + '{experiment}.benchmark'
    threads:
        int(config['tools']['GATK']['realignIndels']['threads'])
    shell:
        ('{config[tools][GATK][call]} ' +
        '-T IndelRealigner ' +
        '-R {params.reference} ' +
        '{params.input} ' +
        '{params.known1} ' +
        '{params.known2} ' +
        '--nWayOut {input.map} ' +
        '-targetIntervals {input.interval} ' +
        '&& touch {output.txt}')

def getExperimentIdFromBam(wildcards):
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            lineSplit = line.strip().split()
            exp = lineSplit[0]
            sample = lineSplit[1]
            sampleType = lineSplit[2]
            tpoint = lineSplit[3]
            sampleMap[sample] = exp
    
    #print("Test1a: ", wildcards.sample)
    if wildcards.sample not in sampleMap.keys():
        #print("Unknown sample: ", wildcards.sample)
        return "UnknownSample"
    #print("Test1b: ",sampleMap[wildcards.sample])
    return sampleMap[wildcards.sample]

def getRealignedExperiment(wildcards):
    return expand(REALIGNINDELSOUT + '{experiment}.realigned.txt', experiment = getExperimentIdFromBam(wildcards))


# This is a "dummy" rule
localrules: getRealignedBam 
rule getRealignedBam:
    input:
        txt = getRealignedExperiment
    output:
        bam = REALIGNINDELSOUT + '{sample}.bam'
    params:
        originalBam = REALIGNINDELSOUT + 'ORIGINAL_{sample}.bam'
    shell:
        'ln -s {params.originalBam} {output.bam}'

# Rule to create the baserecalibration table used for the baserecalibration
# This is a GATK tool
if not 'BASERECALIBRATIONIN' in globals():
    BASERECALIBRATIONIN = REALIGNINDELSOUT
if not 'BASERECALIBRATIONOUT' in globals():
    BASERECALIBRATIONOUT = OUTDIR + 'recalibratedBases/'

rule firstPassCreateRecalibrationTable:
	input:
		bam = BASERECALIBRATIONIN + '{sample}.bam',
		bai = BASERECALIBRATIONIN + '{sample}.bai',
		reference = {config['resources'][ORGANISM]['reference']},
        	referenceDict = {config['resources'][ORGANISM]['referenceDict']},
        	regions = config['resources'][ORGANISM]['regions']
	output:
        	tab = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table',
	params:
        	lsfoutfile = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.lsfout.log',
        	lsferrfile = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.lsferr.log',
        	scratch = config['tools']['GATK']['baseRecalibrator']['scratch'],
        	mem = config['tools']['GATK']['baseRecalibrator']['mem'],
        	time = config['tools']['GATK']['baseRecalibrator']['time'],
        	known1 = config['tools']['GATK']['baseRecalibrator']['known1'],
        	known2 = config['tools']['GATK']['baseRecalibrator']['known2'],
        	known3 = config['tools']['GATK']['baseRecalibrator']['known3']
	benchmark:
		BASERECALIBRATIONOUT + '{sample}_firstPass_reca.benchmark'
	threads:
		int(config['tools']['GATK']['baseRecalibrator']['threads'])
	shell:
		('{config[tools][GATK][call]} ' +
        	'-T BaseRecalibrator ' +
        	'-L {input.regions} ' +
        	'-R {input.reference} ' +
        	'-I {input.bam} ' +
        	'{params.known1} ' + 
        	'{params.known2} ' + 
        	'{params.known3} ' + 
        	'-nct {threads} ' +
        	'-o {output.tab}')

# Rule to create the base-recalibration table used to analyze the effect of the baserecalibration
# This is a GATK tool
#ruleorder: secondPassCreateRecalibrationTable > firstPassCreateRecalibrationTable
rule secondPassCreateRecalibrationTable:
	input:
		bam = BASERECALIBRATIONOUT + '{sample}.bam',
		bai = BASERECALIBRATIONOUT + '{sample}.bai',
		tab = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table',
		reference = {config['resources'][ORGANISM]['reference']},
        	referenceDict = {config['resources'][ORGANISM]['referenceDict']},
        	regions = config['resources'][ORGANISM]['regions']
	output:
    		tab = BASERECALIBRATIONOUT + '{sample}.secondPass_reca.table',
	params:
        	lsfoutfile = BASERECALIBRATIONOUT + '{sample}_secondPass_reca.lsfout.log',
        	lsferrfile = BASERECALIBRATIONOUT + '{sample}_secondPass_reca.lsferr.log',
        	reference = {config['resources'][ORGANISM]['reference']},
        	scratch = config['tools']['GATK']['baseRecalibrator']['scratch'],
        	mem = config['tools']['GATK']['baseRecalibrator']['mem'],
        	time = config['tools']['GATK']['baseRecalibrator']['time'],
        	known1 = config['tools']['GATK']['baseRecalibrator']['known1'],
        	known2 = config['tools']['GATK']['baseRecalibrator']['known2'],
        	known3 = config['tools']['GATK']['baseRecalibrator']['known3']
	benchmark:
        	BASERECALIBRATIONOUT + '{sample}_secondPass_reca.benchmark'
	threads:
        	int(config['tools']['GATK']['baseRecalibrator']['threads'])
	shell:
        	('{config[tools][GATK][call]} ' +
        	'-T BaseRecalibrator ' +
        	'-L {input.regions} ' +
        	'-R {input.reference} ' +
        	'-I {input.bam} ' +
        	'{params.known1} ' + 
        	'{params.known2} ' + 
        	'{params.known3} ' + 
        	'-nct {threads} ' +
        	'-BQSR {input.tab} ' +
        	'-o {output.tab}')
	
# Rule to realing the reads around indels
# This is a GATK tool
rule baseRecalibration:
	input:
		bam = BASERECALIBRATIONIN + '{sample}.bam',
        	bai = BASERECALIBRATIONIN + '{sample}.bai',
        	tab = BASERECALIBRATIONOUT + '{sample}_firstPass_reca.table',
        	reference = {config['resources'][ORGANISM]['reference']},
        	referenceDict = {config['resources'][ORGANISM]['referenceDict']}
	output:
		bam = temp(BASERECALIBRATIONOUT + '{sample}.bam'),
	params:
		lsfoutfile = BASERECALIBRATIONOUT + '{sample}.lsfout.log',
		lsferrfile = BASERECALIBRATIONOUT + '{sample}.lsferr.log',
        	scratch = config['tools']['GATK']['baseRecalibrator']['scratch'],
        	mem = config['tools']['GATK']['baseRecalibrator']['mem'],
        	time = config['tools']['GATK']['baseRecalibrator']['time']
	benchmark:
		BASERECALIBRATIONOUT + '{sample}.benchmark'
	threads:
		int(config['tools']['GATK']['baseRecalibrator']['threads'])
	shell:
    		('{config[tools][GATK][call]} ' +
        	'-T PrintReads ' +
        	'-R {input.reference} ' +
        	'-I {input.bam} ' +
        	'-o {output.bam} ' +
        	'-BQSR {input.tab}')

if not 'MPILEUPIN' in globals():
    MPILEUPIN = BASERECALIBRATIONOUT
if not 'MPILEUPOUT' in globals():
    MPILEUPOUT = OUTDIR + 'mpileup/'

rule mpileupBcf:
    input:
        bam = MPILEUPIN + '/{sample}.bam',
        reference = config['resources'][ORGANISM]['reference'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        bcf = temp(MPILEUPOUT + '/{sample}.bcf')
    priority: 25
    params:
        lsfoutfile = MPILEUPOUT + '/{sample}.bcf.lsfout.log',
        lsferrfile = MPILEUPOUT + '/{sample}.bcf.lsferr.log'
    threads:
        int(config['tools']['samtools']['mpileup']['threads'])
    benchmark:
        MPILEUPOUT + '/{sample}.bcf.benchmark'
    log:
        MPILEUPOUT + '/{sample}.bcf.log'
    shell:
        ('{config[tools][samtools][call]} ' +
        'mpileup ' +
        '--output-tags DP4,DP -C50 -E -g ' +
        '-f {input.reference} ' +
        '-o {output.bcf} ' +
        '-l {input.regions} ' +
        '{input.bam} ' +
        '2>&1 >{log} ' +
        '&& touch {output.suc}')

rule mpileupMpileup:
    input:
        bam = MPILEUPIN + '{sample}.bam',
        reference = config['resources'][ORGANISM]['reference'],
        regions = config['resources'][ORGANISM]['regions']
    output:
        mpileup = MPILEUPOUT + '{sample}.mpileup'
    params:
        lsfoutfile = MPILEUPOUT + '{sample}.mpileup.lsfout.log',
        lsferrfile = MPILEUPOUT + '{sample}.mpileup.lsferr.log',
        scratch = config['tools']['samtools']['mpileup']['scratch'],
        mem = config['tools']['samtools']['mpileup']['mem'],
        time = config['tools']['samtools']['mpileup']['time'],
        params = config['tools']['samtools']['mpileup']['params']
    threads:
        int(config['tools']['samtools']['mpileup']['threads'])
    benchmark:
        MPILEUPOUT + '{sample}.mpileup.benchmark'
    log:
        MPILEUPOUT + '{sample}.mpileup.log'
    shell:
        ('{config[tools][samtools][call]} mpileup ' +
        '{params.params} ' + 
        '-f {input.reference} ' + 
        '-o {output.mpileup} ' + 
        '-l {input.regions} ' +
        '{input.bam} ')


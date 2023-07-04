import ntpath

# BamClipOverlap - NGSbits
# Softclipping overlapping reads.
rule bamClipOverlap:
	input:
		bam = CLIPOVERLAPIN + '{sample}.bam'
	output:
		bam_clipped = temp(CLIPOVERLAPOUT + '{sample}.bam')
	params:
		lsfoutfile = CLIPOVERLAPOUT + '{sample}.bam.lsfout.log',
        	lsferrfile = CLIPOVERLAPOUT + '{sample}.bam.lsferr.log',
        	scratch = config['tools']['clipoverlap']['scratch'],
        	mem = config['tools']['clipoverlap']['mem'],
        	time = config['tools']['clipoverlap']['time'],
                params = config['tools']['clipoverlap']['params']
	threads:
		config['tools']['clipoverlap']['threads']
	benchmark:
		CLIPOVERLAPOUT + '{sample}.clipoverlap.benchmark'
	shell:
		'{config[tools][clipoverlap][call]} -in {input.bam} -out {output.bam_clipped} {params.params}'


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
        	config['tools']['picard']['SortSam']['threads']
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
        config['tools']['picard']['SortSam']['threads']
    shell:
        ('{config[tools][picard][call]} SortSam ' +
        'INPUT={input.bam} ' +
        'OUTPUT={output.bam} ' +
        'SORT_ORDER= coordinate ' +
	'VALIDATION_STRINGENCY=LENIENT ' +
        'MAX_RECORDS_IN_RAM={params.max_records_in_ram} ' +
        'TMP_DIR={TMPDIR} ')



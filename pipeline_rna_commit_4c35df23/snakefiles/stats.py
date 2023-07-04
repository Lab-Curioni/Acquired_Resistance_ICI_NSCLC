###
#FastQC
###

rule fastqc:
		input:
			fastq =  FASTQDIR + '{sample}.fastq.gz'
		output:
			out = STATSOUT + '{sample}_fastqc.zip'
		params:
			out = STATSOUT,
			lsfoutfile = STATSOUT + '{sample}.fastqc.lsfout.log',
			lsferrfile = STATSOUT + '{sample}.fastqc.lsferr.log',
			scratch = config['tools']['fastqc']['scratch'],
			mem = config['tools']['fastqc']['mem'],
			time = config['tools']['fastqc']['time']
		benchmark:
			STATSOUT + '{sample}.fastqc.benchmark'
		shell:
			('{config[tools][fastqc][call]} -o {params.out} {input.fastq}')

###
#Qorts
###
if PAIREDEND == True:
	rule qorts:
		input:
			bam = STATSIN + '{sample}_Aligned.out.bam'
		output:
			out = STATSOUT + '{sample}.qorts/{sample}_QC.multiPlot.png'
#			out = STATSOUT + '{sample}.qorts.success.txt'
		params:
			lsfoutfile = STATSOUT + '{sample}.qorts.lsfout.log',
	                lsferrfile = STATSOUT + '{sample}.qorts.lsferr.log',
	                scratch = config['tools']['qorts']['scratch'],
	                mem = config['tools']['qorts']['mem'],
	                time = config['tools']['qorts']['time'],
			gtf = config['resources'][ORGANISM]['annotation'],
			outdir = STATSOUT +'{sample}.qorts/',
			strandedness = config['general']['strandedness_qorts'],
			name = '{sample}'
		benchmark:
	                STATSOUT + '{sample}.qorts.benchmark'
		shell:
			('{config[tools][samtools][call]} view -b -F4 -F8 {input.bam} | {config[tools][qorts][call]} QC --runFunctions writeGeneBody,writeGenewiseGeneBody --generatePdfReport --outfilePrefix {params.name}. --generatePlots {params.strandedness} - {params.gtf} {params.outdir}; cp {params.outdir}QC.multiPlot.png {params.outdir}{params.name}.QC.multiPlot.png')
else:
	rule qorts:
                input:
                        bam = STATSIN + '{sample}_Aligned.out.bam'
                output:
                        out = STATSOUT + '{sample}.qorts/{sample}_QC.multiPlot.png'
#                        out = STATSOUT + '{sample}.qorts/QC.summary.txt'
                params:
                        lsfoutfile = STATSOUT + '{sample}.qorts.lsfout.log',
                        lsferrfile = STATSOUT + '{sample}.qorts.lsferr.log',
                        scratch = config['tools']['qorts']['scratch'],
                        mem = config['tools']['qorts']['mem'],
                        time = config['tools']['qorts']['time'],
                        gtf = config['resources'][ORGANISM]['annotation'],
                        outdir = STATSOUT +'{sample}.qorts/',
			strandedness = config['general']['strandedness_qorts'],
			name = '{sample}'
                benchmark:
                        STATSOUT + '{sample}.qorts.benchmark'
                shell:
                        ('{config[tools][samtools][call]} view -b -F4 -F8 {input.bam} | {config[tools][qorts][call]} QC --singleEnded --runFunctions writeGeneBody,writeGenewiseGeneBody --generatePdfReport --outfilePrefix {params.name}. --generatePlots {params.strandedness} - {params.gtf} {params.outdir}; cp {params.outdir}QC.multiPlot.png {params.outdir}{params.name}.QC.multiPlot.png')

###
#Htseq
###

rule HTseq:
        input:
                bam = STATSIN + '{sample}_Aligned.out.bam',
                bai = STATSIN + '{sample}_Aligned.out.bam.bai'
        output:
                txt = HTSEQOUT + '{sample}.htseq_counts.txt',
        params:
                lsfoutfile = HTSEQOUT + '{sample}.htseq.lsfout.log',
                lsferrfile = HTSEQOUT + '{sample}.htseq.lsferr.log',
                scratch = config['tools']['htseq']['scratch'],
                mem = config['tools']['htseq']['mem'],
                time = config['tools']['htseq']['time'],
		strandedness = config['general']['strandedness'],
		annotation = config['resources'][ORGANISM]['annotation']
        benchmark:
                HTSEQOUT + '{sample}.htseq.benchmark'
        threads:
                int(config['tools']['htseq']['threads'])
        shell:
                ('{config[tools][samtools][call]} view -F4 {input.bam} | {config[tools][htseq][call]} -m intersection-nonempty --stranded={params.strandedness} --idattr gene_id - {params.annotation} > {output.txt}'
                )


if not 'GETHGNCSYMBOLSIN' in globals():
    GETHGNCSYMBOLSIN = HTSEQOUT
if not 'GETHGNCSYMBOLSOUT' in globals():
    GETHGNCSYMBOLSOUT = HTSEQOUT

# This rule adds to the output of HTseq, which has to columns: ensemblID and count, the hgnc symbols of the genes

rule getHgncSymbols:
    input:
        htseq = GETHGNCSYMBOLSIN + '{sample}.htseq_counts.txt'
    output:
        out = GETHGNCSYMBOLSOUT + '{sample}.htseq_counts_hgnc.txt'
    params:
        mem = config['tools']['getHgncSymbols']['mem'],
        time = config['tools']['getHgncSymbols']['time'],
        scratch = config['tools']['getHgncSymbols']['scratch'],
        lsferrfile = GETHGNCSYMBOLSOUT + "/log/{sample}.getHgncSymbols.lsferr.log",
        lsfoutfile = GETHGNCSYMBOLSOUT + "/log/{sample}.getHgncSymbols.lsfout.log",
    benchmark:
        GETHGNCSYMBOLSOUT + '/log/{sample}.getHgncSymbols.benchmark'
    threads:
        config['tools']['getHgncSymbols']['threads']
    shell:
        '{config[tools][getHgncSymbols][call]} --inFile {input.htseq} --outFile {output.out}'


if not 'FILTERGENESIN' in globals():
    FILTERGENESIN = HTSEQOUT
if not 'FILTERGENESOUT' in globals():
    FILTERGENESOUT = HTSEQOUT

# This rule filters the output table of HTseq. It returns only those rows where the genes are contained in a list that is given to this rule, the list normally contains clinically relevant genes

rule filter_genesOfInterest:
    input:
        finalTable = FILTERGENESIN + '{sample}.htseq_counts_hgnc.txt'
    output:
        out = FILTERGENESOUT + '{sample}.htseq_counts_hgnc_filtered.txt'
    params:
        mem = config['tools']['filter_genesOfInterest']['mem'],
        time = config['tools']['filter_genesOfInterest']['time'],
        scratch = config['tools']['filter_genesOfInterest']['scratch'],
        lsferrfile = FILTERGENESOUT + "/log/{sample}.filter_genesOfInterest.lsferr.log",
        lsfoutfile = FILTERGENESOUT + "/log/{sample}.filter_genesOfInterest.lsfout.log",
        genelist = config['tools']['filter_genesOfInterest']['genelist']
    benchmark:
        FILTERGENESOUT + '/log/{sample}.filter_genesOfInterest.benchmark'
    threads:
        config['tools']['filter_genesOfInterest']['threads']
    shell:
        '{config[tools][filter_genesOfInterest][call]} --inFile {input.finalTable} --outFile {output.out} --geneList {params.genelist}'


###
#RSeQC
###

rule bam_stat:
        input:
                bam = STATSIN + '{sample}_Aligned.out.bam',
                bai = STATSIN + '{sample}_Aligned.out.bam.bai'
        output:
                txt = STATSOUT + '{sample}.rseqc.bam_stat.txt',
        params:
                lsfoutfile = STATSOUT + '{sample}.rseqc.bamstat.lsfout.log',
                lsferrfile = STATSOUT + '{sample}.rseqc.bamstat.lsferr.log',
                scratch = config['tools']['rseqc']['bamstat']['scratch'],
                mem = config['tools']['rseqc']['bamstat']['mem'],
                time = config['tools']['rseqc']['bamstat']['time']
        benchmark:
                STATSOUT + '{sample}.rseqc.bamstat.benchmark'
        threads:
                int(config['tools']['rseqc']['bamstat']['threads'])
        shell:
                ('{config[tools][rseqc][bamstat][call]} -i {input.bam} 2> {output.txt}'
                )

rule read_GC:
        input:
                bam = STATSIN + '{sample}_Aligned.out.bam',
                bai = STATSIN + '{sample}_Aligned.out.bam.bai'
        output:
                txt = STATSOUT + '{sample}.rseqc.read_gc.GC_plot.pdf',
        params:
                lsfoutfile = STATSOUT + '{sample}.rseqc.read_gc.lsfout.log',
                lsferrfile = STATSOUT + '{sample}.rseqc.read_gc.lsferr.log',
                scratch = config['tools']['rseqc']['read_gc']['scratch'],
                mem = config['tools']['rseqc']['read_gc']['mem'],
                time = config['tools']['rseqc']['read_gc']['time'],
		txt = STATSOUT + '{sample}.rseqc.read_gc'
        benchmark:
                STATSOUT + '{sample}.rseqc.read_gc.benchmark'
        threads:
                int(config['tools']['rseqc']['read_gc']['threads'])
        shell:
                ('{config[tools][rseqc][read_gc][call]} -i {input.bam} -o {params.txt}'
                )

rule geneBodyCoverage:
        input:
                bam = STATSIN + '{sample}_Aligned.sortedByCoord.out.bam',
                bai = STATSIN + '{sample}_Aligned.sortedByCoord.out.bam.bai'
        output:
                txt = STATSOUT + '{sample}.rseqc.geneBodyCoverage.txt',
        params:
                lsfoutfile = STATSOUT + '{sample}.geneBodyCoverage.lsfout.log',
                lsferrfile = STATSOUT + '{sample}.geneBodyCoverage.lsferr.log',
                scratch = config['tools']['rseqc']['geneBodyCoverage']['scratch'],
                mem = config['tools']['rseqc']['geneBodyCoverage']['mem'],
                time = config['tools']['rseqc']['geneBodyCoverage']['time'],
                txt = STATSOUT + '{sample}.rseqc',
		ref = config['resources'][ORGANISM]['refgene']
        benchmark:
                STATSOUT + '{sample}.rseqc.geneBodyCoverage.benchmark'
        threads:
                int(config['tools']['rseqc']['geneBodyCoverage']['threads'])
        shell:
                ('{config[tools][rseqc][geneBodyCoverage][call]} -i {input.bam} -o {params.txt} -r {params.ref}'
                )

###
#Picard CollectRnaSeqMetrics
###

rule CollectRnaSeqMetrics:
	input:
		bam = STATSIN + '{sample}_Aligned.out.bam',
		bai = STATSIN + '{sample}_Aligned.out.bam.bai'
	output:
		txt = STATSOUT + '{sample}.RNAMetrics.txt',
	params:
		ref_flat = config['resources'][ORGANISM]['ref_flat'],
		strand = config['general']['picard_strand'],
		lsfoutfile = STATSOUT + '{sample}.CollectRnaSeqMetrics.lsfout.log',
	        lsferrfile = STATSOUT + '{sample}.CollectRnaSeqMetrics.lsferr.log',
	        scratch = config['tools']['picard']['CollectRnaSeqMetrics']['scratch'],
	        mem = config['tools']['picard']['CollectRnaSeqMetrics']['mem'],
	        time = config['tools']['picard']['CollectRnaSeqMetrics']['time']
	benchmark:
		STATSOUT + '{sample}.CollectRnaSeqMetrics.benchmark'
	threads:
		int(config['tools']['picard']['CollectRnaSeqMetrics']['threads'])
	shell:
		('{config[tools][picard][call]} CollectRnaSeqMetrics ' + 
		'INPUT={input.bam} ' +
		'REF_FLAT={params.ref_flat} '+
		'OUTPUT={output.txt} '+
		'VALIDATION_STRINGENCY=LENIENT ' +
		'STRAND_SPECIFICITY={params.strand}'
		)



	
###
#Samtools index
###

rule createIndex:
    input:
        bam = '{samples}_Aligned.sortedByCoord.out.bam',
    output:
        idx = '{samples}_Aligned.sortedByCoord.out.bam.bai',
	idx2 = '{samples}_Aligned.out.bam.bai'
    params:
        lsfoutfile = '{samples}.bai.lsfout.log',
        lsferrfile = '{samples}.bai.lsferr.log',
        scratch = config['tools']['samtools']['index']['scratch'],
        mem = config['tools']['samtools']['index']['mem'],
        time = config['tools']['samtools']['index']['time'],
	sample = '{samples}'
    threads:
        int(config['tools']['samtools']['index']['threads'])
    benchmark:
        '{samples}.bai.benchmark'
    shell:
        ('{config[tools][samtools][call]} index {input.bam} ; cp {params.sample}_Aligned.sortedByCoord.out.bam.bai {params.sample}_Aligned.out.bam.bai')




import os, glob, sys
include: "./misc_snake_neu.py"


# Check if the uses specified the proper input and output directories
if not 'FASTQDIR' in globals():
    print('You have to specify the root directory of the fastq files!')
    sys.exit(1)
if not 'OUTDIR' in globals():
    print('You have to specify the root directory where the results will be generated!')
    sys.exit(1)
if not 'PAIREDEND' in globals():
    print('You have to specify if you have Pairedend-data or not!')
    sys.exit(1)

if not 'STARALIGNOUT' in globals():
    STARALIGNOUT = OUTDIR + 'alignment/'
if not 'STATSIN' in globals():
    STATSIN = STARALIGNOUT
if not 'STATSOUT' in globals():
    STATSOUT = OUTDIR + 'stats/'
if not 'HTSEQIN' in globals():
    HTSEQIN = STARALIGNOUT
if not 'HTSEQOUT' in globals():
    HTSEQOUT = OUTDIR + 'gene_counts/'
if not 'WARNINGSIN' in globals():
    WARNINGSIN = STATSOUT
if not 'WARNINGSOUT' in globals():
    WARNINGSOUT = OUTDIR + 'warnings/'

include: './align.py'
include: './stats.py'
include: './warnings.py'
include: './fusion_detection_rules.py'

SAMPLES=getSampleNames()
MATEPAIR=[]
if PAIREDEND==True:
	MATEPAIR=['_R1', '_R2']
else:
	MATEPAIR=['_R1']

rule rnaseq:
	input:
		expand(STARALIGNOUT + '{samples}_Aligned.out.bam', samples=SAMPLES),
		expand(STATSIN + '{samples}_Aligned.sortedByCoord.out.bam.bai', samples=SAMPLES),
		expand(STATSOUT + '{samples}.RNAMetrics.txt', samples=SAMPLES),
		expand(STATSOUT + '{samples}.rseqc.bam_stat.txt', samples=SAMPLES),
		expand(STATSOUT + '{samples}.rseqc.read_gc.GC_plot.pdf', samples=SAMPLES),
	        expand(HTSEQOUT + '{samples}.htseq_counts.txt', samples=SAMPLES),
		expand(STATSOUT + '{samples}.rseqc.geneBodyCoverage.txt', samples=SAMPLES),
#		expand(HTSEQOUT + '{samples}.htseq_counts_hgnc.txt', samples=SAMPLES),
#		expand(HTSEQOUT + '{samples}.htseq_counts_hgnc_filtered.txt', samples=SAMPLES),
#		expand(STATSOUT + '{samples}.qorts/QC.summary.txt', samples=SAMPLES),
#		expand(STATSOUT + '{samples}.qorts.success.txt', samples=SAMPLES),
		expand(STATSOUT + '{samples}.qorts/{samples}_QC.multiPlot.png', samples=SAMPLES),
		expand(STATSOUT + '{samples}{mate}_fastqc.zip', samples=SAMPLES, mate=MATEPAIR),
		expand(WARNINGSOUT + '{samples}{mate}_fastqc.warning.txt', samples = SAMPLES, mate=MATEPAIR),
###		expand(FUSIONCATCHEROUT + '{samples}_final-list_candidate-fusion-genes.GRCh37.txt', samples=SAMPLES),
		expand(FUSIONCATCHEROUT + '{samples}/{samples}_fusioncatcher.done', samples=SAMPLES),
		expand(STARFUSIONOUT + '{samples}/{samples}.star-fusion.done', samples=SAMPLES),
###		expand(ARRIBAOUT + '{samples}_fusions.tsv', samples=SAMPLES),
		expand(ARRIBAOUT + '{samples}/{samples}_arriba.done', samples=SAMPLES)
	output:
		OUTDIR + 'icgc_complete.txt'
	params:
		lsfoutfile = OUTDIR + 'complete.lsfout.log',
		lsferrfile = OUTDIR + 'complete.lsferr.log',
		out = STATSOUT,
		indir = OUTDIR,
		samples = SAMPLEMAPPING,
		plot_genebodycoverage = config['tools']['plots_stats']['genebody_coverage'],
		plot_gc = config['tools']['plots_stats']['gc_plot'],
		generate_readStats = config['tools']['plots_stats']['generate_readStats'],
		qorts = config['tools']['qorts']['plot'],
	shell:
		'R --vanilla --args {params.out} < {params.plot_gc} > {params.out}out_gc_plot.txt &' +
		'python {params.generate_readStats} {params.indir} {params.samples} &'+
		'R --vanilla --args {params.out} < {params.plot_genebodycoverage} > {params.out}out_genebody_coverage.txt ;'+
		'echo -e "unique.ID\tsample.ID\tgroup.ID\tqc.data.dir" >{params.out}decoder.txt;' +
		'while read sample; do echo -e "$sample.fastq.gz\t$sample.fastq.gz\t$sample\t{params.out}/$sample.qorts">>{params.out}decoder.txt; done <{params.samples};' +
		'Rscript {params.qorts} {params.out}decoder.txt {params.out}qorts.genebody.pdf ; date > {output}'

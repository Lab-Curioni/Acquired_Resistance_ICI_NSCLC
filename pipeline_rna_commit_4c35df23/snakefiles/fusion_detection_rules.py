if not 'FUSIONCATCHERIN' in globals():
    FUSIONCATCHERIN = FASTQDIR
if not 'FUSIONCATCHEROUT' in globals():
    FUSIONCATCHEROUT = OUTDIR + 'fusions/fusioncatcher/'

# This rule runs the gene fusion detection pipeline FusionCatcher from the start
rule fusioncatcher:
    input:
        fastq = FUSIONCATCHERIN + '{sample}_R1.fastq.gz',
        fastq2 = FUSIONCATCHERIN + '{sample}_R2.fastq.gz'
    output:
        #hg19_list = FUSIONCATCHEROUT + 'final-list_candidate-fusion-genes.GRCh37.txt',
        done = FUSIONCATCHEROUT + '{sample}/{sample}_fusioncatcher.done'
    params:
        mem = config['tools']['fusioncatcher']['mem'],
        time = config['tools']['fusioncatcher']['time'],
        scratch = config['tools']['fusioncatcher']['scratch'],
        lsferrfile = FUSIONCATCHEROUT + "log/{sample}.fusioncatcher.lsferr.log",
        lsfoutfile = FUSIONCATCHEROUT + "log/{sample}.fusioncatcher.lsfout.log",
        outdir = FUSIONCATCHEROUT + "/{sample}" ,
        sample = '{sample}',
        tmpdir = FUSIONCATCHEROUT + 'fusioncatcher_tmp/',
        configfile = config['tools']['fusioncatcher']['configfile']
    benchmark:
        FUSIONCATCHEROUT + 'log/{sample}.fusioncatcher.benchmark'
    threads:
        config['tools']['fusioncatcher']['threads']
    shell:
        '{config[tools][fusioncatcher][call]} -i {input.fastq},{input.fastq2} -o {params.outdir} -p {threads} --config {params.configfile} --tmp {params.tmpdir}; touch {params.outdir}/{params.sample}_fusioncatcher.done'
# cp {params.outdir}/final-list_candidate-fusion-genes.GRCh37.txt {params.outdir}/{params.sample}_final-list_candidate-fusion-genes.GRCh37.txt; cp {params.outdir}/summary_candidate_fusions.txt {params.outdir}/{params.sample}_summary_candidate_fusions.txt; cp {params.outdir}/final-list_candidate-fusion-genes.txt {params.outdir}/{params.sample}_final-list_candidate-fusion-genes.txt'


if not 'ARRIBAIN' in globals():
    ARRIBAIN = FASTQDIR
if not 'ARRIBAOUT' in globals():
    ARRIBAOUT = OUTDIR + 'fusions/arriba/'

# This rule runs the gene fusion detection tool Arriba from the start
rule arriba:
    input:
        fastq = ARRIBAIN + '{sample}_R1.fastq.gz',
        fastq2 = ARRIBAIN + '{sample}_R2.fastq.gz'
    output:
        #outtable = ARRIBAOUT + 'fusions.tsv',
        done = ARRIBAOUT + '{sample}/{sample}_arriba.done'
    params:
        mem = config['tools']['arriba']['mem'],
        time = config['tools']['arriba']['time'],
        scratch = config['tools']['arriba']['scratch'],
        lsferrfile = ARRIBAOUT + "log/{sample}.arriba.lsferr.log",
        lsfoutfile = ARRIBAOUT + "log/{sample}.arriba.lsfout.log",
        star_genomedir = config['resources'][ORGANISM]['genomedir'],
        annotation = config['resources'][ORGANISM]['annotation'],
        assembly = config['resources'][ORGANISM]['genomefasta'],
        blacklist = config['tools']['arriba']['blacklist'],
        sample = '{sample}',
        outdir = ARRIBAOUT + '/{sample}'
    benchmark:
        ARRIBAOUT + 'log/{sample}.arriba.benchmark'
    threads:
        config['tools']['arriba']['threads']
    shell:
        'cd {params.outdir}/; {config[tools][arriba][call]} {params.star_genomedir} {params.annotation} {params.assembly} {params.blacklist} {input.fastq} {input.fastq2} {threads}; touch {params.outdir}/{params.sample}_arriba.done'
#; cp {params.outdir}/fusions.tsv {params.outdir}/{params.sample}.fusions.tsv; cp {params.outdir}/fusions.discarded.tsv {params.outdir}/{params.sample}.fusions.discarded.tsv'


if not 'STARFUSIONIN' in globals():
    STARFUSIONIN = FASTQDIR
if not 'STARFUSIONOUT' in globals():
    STARFUSIONOUT = OUTDIR + 'fusions/starfusion/'

# This rule runs the gene fusion prediction tool STAR-fusion from the start
rule star_fusion:
    input:
        fastq = STARFUSIONIN + '{sample}_R1.fastq.gz',
        fastq2 = STARFUSIONIN + '{sample}_R2.fastq.gz'
    output:
#        outtable = STARFUSIONOUT + '{sample}.star-fusion.fusion_predictions.tsv',
#        outtable_abridged = STARFUSIONOUT + '{sample}.star-fusion.fusion_predictions.abridged.tsv'
        done = STARFUSIONOUT + '{sample}/{sample}.star-fusion.done'
    params:
        mem = config['tools']['star_fusion']['mem'],
        time = config['tools']['star_fusion']['time'],
        scratch = config['tools']['star_fusion']['scratch'],
        lsferrfile = STARFUSIONOUT + "log/{sample}.star_fusion.lsferr.log",
        lsfoutfile = STARFUSIONOUT + "log/{sample}.star_fusion.lsfout.log",
        star_genomedir = config['tools']['star_fusion']['star_genomedir'],
        outdir = STARFUSIONOUT + "/{sample}",
        sample = '{sample}'
    benchmark:
        STARFUSIONOUT + 'log/{sample}.star_fusion.benchmark'
    threads:
        config['tools']['star_fusion']['threads']
    shell:
        '{config[tools][star_fusion][call]} --left_fq {input.fastq} --right_fq {input.fastq2} --genome_lib_dir {params.star_genomedir} --CPU {threads} --output_dir {params.outdir}; touch {params.outdir}/{params.sample}_star-fusion.done'

#; cp {params.outdir}/star-fusion.fusion_predictions.abridged.tsv {params.outdir}/{params.sample}.star-fusion.fusion_predictions.abridged.tsv; cp {params.outdir}/star-fusion.fusion_predictions.tsv {params.outdir}/{params.sample}.star-fusion.fusion_predictions.tsv'



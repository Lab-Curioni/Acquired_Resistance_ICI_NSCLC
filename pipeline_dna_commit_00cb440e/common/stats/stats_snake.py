# This rule uses Picard Calculate HsMetrics to calculate statistics of a BAM file.

rule hsmetrics:
    input:
        bam = '{sample}.bam',
        ref = config['resources'][ORGANISM]['reference'],
        bait = config['resources'][ORGANISM]['intervalsHsMetrics'],
        target = config['resources'][ORGANISM]['intervalsHsMetrics']
    output:
        file = '{sample}.hsmetrics.txt'
    priority:
        50
    params:
        lsfoutfile = '{sample}.hsmetrics.lsfout.log',
        lsferrfile = '{sample}.hsmetrics.lsferr.log',
        scratch = config['tools']['picard']['calculateHsMetrics']['scratch'],
        mem = config['tools']['picard']['calculateHsMetrics']['mem'],
        time = config['tools']['picard']['calculateHsMetrics']['time']
    benchmark:
        '{sample}.hsmetrics.benchmark'
    threads:
        config['tools']['picard']['calculateHsMetrics']['threads']
    shell:
        '{config[tools][picard][call]} CalculateHsMetrics I={input.bam} O={output.file} R={input.ref} BAIT_INTERVALS={input.bait} TARGET_INTERVALS={input.target}'



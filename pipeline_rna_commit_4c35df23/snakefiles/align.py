if not 'FASTQIN' in globals():
    FASTQIN = FASTQDIR
if not 'STARALIGNOUT' in globals():
    STARALIGNOUT = OUTDIR + 'star/'


if PAIREDEND == False:
    rule align_star:
        input:
            fastq = FASTQIN + '{fastq}_R1.fastq.gz',
        output:
            bam = STARALIGNOUT + '{fastq}_Aligned.sortedByCoord.out.bam',
            outSJ = STARALIGNOUT + "{fastq}_Aligned.out.bam"
        params:
            mem = config['tools']['STAR']['mem'],
            time = config['tools']['STAR']['time'],
            scratch = config['tools']['STAR']['scratch'],
            lsferrfile = STARALIGNOUT + '{fastq}_align_star.lsferr.log',
            lsfoutfile = STARALIGNOUT + '{fastq}_align_star.lsfout.log',
            genomedir = config['resources'][ORGANISM]['genomedir'],
            variousParams = config['tools']['STAR']['variousParams'],
            out = STARALIGNOUT + '{fastq}_',
            fastq = '{fastq}',
            Star_sjdbOverhang = config['general']['Star_sjdbOverhang']
        benchmark:
            STARALIGNOUT + '{fastq}_align_star.benchmark'
        threads:
            config['tools']['STAR']['threads']
        shell:
            ('{config[tools][STAR][call]} ' +
            '--genomeDir {params.genomedir} ' +
            '--readFilesIn {input.fastq} ' +
            '--runThreadN {threads} ' +
            '--outFileNamePrefix {params.out} ' +
            '--sjdbOverhang {params.Star_sjdbOverhang} ' +
            '{params.variousParams}')

    rule align_star_trans:
        input:
            fastq = FASTQIN + '{fastq}_R1.fastq.gz',
        output:
            trans_out = STARALIGNOUT + '{fastq}_Aligned.toTranscriptome.out.bam'
        params:
            mem = config['tools']['STAR']['mem'],
            time = config['tools']['STAR']['time'],
            scratch = config['tools']['STAR']['scratch'],
            lsferrfile = STARALIGNOUT + '{fastq}_align_star_trans.lsferr.log',
            lsfoutfile = STARALIGNOUT + '{fastq}_align_star_trans.lsfout.log',
            genomedir = config['resources'][ORGANISM]['genomedir'],
            variousParams = config['tools']['STAR']['variousParams'],
            out = STARALIGNOUT + '{fastq}_',
            fastq = '{fastq}',
            Star_sjdbOverhang = config['general']['Star_sjdbOverhang']
        benchmark:
            STARALIGNOUT + '{fastq}_align_star_trans.benchmark'
        threads:
            config['tools']['STAR']['threads']
        shell:
            ('{config[tools][STAR][call]} ' +
            '--genomeDir {params.genomedir} ' +
            '--readFilesIn {input.fastq} ' +
            '--runThreadN {threads} ' +
            '--outFileNamePrefix {params.out} ' +
            '--sjdbOverhang {params.Star_sjdbOverhang} ' +
            '{params.variousParams}')


elif PAIREDEND == True:
    rule align_star_pe:
            input:
                fastq = FASTQIN + '{fastq}_R1.fastq.gz',
                fastq2 = FASTQIN + '{fastq}_R2.fastq.gz'
            output:
                bam = STARALIGNOUT + '{fastq}_Aligned.sortedByCoord.out.bam',
                outSJ = STARALIGNOUT + '{fastq}_Aligned.out.bam',
            params:
                mem = config['tools']['STAR']['mem'],
                time = config['tools']['STAR']['time'],
                scratch = config['tools']['STAR']['scratch'],
                lsferrfile = STARALIGNOUT + '{fastq}_align_star_pe.lsferr.log',
                lsfoutfile = STARALIGNOUT + '{fastq}_align_star_pe.lsfout.log',
                genomedir = config['resources'][ORGANISM]['genomedir'],
                variousParams = config['tools']['STAR']['variousParams'],
                out = STARALIGNOUT + '{fastq}_',
                fastq = '{fastq}',
                Star_sjdbOverhang = config['general']['Star_sjdbOverhang']
            benchmark:
                STARALIGNOUT + '{fastq}_align_star_pe.benchmark'
            threads:
                config['tools']['STAR']['threads']
            shell:
                ('{config[tools][STAR][call]} ' +
                '--genomeDir {params.genomedir} ' +
                '--readFilesIn {input.fastq} {input.fastq2} ' +
                '--runThreadN {threads} ' +
                '--outFileNamePrefix {params.out} ' +
                '--sjdbOverhang {params.Star_sjdbOverhang} ' +
                '{params.variousParams}')

    rule align_star_trans_pe:
            input:
                fastq = FASTQIN + '{fastq}_R1.fastq.gz',
                fastq2 = FASTQIN + '{fastq}_R2.fastq.gz'
            output:
                trans_out = STARALIGNOUT + '{fastq}_Aligned.toTranscriptome.out.bam'
            params:
                mem = config['tools']['STAR']['mem'],
                time = config['tools']['STAR']['time'],
                scratch = config['tools']['STAR']['scratch'],
                lsferrfile = STARALIGNOUT + '{fastq}_align_star_trans_pe.lsferr.log',
                lsfoutfile = STARALIGNOUT + '{fastq}_align_star_trans_pe.lsfout.log',
                genomedir = config['resources'][ORGANISM]['genomedir'],
                variousParams = config['tools']['STAR']['variousParams'],
                out = STARALIGNOUT + '{fastq}_',
                fastq = '{fastq}',
                Star_sjdbOverhang = config['general']['Star_sjdbOverhang']
            benchmark:
                STARALIGNOUT + '{fastq}_align_star_trans_pe.benchmark'
            threads:
                config['tools']['STAR']['threads']
            shell:
                ('{config[tools][STAR][call]} ' +
                '--genomeDir {params.genomedir} ' +
                '--readFilesIn {input.fastq} {input.fastq2} ' +
                '--runThreadN {threads} ' +
                '--outFileNamePrefix {params.out} ' +
                '--sjdbOverhang {params.Star_sjdbOverhang} ' +
                '{params.variousParams}')

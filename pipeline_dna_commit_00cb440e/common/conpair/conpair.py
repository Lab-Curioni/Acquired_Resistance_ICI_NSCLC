rule gatk_pileup:
    input:
        bam = CONPAIRIN+ '{sample}.bam',
    	bai = CONPAIRIN+ '{sample}.bai'
    output:
        pileup=temp(CONPAIROUT + '{sample}_pileup'),
    params:
        lsfoutfile = CONPAIROUT + '{sample}_pileup.lsfout.log',
        lsferrfile = CONPAIROUT + '{sample}_pileup.lsferr.log',
        scratch = config['tools']['conpair']['scratch'],
        mem = config['tools']['conpair']['mem'],
        time = config['tools']['conpair']['time'],
	GATK = config['tools']['conpair']['GATK'], 
	marker = config['resources']['H_sapiens_hg19']['marker'],
	ref =  config['resources']['H_sapiens_hg19']['reference']
    benchmark:
        CONPAIROUT + '{sample}_pileup.benchmark'
    threads:
        int(config['tools']['conpair']['threads'])
    shell:
        'export CONPAIR_DIR={config[tools][conpair][conpair_dir]} ; {config[tools][conpair][conpair_dir]}/scripts/run_gatk_pileup_for_sample.py -G {params.GATK} --markers {params.marker} --reference {params.ref} -B {input.bam} -O {output.pileup}'

rule concordance:
	input:
		pileT = CONPAIROUT + '{tumor}_pileup',
		pileN = CONPAIROUT + '{normal}_pileup'
	output:
		concordance = CONPAIROUT + '{tumor}_vs_{normal}_concordance.txt'
	params:
		lsfoutfile = CONPAIROUT + '{tumor}_{normal}_concordance.lsfout.log',
        	lsferrfile = CONPAIROUT + '{tumor}_{normal}_concordance.lsferr.log',
        	scratch = config['tools']['conpair']['scratch'],
        	mem = config['tools']['conpair']['mem'],
        	time = config['tools']['conpair']['time'],
		markers = config['resources']['H_sapiens_hg19']['marker_txt']
	benchmark:
		CONPAIROUT + '{tumor}_{normal}_concordance.benchmark'
	threads:
		int(config['tools']['conpair']['threads'])
	shell:
		'export CONPAIR_DIR={config[tools][conpair][conpair_dir]} ; export PYTHONPATH=${{PYTHONPATH}}:{config[tools][conpair][conpair_dir]}/modules ; {config[tools][conpair][conpair_dir]}/scripts/verify_concordance.py --markers {params.markers} -T {input.pileT} -N {input.pileN} -H --outfile {output.concordance}'	

rule contamination:
        input:
                pileT = CONPAIROUT + '{tumor}_pileup',
                pileN = CONPAIROUT + '{normal}_pileup'
        output:
                contamination = CONPAIROUT + '{tumor}_vs_{normal}_contamination.txt'
        params:
                lsfoutfile = CONPAIROUT + '{tumor}_{normal}_contamination.lsfout.log',
                lsferrfile = CONPAIROUT + '{tumor}_{normal}_contamination.lsferr.log',
                scratch = config['tools']['conpair']['scratch'],
                mem = config['tools']['conpair']['mem'],
                time = config['tools']['conpair']['time'],
                markers = config['resources']['H_sapiens_hg19']['marker_txt']
        benchmark:
                CONPAIROUT + '{tumor}_{normal}_contamination.benchmark'
        threads:
                int(config['tools']['conpair']['threads'])
        shell:
                'export CONPAIR_DIR={config[tools][conpair][conpair_dir]} ; export PYTHONPATH=${{PYTHONPATH}}:{config[tools][conpair][conpair_dir]}/modules ; {config[tools][conpair][conpair_dir]}/scripts/estimate_tumor_normal_contamination.py --markers {params.markers} -T {input.pileT} -N {input.pileN} --outfile {output.contamination}'

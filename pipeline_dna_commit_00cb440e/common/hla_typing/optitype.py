
rule optitype:
	input:
		forward = HLATYPEIN + '{normal}/PAIREDEND/{normal}_R1.fastq.gz',
		reverse = HLATYPEIN + '{normal}/PAIREDEND/{normal}_R2.fastq.gz'
	output:
		outcov = HLATYPEOUT + '{normal}/{normal}_coverage_plot.pdf',
		outres = HLATYPEOUT + '{normal}/{normal}_result.tsv'
	params:
		mem = config['tools']['optitype']['mem'],
		time = config['tools']['optitype']['time'],
		threads = config['tools']['optitype']['threads'],
		scratch = config['tools']['optitype']['scratch'],
		lsfoutfile = HLATYPEOUT + 'optitype.lsfout.log',
		lsferrfile = HLATYPEOUT + 'optitype.lsferr.log',
		out = HLATYPEOUT,
		normal = '{normal}',
		outdir = HLATYPEOUT + '{normal}/'
	shell:
		"""
		cd {params.out}
		mkdir {params.normal} || true
		{config[tools][optitype][call]} {input.forward} {input.reverse} {config[tools][optitype][config_optitype]} {params.outdir}
		cd {params.outdir}
		mv */*_coverage_plot.pdf ./{params.normal}_coverage_plot.pdf || true
        mv */*_result.tsv ./{params.normal}_result.tsv || true
		"""

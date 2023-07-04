###
# STAR Warning
# Produces a warning if the number of reads is < 1M and or less then 50% of the reads are aligned.
###
"""
rule star_warning:
		input:
			alignment = SECONDPASSOUT + '{sample}_Log.final.out'
		output:
			out = WARNINGSOUT + '{sample}_star.warning.txt'
		params:
			lsfoutfile = WARNINGSOUT + '{sample}_star.warning.lsfout.log',
			lsferrfile = WARNINGSOUT + '{sample}_star.warning.lsferr.log',
			scratch = config['tools']['warnings']['scratch'],
			mem = config['tools']['warnings']['mem'],
			time = config['tools']['warnings']['time'],
			summ = WARNINGSIN + '{sample}_star/summary.txt',
			warn = WARNINGSIN,
			out  = OUTDIR + 'check_{sample}_star.warnings.txt'
		benchmark:
			WARNINGSOUT + '{sample}_star.warning.benchmark'
		run:
			f=open(input[0])
			for line in f:
				if "Number of input reads" in line:
					reads=line.split("\t")[1]
					if reads < 1000000:
						print("test")
"""
###
# FastQC Warning
# Reads out the Flags for overrepresented sequences, gc content and basewise quality.
# If one of these Flags is not PASS a warning file is produced.
###

rule fastqc_warning:
                input:
                        fastq = WARNINGSIN + '{sample}_fastqc.zip'
                output:
                        out = WARNINGSOUT + '{sample}_fastqc.warning.txt'
                params:
                        lsfoutfile = WARNINGSOUT + '{sample}_fastqc.warning.lsfout.log',
                        lsferrfile = WARNINGSOUT + '{sample}_fastqc.warning.lsferr.log',
                        scratch = config['tools']['warnings']['scratch'],
                        mem = config['tools']['warnings']['mem'],
                        time = config['tools']['warnings']['time'],
			summ = WARNINGSIN + '{sample}_fastqc/summary.txt',
			warn = WARNINGSIN,
			out  = OUTDIR + 'check_{sample}_fastqc.warnings.txt'
                benchmark:
                        WARNINGSOUT + '{sample}_fastqc.warning.benchmark'
                shell:
                        ("cd {params.warn} ; unzip -o {input.fastq} ; "+
			"over_seq=$(awk 'FNR == 10 {{print $1}}' {params.summ}) ;"+
			"gc=$(awk 'FNR == 6 {{print $1}}' {params.summ}) ;"+
			"seq_qual=$(awk 'FNR == 2 {{print $1}}' {params.summ}) ;" +
			" touch {output.out} ;" +
			"if [ $over_seq != 'PASS' ] ; then echo 'The module overrepresented sequences has a '$over_seq' flag' >> {params.out}; touch {params.out} ;fi;"+
			"if [ $gc != 'PASS' ] ; then echo 'The module gc content has a '$gc' flag' >> {params.out} ; fi;"+
			"if [ $seq_qual != 'PASS' ] ; then echo 'The module basewise quality has a '$seq_qual' flag' >> {params.out} ; fi;" )

rule haplotypeCaller:
        input:
                bam = HAPLOTYPECALLERIN + '{sample}.bam',
                bai = HAPLOTYPECALLERIN + '{sample}.bai',
                reference = {config['resources'][ORGANISM]['reference']},
                dbsnp = {config['resources'][ORGANISM]['dbSNP']},
		intervals = {config['resources'][ORGANISM]['intervals']}
        output:
                bam = HAPLOTYPECALLEROUT + '{sample}.g.vcf'
        params:
                lsfoutfile = HAPLOTYPECALLEROUT + '{sample}.lsfout.log',
                lsferrfile = HAPLOTYPECALLEROUT + '{sample}.lsferr.log',
                scratch = config['tools']['GATK']['haplotypeCaller']['scratch'],
                mem = config['tools']['GATK']['haplotypeCaller']['mem'],
                time = config['tools']['GATK']['haplotypeCaller']['time']
        benchmark:
                HAPLOTYPECALLEROUT + '{sample}.benchmark'
        threads:
                int(config['tools']['GATK']['haplotypeCaller']['threads'])
        shell:
                ('{config[tools][GATK][call]} ' +
                '-T HaplotypeCaller ' +
                '-R {input.reference} ' +
                '-I {input.bam} ' +
		'--genotyping_mode DISCOVERY ' +
		'--emitRefConfidence GVCF ' +
                '--dbsnp {input.dbsnp} ' +
		'-L {input.intervals} ' +
		'-o {output.bam} ')


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
            if exp not in expMap.keys():
                expMap[exp] = []
            expMap[exp].append(sample)
    if wildcards.experiment not in expMap.keys():
        return "UnknownExperiment"
    return expMap[wildcards.experiment]

def get_gVCFsFromExperimentId(wildcards):
    return expand('{sample}.g.vcf', sample = getSamplesFromExperimentId(wildcards))


def get_gVCFstoGenotypeFromExperimentId(wildcards):
    return expand(GENOTYPERIN + '{gvcf}', gvcf = get_gVCFsFromExperimentId(wildcards))

def prepend_V_GVCFs(wildcards):
    gVcfs = get_gVCFstoGenotypeFromExperimentId(wildcards)
    return ''.join(['--variant '+vcf+' ' for vcf in gVcfs])

rule joint_genotype:
        input:
                gvcf = get_gVCFstoGenotypeFromExperimentId,
        	reference = {config['resources'][ORGANISM]['reference']},
                dbsnp = {config['resources'][ORGANISM]['dbSNP']}
	output:
                vcf = GENOTYPEROUT + '{experiment}.vcf'
	params:
		variant = prepend_V_GVCFs,
		lsfoutfile = GENOTYPEROUT + '{experiment}.joint_genotype.lsfout.log',
		lsferrfile = GENOTYPEROUT + '{experiment}.joint_genotype.lsferr.log',
		scratch = config['tools']['GATK']['jointGenotype']['scratch'],
                mem = config['tools']['GATK']['jointGenotype']['mem'],
                time = config['tools']['GATK']['jointGenotype']['time']
	benchmark:
		GENOTYPEROUT + '{experiment}.joint_genotype.benchmark'
	threads:
		int(config['tools']['GATK']['jointGenotype']['threads'])
	shell:
		('{config[tools][GATK][call]} ' +
                '-T GenotypeGVCFs ' +
                '-R {input.reference} ' +
                '--dbsnp {input.dbsnp} ' +
                '{params.variant} ' +
                '-o {output.vcf} ')

rule snp_variantRecalibrator:
	input:
		vcf = GENOTYPEROUT + '{experiment}.vcf',
		reference = {config['resources'][ORGANISM]['reference']},
		hapmap = {config['resources'][ORGANISM]['hapmap']},
                omni = {config['resources'][ORGANISM]['omni']},
                G1000 = {config['resources'][ORGANISM]['1000G']},
                dbsnp = {config['resources'][ORGANISM]['dbSNP']}
	output:
		recal = GENOTYPEROUT + '{experiment}.snps.recal',
		tranches = GENOTYPEROUT + '{experiment}.snps.tranches'
	params:
                lsfoutfile = GENOTYPEROUT + '{experiment}.snps.varrecal.lsfout.log',
                lsferrfile = GENOTYPEROUT + '{experiment}.snps.varrecal.lsferr.log',
                scratch = config['tools']['GATK']['varRecalibrator']['scratch'],
                mem = config['tools']['GATK']['varRecalibrator']['mem'],
                time = config['tools']['GATK']['varRecalibrator']['time']
	benchmark:
		GENOTYPEROUT + '{experiment}.snps.varrecal.benchmark'
	threads:
		int(config['tools']['GATK']['varRecalibrator']['threads'])		
	shell:
		('{config[tools][GATK][call]} ' +
                '-T VariantRecalibrator ' +
                '-R {input.reference} ' +
                '-input {input.vcf} ' +
		'-recalFile {output.recal} ' +
		'-tranchesFile {output.tranches} ' +
		'-nt {threads} ' +
		'-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} ' +
		'-resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni} ' +
		'-resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.G1000} ' +
		'-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} ' +
		'-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP ' +
		'-mode SNP '
		)

rule indel_variantRecalibrator:
        input:
                vcf = GENOTYPEROUT + '{experiment}.snps.recal.vcf',
		reference = {config['resources'][ORGANISM]['reference']},
                mills = {config['resources'][ORGANISM]['Mills_indels']},
                dbsnp = {config['resources'][ORGANISM]['dbSNP']}
        output:
                recal = GENOTYPEROUT + '{experiment}.indels.recal',
                tranches = GENOTYPEROUT + '{experiment}.indels.tranches'
        params:
                lsfoutfile = GENOTYPEROUT + '{experiment}.indels.varrecal.lsfout.log',
                lsferrfile = GENOTYPEROUT + '{experiment}.indels.varrecal.lsferr.log',
                scratch = config['tools']['GATK']['varRecalibrator']['scratch'],
                mem = config['tools']['GATK']['varRecalibrator']['mem'],
                time = config['tools']['GATK']['varRecalibrator']['time']
        benchmark:
                GENOTYPEROUT + '{experiment}.indels.varrecal.benchmark'
        threads:
                int(config['tools']['GATK']['varRecalibrator']['threads'])
        shell:
                ('{config[tools][GATK][call]} ' +
                '-T VariantRecalibrator ' +
                '-R {input.reference} ' +
                '-input {input.vcf} ' +
                '-recalFile {output.recal} ' +
                '-tranchesFile {output.tranches} ' +
                '-nt {threads} ' +
		'--maxGaussians 4 ' +
                '-resource:mills,known=false,training=true,truth=true,prior=12.0 {input.mills} ' +
                '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} ' +
                '-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum ' +
                '-mode INDEL '
                )

rule snp_applyRecalibration:
	input:
                vcf = GENOTYPEROUT + '{experiment}.vcf',
                recal = GENOTYPEROUT + '{experiment}.snps.recal',
                tranches = GENOTYPEROUT + '{experiment}.snps.tranches',
		reference = {config['resources'][ORGANISM]['reference']}
	output:
                recal_vcf = GENOTYPEROUT + '{experiment}.snps.recal.vcf'
	params:
                lsfoutfile = GENOTYPEROUT + '{experiment}.snps.apply_recal.lsfout.log',
                lsferrfile = GENOTYPEROUT + '{experiment}.snps.apply_recal.lsferr.log',
                scratch = config['tools']['GATK']['applyRecalibrator']['scratch'],
                mem = config['tools']['GATK']['applyRecalibrator']['mem'],
                time = config['tools']['GATK']['applyRecalibrator']['time']
	benchmark:
                GENOTYPEROUT + '{experiment}.snps.apply_recal.benchmark'
	threads:
                int(config['tools']['GATK']['applyRecalibrator']['threads'])
	shell:
                ('{config[tools][GATK][call]} ' +
                '-T ApplyRecalibration ' +
                '-R {input.reference} ' +
                '-input {input.vcf} ' +
                '-recalFile {input.recal} ' +
                '-tranchesFile {input.tranches} ' +
		'-o {output.recal_vcf} ' +
		'--ts_filter_level 99.5 ' +
                '-mode SNP '
                )

rule indel_applyRecalibration:
        input:
                vcf = GENOTYPEROUT + '{experiment}.snps.recal.vcf',
                recal = GENOTYPEROUT + '{experiment}.indels.recal',
                tranches = GENOTYPEROUT + '{experiment}.indels.tranches',
                reference = {config['resources'][ORGANISM]['reference']}
        output:
                recal_vcf = GENOTYPEROUT + '{experiment}.recal.vcf'
        params:
                lsfoutfile = GENOTYPEROUT + '{experiment}.indels.apply_recal.lsfout.log',
                lsferrfile = GENOTYPEROUT + '{experiment}.indels.apply_recal.lsferr.log',
                scratch = config['tools']['GATK']['applyRecalibrator']['scratch'],
                mem = config['tools']['GATK']['applyRecalibrator']['mem'],
                time = config['tools']['GATK']['applyRecalibrator']['time']
        benchmark:
                GENOTYPEROUT + '{experiment}.indels.apply_recal.benchmark'
        threads:
                int(config['tools']['GATK']['applyRecalibrator']['threads'])
        shell:
                ('{config[tools][GATK][call]} ' +
                '-T ApplyRecalibration ' +
                '-R {input.reference} ' +
                '-input {input.vcf} ' +
                '-recalFile {input.recal} ' +
                '-tranchesFile {input.tranches} ' +
                '-o {output.recal_vcf} ' +
                '--ts_filter_level 99.0 ' +
                '-mode INDEL '
                )
	

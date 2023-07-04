# rule to annotate a vcf file using snpEff
rule snpEffAnnotation:
    input:
        vcf = '{sample}.vcf',
        snpEffDB  = {config['resources'][ORGANISM]['pathSnpEffDB']}
    output:
        vcf = '{sample}.snpEff.vcf',
        stats = '{sample}.snpEff.stats.html'
    params:
        lsfoutfile = '{sample}.snpEff.lsfout.log',
        lsferrfile = '{sample}.snpEff.lsferr.log',
        scratch = config['tools']['snpEff']['scratch'],
        mem = config['tools']['snpEff']['mem'],
        time = config['tools']['snpEff']['time'],
        dbName = config['tools']['snpEff']['dbName']
    threads:
        int(config['tools']['snpEff']['threads'])
    benchmark:
        '{sample}.snpEff.benchmark'
    shell: 
        '{config[tools][snpEff][call]} ann {params.dbName} -dataDir {input.snpEffDB} -nodownload -onlyProtein -canon -s {output.stats} {input.vcf} > {output.vcf}'
        
# This rule annotates a vcf file using snpSift and the dbSNP database
rule snpSift_dbSNP_Annotation:
    input:
        vcf = '{sample}.vcf',
        dbSnpDB  = {config['resources'][ORGANISM]['dbSNP']}
    output:
        vcf = '{sample}.dbSNP.vcf'
    params:
        lsfoutfile = '{sample}.dbSNP.vcf.lsfout.log',
        lsferrfile = '{sample}.dbSNP.vcf.lsferr.log',
        scratch = config['tools']['snpSift']['scratch'],
        mem = config['tools']['snpSift']['mem'],
        time = config['tools']['snpSift']['time']
    threads:
        int(config['tools']['snpSift']['threads'])
    benchmark:
        '{sample}.dbSNP.vcf.benchmark'
    shell:
        '{config[tools][snpSift][call]} annotate -noLog -noDownload {input.dbSnpDB} {input.vcf} > {output.vcf}'

# This rule annotates a vcf file using snpSift and the clinVar database
rule snpSift_clinVar_Annotation:
    input:
        vcf = '{sample}.vcf',
        clinVarDB  = {config['resources'][ORGANISM]['clinvar']}
    output:
        vcf = '{sample}.clinVar.vcf'
    params:
        lsfoutfile = '{sample}.clinVar.vcf.lsfout.log',
        lsferrfile = '{sample}.clinVar.vcf.lsferr.log',
        scratch = config['tools']['snpSift']['scratch'],
        mem = config['tools']['snpSift']['mem'],
        time = config['tools']['snpSift']['time']
    threads:
        int(config['tools']['snpSift']['threads'])
    benchmark:
        '{sample}.clinVar.vcf.benchmark'
    shell:
        '{config[tools][snpSift][call]} annotate -noLog -noDownload {input.clinVarDB} {input.vcf} > {output.vcf}'
        
# This rule annotates a vcf file using snpSift and the cosmic database
rule snpSift_COSMIC_Annotation:
    input:
        vcf = '{sample}.vcf',
        cosmicDB  = {config['resources'][ORGANISM]['cosmic']}
    output:
        vcf = '{sample}.cosmic.vcf'
    params:
        lsfoutfile = '{sample}.cosmic.vcf.lsfout.log',
        lsferrfile = '{sample}.cosmic.vcf.lsferr.log',
        scratch = config['tools']['snpSift']['scratch'],
        mem = config['tools']['snpSift']['mem'],
        time = config['tools']['snpSift']['time']
    threads:
        int(config['tools']['snpSift']['threads'])
    benchmark:
        '{sample}.cosmic.vcf.benchmark'
    shell:
        '{config[tools][snpSift][call]} annotate -noLog -noDownload {input.cosmicDB} {input.vcf} > {output.vcf}'
        
# This rule annotates a vcf file using snpSift and the dbNSFP database (functional annotation)
rule snpSift_dbNSFP_Annotation:
    input:
        vcf = '{sample}.vcf',
        dbNSFPDB  = {config['resources'][ORGANISM]['dbnsfp']}
    output:
        vcf = '{sample}.dbnsfp.vcf'
    params:
        lsfoutfile = '{sample}.dbnsfp.vcf.lsfout.log',
        lsferrfile = '{sample}.dbnsfp.vcf.lsferr.log',
        scratch = config['tools']['snpSift']['scratch'],
        mem = config['tools']['snpSift']['mem'],
        time = config['tools']['snpSift']['time']
    threads:
        int(config['tools']['snpSift']['threads'])
    benchmark:
        '{sample}.dbnsfp.vcf.benchmark'
    shell:
        'dbPath=$(readlink {input.dbNSFPDB}) ; {config[tools][snpSift][call]} dbnsfp -db $dbPath {input.vcf} > {output.vcf}'

# This rule calls a python script that converts a variant table with results from a panel sequencing to a vcf file
rule convertPanelFile_toVCF:
    input:
        panelFile = CONVERTPANELVCFIN + '{sample}.txt'
    output:
        vcf = CONVERTPANELVCFOUT + '{sample}.converted.vcf',
    params:
        lsfoutfile = CONVERTPANELVCFOUT + '{sample}.converted.lsfout.log',
        lsferrfile = CONVERTPANELVCFOUT + '{sample}.converted.lsferr.log',
        scratch = config['tools']['convertPanelVariantFile']['scratch'],
        mem = config['tools']['convertPanelVariantFile']['mem'],
        time = config['tools']['convertPanelVariantFile']['time'],
        panelType = config['tools']['convertPanelVariantFile']['panelType']
    threads:
        int(config['tools']['convertPanelVariantFile']['threads'])
    benchmark:
        CONVERTPANELVCFOUT + '{sample}.converted.benchmark'
    shell:
        '{config[tools][convertPanelVariantFile][call]} -i {input.panelFile} -o {output.vcf} -t {params.panelType}'
		
# This rule calls a python script parsing an annotated vcf file in order to create output files necessary for postprocessing, including database queries
rule parseAnnotatedVCF_forPostprocessingSteps_panelAdapted:
    input:
        vcf = PARSEANNOTATEDVCFIN + '{sample}.vcf',
        pathwayDB = {config['resources'][ORGANISM]['pathwayDB']}
    output:
        outList_all = DATABASEQUERY + '{sample}.overview.txt',
        outList_genes = DATABASEQUERY + '{sample}.overview.txt_distinctGeneNames.txt',
        outList_cBioportal = DATABASEQUERY + '{sample}.overview.txt_listForCbioportalQuery.txt'
    params:
        lsfoutfile = DATABASEQUERY + '{sample}.overview.lsfout.log',
        lsferrfile = DATABASEQUERY + '{sample}.overview.lsferr.log',
        scratch = config['tools']['parseAnnotatedGenesFromVCF']['scratch'],
        mem = config['tools']['parseAnnotatedGenesFromVCF']['mem'],
        time = config['tools']['parseAnnotatedGenesFromVCF']['time']
    threads:
        int(config['tools']['parseAnnotatedGenesFromVCF']['threads'])
    benchmark:
        DATABASEQUERY + '{sample}.overview.benchmark'
    shell:
        '{config[tools][parseAnnotatedGenesFromVCF][call]} {input.vcf} {output.outList_all} {input.pathwayDB}'
		
# query identified variants at dgidb
rule dgidbQuery:
    input:
        infile = DATABASEQUERY + '{sample}.txt_distinctGeneNames.txt'
    output:
        outfile = DATABASEQUERY + '{sample}.dgidb.txt',
        outfileCompleteTable = DATABASEQUERY + '{sample}.dgidb.txt.CompleteTable.txt',
        outfileGeneCategory = DATABASEQUERY + '{sample}.dgidb.txt.GeneCategories.txt'
    params:
        lsfoutfile = DATABASEQUERY + '{sample}.dgidbQuery.lsfout.log',
        lsferrfile = DATABASEQUERY + '{sample}.dgidbQuery.lsferr.log',
        scratch = config['tools']['queryDGIDB']['scratch'],
        mem = config['tools']['queryDGIDB']['mem'],
        time = config['tools']['queryDGIDB']['time'],
        minDatabaseNum = config['tools']['queryDGIDB']['minDatabaseNum']
    threads:
        int(config['tools']['queryDGIDB']['threads'])
    benchmark:
        DATABASEQUERY + '{sample}.dgidbQuery.benchmark'
    shell:
         '{config[tools][queryDGIDB][call]} {input.infile} {output.outfile} {params.minDatabaseNum}'
		
# clnical trials query
rule clinicalTrialsQuery:
    input:
        infile = '{sample}.dgidb.txt.CompleteTable.txt',
        downloadSuccess = DOWNLOADCLINICALTRIALSOUT + 'downloadSuccess.txt',
	#clinicalTrialsFolder = {config['resources']['general']['clinicalTrialsFolder']}
    output:
        outfile = '{sample}.dgidb.txt.CompleteTable.ClinicalTrials.txt'
    params:
        lsfoutfile = '{sample}.clinicalTrialsQuery.lsfout.log',
        lsferrfile = '{sample}.clinicalTrialsQuery.lsferr.log',
        scratch = config['tools']['queryClinicalTrials']['scratch'],
        mem = config['tools']['queryClinicalTrials']['mem'],
        time = config['tools']['queryClinicalTrials']['time'],
        cancerType = config['tools']['downloadClinicalTrials']['cancerType'],
        whiteList = config['tools']['queryClinicalTrials']['whiteList'],
        blackList = config['tools']['queryClinicalTrials']['blackList'],
        outDirec = DOWNLOADCLINICALTRIALSOUT
    threads:
        int(config['tools']['queryClinicalTrials']['threads'])
    benchmark:
        '{sample}.clinicalTrialsQuery.benchmark'
    shell:
        '{config[tools][queryClinicalTrials][call]} {input.infile} {output.outfile} {params.outDirec}/{params.cancerType}_clinicalTrials/ "{params.whiteList}" "{params.blackList}"'
        #'{config[tools][queryClinicalTrials][call]} {input.infile} {output.outfile} {input.clinicalTrialsFolder}/'
		
# download the clinical trials necessary for the query
rule downloadClinicalTrials:
    output:
        outfile = DOWNLOADCLINICALTRIALSOUT + 'downloadSuccess.txt'
    params:
        lsfoutfile = DOWNLOADCLINICALTRIALSOUT + 'downloadClinicalTrials.lsfout.log',
        lsferrfile = DOWNLOADCLINICALTRIALSOUT + 'downloadClinicalTrials.lsferr.log',
        scratch = config['tools']['downloadClinicalTrials']['scratch'],
        mem = config['tools']['downloadClinicalTrials']['mem'],
        time = config['tools']['downloadClinicalTrials']['time'],
        cancerType = config['tools']['downloadClinicalTrials']['cancerType'],
        outDirec = DOWNLOADCLINICALTRIALSOUT
    threads:
        int(config['tools']['downloadClinicalTrials']['threads'])
    benchmark:
        DOWNLOADCLINICALTRIALSOUT + 'downloadClinicalTrials.benchmark'
    shell:
        ('wget "https://clinicaltrials.gov/search?term={params.cancerType}&studyxml=true" ' + 
	'-O {params.outDirec}/{params.cancerType}_clinicalTrials.zip ; ' + 
	'unzip {params.outDirec}/{params.cancerType}_clinicalTrials.zip -d {params.outDirec}/{params.cancerType}_clinicalTrials ; ' +
	'touch {params.outDirec}/downloadSuccess.txt')
	
# Combine all database queries for the snv analysis
rule combineDatabaseQueries_snvs:
    input:
        inOverview = '{sample}.txt',
        inDGIDB = '{sample}.dgidb.txt.GeneCategories.txt',
        inClinicalTrials = '{sample}.dgidb.txt.CompleteTable.ClinicalTrials.txt',
        inCbio_cancerType = '{sample}.cbioportalSNV.txt_cancerTypeSpecific.txt',
        inCbio_maximum = '{sample}.cbioportalSNV.txt_maxPerVariant.txt'
    output:
        outTable = '{sample}.databaseQueriesSNV.txt',
        outTable_dgdidbIndependent = '{sample}.databaseQueriesSNV.txt_dgidbIndependent.txt'
    params:
        lsfoutfile = '{sample}.combineDatabaseQueries_snvs.lsfout.log',
        lsferrfile = '{sample}.combineDatabaseQueries_snvs.lsferr.log',
        scratch = config['tools']['combineDatabaseQueries_snvs']['scratch'],
        mem = config['tools']['combineDatabaseQueries_snvs']['mem'],
        time = config['tools']['combineDatabaseQueries_snvs']['time']
    threads:
        int(config['tools']['combineDatabaseQueries_snvs']['threads'])
    benchmark:
        '{sample}.combineDatabaseQueries_snvs.benchmark'
    shell:
        '{config[tools][combineDatabaseQueries_snvs][call]} {input.inOverview} {input.inDGIDB} {input.inCbio_maximum} {input.inCbio_cancerType} {input.inClinicalTrials} {output.outTable}'
		
# Combine all database queries for the cnv analysis
rule combineDatabaseQueries_cnvs:
    input:
        inOverview = '{sample}.txt',
        inDGIDB = '{sample}.dgidb.txt.GeneCategories.txt',
        inClinicalTrials = '{sample}.dgidb.txt.CompleteTable.ClinicalTrials.txt',
        pathwayDB = {config['resources'][ORGANISM]['pathwayDB']},
        inCbio_cancerType = '{sample}.cbioportalCNV.txt_cancerTypeSpecific.txt',
        inCbio_maximum = '{sample}.cbioportalCNV.txt_maxPerVariant.txt'
    output:
        outTable = '{sample}.databaseQueriesCNV.txt',
    params:
        lsfoutfile = '{sample}.combineDatabaseQueries_cnvs.lsfout.log',
        lsferrfile = '{sample}.combineDatabaseQueries_cnvs.lsferr.log',
        scratch = config['tools']['combineDatabaseQueries_cnvs']['scratch'],
        mem = config['tools']['combineDatabaseQueries_cnvs']['mem'],
        time = config['tools']['combineDatabaseQueries_cnvs']['time']
    threads:
        int(config['tools']['combineDatabaseQueries_cnvs']['threads'])
    benchmark:
        '{sample}.combineDatabaseQueries_cnvs.benchmark'
    shell:
        '{config[tools][combineDatabaseQueries_cnvs][call]} {input.inOverview} {input.inDGIDB} {input.inCbio_maximum} {input.inCbio_cancerType} {input.inClinicalTrials} {input.pathwayDB} {output.outTable}'

# query identified variants at cBioportal, SNV
# defined as localrule to avoid crash when several samples
localrules: cbioportalQuery_SNV
rule cbioportalQuery_SNV:
    input:
        infile = '{sample}.txt_listForCbioportalQuery.txt'
    output:
        oufileMaxStudy = '{sample}.cbioportalSNV.txt_maxPerVariant.txt',
	oufileCancerStudy = '{sample}.cbioportalSNV.txt_cancerTypeSpecific.txt',
        outfileAll = '{sample}.cbioportalSNV.txt'
    params:
        lsfoutfile = '{sample}.cbioportalSNV.lsfout.log',
        lsferrfile = '{sample}.cbioportalSNV.lsferr.log',
        scratch = config['tools']['queryCbioportal_snv']['scratch'],
        mem = config['tools']['queryCbioportal_snv']['mem'],
        time = config['tools']['queryCbioportal_snv']['time'],
        cancerType = config['tools']['queryCbioportal_snv']['cancerType']
    threads:
        int(config['tools']['queryCbioportal_snv']['threads'])
    benchmark:
        '{sample}.cbioportalSNV.benchmark'
    shell:
         '{config[tools][queryCbioportal_snv][call]} {input.infile} \"{params.cancerType}\" {output.outfileAll}'


# query identified variants at cBioportal, CNV
# defined as localrule to avoid crash when several samples
localrules: cbioportalQuery_CNV
rule cbioportalQuery_CNV:
    input:
        infile = '{sample}.txt_listForCbioportalQuery.txt'
    output:
        oufileMaxStudy = '{sample}.cbioportalCNV.txt_maxPerVariant.txt',
        oufileCancerStudy = '{sample}.cbioportalCNV.txt_cancerTypeSpecific.txt',
        outfileAll = '{sample}.cbioportalCNV.txt'
    params:
        lsfoutfile = '{sample}.cbioportalCNV.lsfout.log',
        lsferrfile = '{sample}.cbioportalCNV.lsferr.log',
        scratch = config['tools']['queryCbioportal_cnv']['scratch'],
        mem = config['tools']['queryCbioportal_cnv']['mem'],
        time = config['tools']['queryCbioportal_cnv']['time'],
        cancerType = config['tools']['queryCbioportal_cnv']['cancerType']
    threads:
        int(config['tools']['queryCbioportal_cnv']['threads'])
    benchmark:
        '{sample}.cbioportalCNV.benchmark'
    shell:
         '{config[tools][queryCbioportal_cnv][call]} {input.infile} \"{params.cancerType}\" {output.outfileAll}'


# extract all non-synonymous mutations and damaging mutations from overview table generated with combine_DatabaseQueries_snvs
rule extractProteinCodingMutations_fromSNVOverviewTable:
    input:
        inOverview = '{sample}.databaseQueriesSNV.txt_dgidbIndependent.txt'
    output:
        outTable = '{sample}.databaseQueriesSNV.txt_dgidbIndependent.txt_damaging.txt'
    params:
        lsfoutfile = '{sample}.extractProteinCodingMutations.lsfout.log',
        lsferrfile = '{sample}.extractProteinCodingMutations.lsferr.log',
        scratch = config['tools']['extractProteinCodingMutations']['scratch'],
        mem = config['tools']['extractProteinCodingMutations']['mem'],
        time = config['tools']['extractProteinCodingMutations']['time']
    threads:
        int(config['tools']['extractProteinCodingMutations']['threads'])
    benchmark:
        '{sample}.extractProteinCodingMutations.benchmark'
    shell:
        '{config[tools][extractProteinCodingMutations][call]} {input.inOverview}'

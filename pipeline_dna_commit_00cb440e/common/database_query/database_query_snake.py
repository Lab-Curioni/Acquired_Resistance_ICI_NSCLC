if not 'DATABASEQUERY' in globals():
    DATABASEQUERY = OUTDIR + 'databaseQuery/'

if not 'DGIDB_IN' in globals():
    DGIDB_IN = DATABASEQUERY
if not 'DGIDB_OUT' in globals():
    DGIDB_OUT = OUTDIR + 'dgidb/'

# query identified variants at dgidb
rule dgidbQuery:
    input:
        infile = DGIDB_IN + '{sample}.txt_distinctGeneNames.txt'
        #infile = DGIDB_IN + '{sample}.overview.txt'
    output:
        outfile = DGIDB_OUT + '{sample}.dgidb.txt',
	outfileCompleteTable = DGIDB_OUT + '{sample}.dgidb.txt.CompleteTable.txt',
        outfileGeneCategory = DGIDB_OUT + '{sample}.dgidb.txt.GeneCategories.txt'
    params:
        lsfoutfile = DGIDB_OUT + '{sample}.dgidbQuery.lsfout.log',
        lsferrfile = DGIDB_OUT + '{sample}.dgidbQuery.lsferr.log',
        scratch = config['tools']['queryDGIDB']['scratch'],
        mem = config['tools']['queryDGIDB']['mem'],
        time = config['tools']['queryDGIDB']['time'],
        minDatabaseNum = config['tools']['queryDGIDB']['minDatabaseNum'],
        colName_genes = config['tools']['queryDGIDB']['colName_genes']
    threads:
        config['tools']['queryDGIDB']['threads']
    benchmark:
        DGIDB_OUT + '{sample}.dgidbQuery.benchmark'
    shell:
         '{config[tools][queryDGIDB][call]} {input.infile} {output.outfile} {params.minDatabaseNum} {params.colName_genes}'

if not 'SL_DIR' in globals():
    SL_DIR = OUTDIR + 'spezialitaetenliste/'
# query the spezialitaetenliste
# ToDo: needs to be revised in order to respect the full Spezialitaetenliste
rule spezialitaetenlisteQuery:
    input:
        infile = DGIDB_OUT + '{sample}.dgidb.txt'
    output:
        outfile = SL_DIR + '{sample}.dgidb.txt.Spezialitaetenliste.txt'
    params:
        out = SL_DIR + 'Preparations.xml',
        outDir_dq = DATABASEQUERY,
        outputDir = SL_DIR,
        lsfoutfile = SL_DIR + 'Spezialitaetenliste.lsfout.log',
        lsferrfile = SL_DIR + 'Spezialitaetenliste.lsferr.log',
        scratch = config['tools']['querySpezialitaetenliste']['scratch'],
        mem = config['tools']['querySpezialitaetenliste']['mem'],
        time = config['tools']['querySpezialitaetenliste']['time']
    threads:
         config['tools']['querySpezialitaetenliste']['threads']
    benchmark:
        SL_DIR + 'spezialitaetenliste.benchmark'
    shell:
        """
        cd {params.outputDir}
        wget {config[tools][querySpezialitaetenliste][database_url]} -O XMLPublications.zip || true
        unzip -o XMLPublications.zip || true
        cd {params.outDir_dq}
        {config[tools][querySpezialitaetenliste][call]} {input.infile} {params.out} {output.outfile}
        """
        
if not 'CLINICALTRIALS_IN' in globals():
    CLINICALTRIALS_IN = DGIDB_OUT		
if not 'CLINICALTRIALS_OUT' in globals():
    CLINICALTRIALS_OUT = OUTDIR + 'clinicalTrials/'		
# clnical trials query
rule clinicalTrialsQuery:
    input:
        infile = CLINICALTRIALS_IN + '{sample}.dgidb.txt.CompleteTable.txt',
        downloadSuccess = CLINICALTRIALS_OUT + 'downloadSuccess.txt',
        #clinicalTrialsFolder = {config['resources']['general']['clinicalTrialsFolder']}
    output:
        outfile = CLINICALTRIALS_OUT + '{sample}.dgidb.txt.CompleteTable.ClinicalTrials.txt'
    params:
        lsfoutfile = CLINICALTRIALS_OUT + '{sample}.clinicalTrialsQuery.lsfout.log',
        lsferrfile = CLINICALTRIALS_OUT + '{sample}.clinicalTrialsQuery.lsferr.log',
        scratch = config['tools']['queryClinicalTrials']['scratch'],
        mem = config['tools']['queryClinicalTrials']['mem'],
        time = config['tools']['queryClinicalTrials']['time'],
        cancerType = config['tools']['downloadClinicalTrials']['cancerType'],
        whiteList = config['tools']['queryClinicalTrials']['whiteList'],
        blackList = config['tools']['queryClinicalTrials']['blackList'],
        outDirec = CLINICALTRIALS_OUT
    threads:
        config['tools']['queryClinicalTrials']['threads']
    benchmark:
        CLINICALTRIALS_OUT + '{sample}.clinicalTrialsQuery.benchmark'
    shell:
        '{config[tools][queryClinicalTrials][call]} {input.infile} {output.outfile} {params.outDirec}/{params.cancerType}_clinicalTrials/ "{params.whiteList}" "{params.blackList}"'
        #'{config[tools][queryClinicalTrials][call]} {input.infile} {output.outfile} {input.clinicalTrialsFolder}/'

# download the clinical trials necessary for the query
rule downloadClinicalTrials:
    output:
        outfile = CLINICALTRIALS_OUT + 'downloadSuccess.txt'
    params:
        lsfoutfile = CLINICALTRIALS_OUT + 'downloadClinicalTrials.lsfout.log',
        lsferrfile = CLINICALTRIALS_OUT + 'downloadClinicalTrials.lsferr.log',
        scratch = config['tools']['downloadClinicalTrials']['scratch'],
        mem = config['tools']['downloadClinicalTrials']['mem'],
        time = config['tools']['downloadClinicalTrials']['time'],
        cancerType = config['tools']['downloadClinicalTrials']['cancerType'],
        outDirec = CLINICALTRIALS_OUT
    threads:
        config['tools']['downloadClinicalTrials']['threads']
    benchmark:
        CLINICALTRIALS_OUT + 'downloadClinicalTrials.benchmark'
    shell:
        ('wget "https://clinicaltrials.gov/search?term={params.cancerType}&studyxml=true" ' +
	'-O {params.outDirec}/{params.cancerType}_clinicalTrials.zip ; ' +
	'unzip -o {params.outDirec}/{params.cancerType}_clinicalTrials.zip -d {params.outDirec}/{params.cancerType}_clinicalTrials ; ' +
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
        config['tools']['combineDatabaseQueries_snvs']['threads']
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
        config['tools']['combineDatabaseQueries_cnvs']['threads']
    benchmark:
        '{sample}.combineDatabaseQueries_cnvs.benchmark'
    shell:
        '{config[tools][combineDatabaseQueries_cnvs][call]} {input.inOverview} {input.inDGIDB} {input.inCbio_maximum} {input.inCbio_cancerType} {input.inClinicalTrials} {input.pathwayDB} {output.outTable}'

# query identified variants at cBioportal, SNV
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
        config['tools']['queryCbioportal_snv']['threads']
    benchmark:
        '{sample}.cbioportalSNV.benchmark'
    shell:
         '{config[tools][queryCbioportal_snv][call]} {input.infile} \"{params.cancerType}\" {output.outfileAll}'


# query identified variants at cBioportal, CNV
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
        config['tools']['queryCbioportal_cnv']['threads']
    benchmark:
        '{sample}.cbioportalCNV.benchmark'
    shell:
         '{config[tools][queryCbioportal_cnv][call]} {input.infile} \"{params.cancerType}\" {output.outfileAll}'

# extract all non-synonymous mutations and damaging mutations from overview table generated with combine_DatabaseQueries_snvs_woCbio
rule extractProteinCodingMutations_fromSNVOverviewTable_woCbio:
    input:
        inOverview = '{sample}.databaseQueriesSNV_woCbio.txt'
    output:
        outTable = '{sample}.databaseQueriesSNV_woCbio.txt_dgidbIndependent.txt_damaging.txt'
    params:
        lsfoutfile = '{sample}.extractProteinCodingMutations.lsfout.log',
        lsferrfile = '{sample}.extractProteinCodingMutations.lsferr.log',
        scratch = config['tools']['extractProteinCodingMutations']['scratch'],
        mem = config['tools']['extractProteinCodingMutations']['mem'],
        time = config['tools']['extractProteinCodingMutations']['time']
    threads:
        config['tools']['extractProteinCodingMutations']['threads']
    benchmark:
        '{sample}.extractProteinCodingMutations.benchmark'
    shell:
        '{config[tools][extractProteinCodingMutations][call]} {input.inOverview}'

if not 'ANNOTATECLINICAL_IN' in globals():
    ANNOTATECLINICAL_IN = PARSEANNOTATEDVCFOUT	
if not 'ANNOTATECLINICAL_OUT' in globals():
    ANNOTATECLINICAL_OUT = OUTDIR + 'clinicalAnnotation/'		

# Combine different database queries, annotate input table with clinical information
rule annotate_DE_clinicalInformation:
    input:
        infile = ANNOTATECLINICAL_IN + '{sample}.txt',
        #pathwayDB = config['resources'][ORGANISM]['pathwayDB'],
        inDGIDB = DGIDB_OUT + '{sample}.dgidb.txt.GeneCategories.txt',
        inClinicalTrials = CLINICALTRIALS_OUT + '{sample}.dgidb.txt.CompleteTable.ClinicalTrials.txt'
    output:
        outTable = ANNOTATECLINICAL_OUT + '{sample}.clinicalAnnotation.txt',
        outTable_dgdidbIndependent = ANNOTATECLINICAL_OUT + '{sample}.clinicalAnnotation.txt_dgidbIndependent.txt'
    params:
        lsfoutfile = ANNOTATECLINICAL_OUT + '{sample}.annotateClinical.lsfout.log',
        lsferrfile = ANNOTATECLINICAL_OUT + '{sample}.annotateClinical.lsferr.log',
        scratch = config['tools']['annotateClinical']['scratch'],
        mem = config['tools']['annotateClinical']['mem'],
        time = config['tools']['annotateClinical']['time'],
        variousParams = config['tools']['annotateClinical']['variousParams']
    threads:
        config['tools']['annotateClinical']['threads']
    benchmark:
        ANNOTATECLINICAL_OUT + '{sample}.annotateClinical.benchmark'
    shell:
        '{config[tools][annotateClinical][call]} --inputTable {input.infile} --outFile {output.outTable} --dgidb_categ {input.inDGIDB} --clinTrials {input.inClinicalTrials} {params.variousParams}'

if not 'MUTATIONALBURDEN_OUT' in globals():
    MUTATIONALBURDEN_OUT = OUTDIR + 'variants/mutationalBurden/'

# calculate the length of the regions bed file, once for each regions file
rule calculateRegionsLength:
    input:
        regions = config['resources'][ORGANISM]['regionsMutationalBurden'] 
    output:
        outFile= MUTATIONALBURDEN_OUT + config['tools']['calculateRegionsLength']['regionsName'] + '.regionsFile_length.txt'
    params:
        lsfoutfile = MUTATIONALBURDEN_OUT + config['tools']['calculateRegionsLength']['regionsName'] + '.regionsFile_length.lsfout.log',
        lsferrfile = MUTATIONALBURDEN_OUT + config['tools']['calculateRegionsLength']['regionsName'] + '.regionsFile_length.lsferr.log',
        scratch = config['tools']['calculateRegionsLength']['scratch'],
        mem = config['tools']['calculateRegionsLength']['mem'],
        time = config['tools']['calculateRegionsLength']['time']
    threads:
        config['tools']['calculateRegionsLength']['threads']
    benchmark:
        MUTATIONALBURDEN_OUT + config['tools']['calculateRegionsLength']['regionsName'] + '.regionsFile_length.benchmark'
    shell:
        '{config[tools][calculateRegionsLength][call]} {input.regions} {output.outFile}'


# extract all non-synonymous mutations and damaging mutations from overview table generated with combine_DatabaseQueries_snvs
# or annotate_DE_clinicalInformation
rule extractProteinCodingMutations_fromSNVOverviewTable:
    input:
        inOverview = ANNOTATECLINICAL_OUT + '{sample}.{placeholder}.txt',
        regionsLength = MUTATIONALBURDEN_OUT + config['tools']['calculateRegionsLength']['regionsName'] + '.regionsFile_length.txt'
    output:
        outTable = MUTATIONALBURDEN_OUT + '{sample}.{placeholder}.damaging.txt',
        mutLoad = MUTATIONALBURDEN_OUT + '{sample}.{placeholder}.mutationalBurden.txt'
    params:
        lsfoutfile = MUTATIONALBURDEN_OUT + '{sample}.{placeholder}.extractProteinCodingMutations.lsfout.log',
        lsferrfile = MUTATIONALBURDEN_OUT + '{sample}.{placeholder}.extractProteinCodingMutations.lsferr.log',
        scratch = config['tools']['extractProteinCodingMutations']['scratch'],
        mem = config['tools']['extractProteinCodingMutations']['mem'],
        time = config['tools']['extractProteinCodingMutations']['time'],
        outfileNameTag = MUTATIONALBURDEN_OUT + '{sample}.{placeholder}'
    threads:
        config['tools']['extractProteinCodingMutations']['threads']
    benchmark:
        MUTATIONALBURDEN_OUT + '{sample}.{placeholder}.extractProteinCodingMutations.benchmark'
    shell:
        '{config[tools][extractProteinCodingMutations][call]} {input.inOverview} {input.regionsLength} {params.outfileNameTag} {output.mutLoad}'



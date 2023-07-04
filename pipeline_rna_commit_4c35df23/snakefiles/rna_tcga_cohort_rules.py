if not 'RSEMCALCEXPRIN' in globals():
    RSEMCALCEXPRIN = STARALIGNOUT
if not 'RSEMCALCEXPROUT' in globals():
    RSEMCALCEXPROUT = OUTDIR + 'rsem/'

# This rule performs transcript quantification using rsem-calculate-expression
rule rsem_calculate_expression:
    input:
        bam = RSEMCALCEXPRIN + '{sample}_Aligned.toTranscriptome.out.bam'
    output:
        iso = RSEMCALCEXPROUT + '{sample}.isoforms.results',
        gene = RSEMCALCEXPROUT + '{sample}.genes.results'
    params:
        mem = config['tools']['RSEM_calc_expr']['mem'],
        time = config['tools']['RSEM_calc_expr']['time'],
        scratch = config['tools']['RSEM_calc_expr']['scratch'],
        lsferrfile = RSEMCALCEXPROUT + "log/{sample}.rsem_calc_expr.lsferr.log",
        lsfoutfile = RSEMCALCEXPROUT + "log/{sample}.rsem_calc_expr.lsfout.log",
        reference = config['tools']['RSEM_calc_expr']['reference'],
        prefix = RSEMCALCEXPROUT + '{sample}',
        variousParams = config['tools']['RSEM_calc_expr']['variousParams']
    benchmark:
        RSEMCALCEXPROUT + 'log/{sample}.rsem_calc_expr.benchmark'
    threads:
        config['tools']['RSEM_calc_expr']['threads']
    shell:
        '{config[tools][RSEM_calc_expr][call]} {params.variousParams} -p {threads} --bam {input.bam} {params.reference} {params.prefix}'


if not 'STRIPTABSIN' in globals():
    STRIPTABSIN = RSEMCALCEXPROUT
if not 'STRIPTABSOUT' in globals():
    STRIPTABSOUT = OUTDIR + 'rsem/'

# This rule strips trailing tabs from the file {sample}.isoforms.results that is an output file of RSEM
rule strip_trailing_tabs:
    input:
        infile = STRIPTABSIN + '{sample}.isoforms.results'
    output:
        out = STRIPTABSOUT + '{sample}.isoforms.results.bak'
    params:
        mem = config['tools']['strip_trailing_tabs']['mem'],
        time = config['tools']['strip_trailing_tabs']['time'],
        scratch = config['tools']['strip_trailing_tabs']['scratch'],
        lsferrfile = STRIPTABSOUT + "log/{sample}.striptabs.lsferr.log",
        lsfoutfile = STRIPTABSOUT + "log/{sample}.striptabs.lsfout.log",
    benchmark:
        STRIPTABSOUT + 'log/{sample}.striptabs.benchmark'
    threads:
        config['tools']['strip_trailing_tabs']['threads']
    shell:
        '{config[tools][strip_trailing_tabs][call]} --input {input.infile} --temp {output.out}'


if not 'PRUNEISOIN' in globals():
    PRUNEISOIN = STRIPTABSOUT
if not 'PRUNEISOOUT' in globals():
    PRUNEISOOUT = OUTDIR + 'rsem/'

# This rule prunes the isoforms from the gene quant file which is of type {sample}.genes.results
# and is given out by rsem_calculate_expression
rule prune_isoforms:
    input:
        infile = PRUNEISOIN + '{sample}.genes.results'
    output:
        out = PRUNEISOOUT + '{sample}.genes.results.pruned.txt'
    params:
        mem = config['tools']['prune_isoforms']['mem'],
        time = config['tools']['prune_isoforms']['time'],
        scratch = config['tools']['prune_isoforms']['scratch'],
        lsferrfile = PRUNEISOOUT + "log/{sample}.prune_isoforms.lsferr.log",
        lsfoutfile = PRUNEISOOUT + "log/{sample}.prune_isoforms.lsfout.log"
    benchmark:
        PRUNEISOOUT + 'log/{sample}.prune_isoforms.benchmark'
    threads:
        config['tools']['prune_isoforms']['threads']
    shell:
        'sed /^uc0/d {input.infile} > {output.out}'

if not 'NORMQUANTIN' in globals():
    NORMQUANTIN = STRIPTABSOUT
if not 'NORMQUANTOUT' in globals():
    NORMQUANTOUT = OUTDIR + 'rsem/'

# This rule normalizes the gene expression quantification data of files of type {sample}.genes.results that are given out by rsem_calculate_expression
rule normalize_gene_quant:
    input:
        ingene = NORMQUANTIN + '{sample}.genes.results.pruned.txt'
    output:
        outgene = NORMQUANTOUT + '{sample}.genes.results.pruned.normalized.txt'
    params:
        mem = config['tools']['normalize_gene_quant']['mem'],
        time = config['tools']['normalize_gene_quant']['time'],
        scratch = config['tools']['normalize_gene_quant']['scratch'],
        lsferrfile = NORMQUANTOUT + "log/{sample}.normalize_gene_quant.lsferr.log",
        lsfoutfile = NORMQUANTOUT + "log/{sample}.normalize_gene_quant.lsfout.log",
        variousParams = config['tools']['normalize_gene_quant']['variousParams']
    benchmark:
        NORMQUANTOUT + 'log/{sample}.normalize_gene_quant.benchmark'
    threads:
        config['tools']['normalize_gene_quant']['threads']
    shell:
        '{config[tools][normalize_gene_quant][call]} {params.variousParams} -t 1000 -o {output.outgene} {input.ingene} '

# This rule normalizes the gene expression quantification data of files of type {sample}.isoforms.results that are given out by rsem_calculate_expression and 
rule normalize_iso_quant:
    input:
        iniso = NORMQUANTIN + '{sample}.isoforms.results',
        bakiso = NORMQUANTIN + '{sample}.isoforms.results.bak'
    output:
        outiso = NORMQUANTOUT + '{sample}.isoforms.results.normalized.txt'
    params:
        mem = config['tools']['normalize_gene_quant']['mem'],
        time = config['tools']['normalize_gene_quant']['time'],
        scratch = config['tools']['normalize_gene_quant']['scratch'],
        lsferrfile = NORMQUANTOUT + "log/{sample}.normalize_iso_quant.lsferr.log",
        lsfoutfile = NORMQUANTOUT + "log/{sample}.normalize_iso_quant.lsfout.log",
        variousParams = config['tools']['normalize_gene_quant']['variousParams']
    benchmark:
        NORMQUANTOUT + 'log/{sample}.normalize_iso_quant.benchmark'
    threads:
        config['tools']['normalize_gene_quant']['threads']
    shell:
        '{config[tools][normalize_gene_quant][call]} {params.variousParams} -t 300 -o {output.outiso} {input.iniso}'


if not 'CHANGEHEADIN' in globals():
    CHANGEHEADIN = NORMQUANTOUT
if not 'CHANGEHEADOUT' in globals():
    CHANGEHEADOUT = OUTDIR + 'rsem/'

# This rule changes the header of the fil
rule change_header:
    input:
        infile = CHANGEHEADIN + '{sample}.genes.results.pruned.normalized.txt'
    output:
        out = CHANGEHEADOUT + '{sample}.genes.results.pruned.normalized.header.txt'
    params:
        mem = config['tools']['change_header']['mem'],
        time = config['tools']['change_header']['time'],
        scratch = config['tools']['change_header']['scratch'],
        lsferrfile = CHANGEHEADOUT + "log/{sample}.change_header.lsferr.log",
        lsfoutfile = CHANGEHEADOUT + "log/{sample}.change_header.lsfout.log",
        sampleName = '{sample}'
    benchmark:
        CHANGEHEADOUT + 'log/{sample}.change_header.benchmark'
    threads:
        config['tools']['change_header']['threads']
    shell:
        'sed "1 s/0.0000/{params.sampleName}/" {input.infile} > {output.out}'


if not 'PARSECOHORTIN' in globals():
    PARSECOHORTIN = CHANGEHEADOUT
if not 'PARSECOHORTOUT' in globals():
    PARSECOHORTOUT = OUTDIR + 'output_tcga/'

# This rule parses the tcga pipeline results and cohort values to a useful table format
rule parse_cohort_comparison:
    input:
        infile = PARSECOHORTIN + '{sample}.genes.results.pruned.normalized.header.txt'
    output:
        out = PARSECOHORTOUT + '{sample}.final_rna_table.txt'
    params:
        mem = config['tools']['parse_cohort_comparison']['mem'],
        time = config['tools']['parse_cohort_comparison']['time'],
        scratch = config['tools']['parse_cohort_comparison']['scratch'],
        lsferrfile = PARSECOHORTOUT + "log/{sample}.parse_cohort_comparison.lsferr.log",
        lsfoutfile = PARSECOHORTOUT + "log/{sample}.parse_cohort_comparison.lsfout.log",
        cohort = config['tools']['parse_cohort_comparison']['cohort']
    benchmark:
        PARSECOHORTOUT + 'log/{sample}.parse_cohort_comparison.benchmark'
    threads:
        config['tools']['parse_cohort_comparison']['threads']
    shell:
        '{config[tools][parse_cohort_comparison][call]} {params.cohort} {input.infile} {output.out}'


if not 'CLEANTABLEIN' in globals():
    CLEANTABLEIN = PARSECOHORTOUT
if not 'CLEANTABLEOUT' in globals():
    CLEANTABLEOUT = OUTDIR = "cleaned/"

# This rule removes lines in the file PARSECOHORTOUT/{sample}.final_rna_table.txt that contain "?" for a gene name
rule remove_questionmark:
    input:
        infile = CLEANTABLEIN + '{sample}.final_rna_table.txt'
    output:
        out = CLEANTABLEOUT + '{sample}.final_rna_table.txt'
    params:
        mem = config['tools']['remove_questionmark']['mem'],
        time = config['tools']['remove_questionmark']['time'],
        scratch = config['tools']['remove_questionmark']['scratch'],
        lsferrfile = CLEANTABLEOUT + "log/{sample}.remove_questionmark.lsferr.log",
        lsfoutfile = CLEANTABLEOUT + "log/{sample}.remove_questionmark.lsfout.log"
    benchmark:
        CLEANTABLEOUT + 'log/{sample}.remove_questionmark.benchmark'
    threads:
        config['tools']['remove_questionmark']['threads']
    shell:
        'grep -vP "^\?\\t" {input.infile} > {output.out}'


if not 'FILTERGENESIN' in globals():
    FILTERGENESIN = ADDDGIDBOUT
if not 'FILTERGENESOUT' in globals():
    FILTERGENESOUT = ADDDGIDBOUT

# This rule filters the output table of the rule parse_cohort_comparison of type {sample}.final_rna_table.txt. It returns all rows of {sample}.final_rna_table.txt of genes that are on a list of clinically interesting genes
rule filter_genesOfInterest:
    input:
        finalTable = FILTERGENESIN + '{sample}.final_rna_table_dgidb_results.txt'
    output:
        out = FILTERGENESOUT + '{sample}.final_rna_table_dgidb_results_filtered_clinical_genes.txt'
    params:
        mem = config['tools']['filter_genesOfInterest']['mem'],
        time = config['tools']['filter_genesOfInterest']['time'],
        scratch = config['tools']['filter_genesOfInterest']['scratch'],
        lsferrfile = FILTERGENESOUT + "/log/{sample}.filter_genesOfInterest.lsferr.log",
        lsfoutfile = FILTERGENESOUT + "/log/{sample}.filter_genesOfInterest.lsfout.log",
        genelist = config['tools']['filter_genesOfInterest']['genelist']
    benchmark:
        FILTERGENESOUT + '/log/{sample}.filter_genesOfInterest.benchmark'
    threads:
        config['tools']['filter_genesOfInterest']['threads']
    shell:
        '{config[tools][filter_genesOfInterest][call]} --inFile {input.finalTable} --outFile {output.out} --geneList {params.genelist}'



if not 'BOXPLOTIN' in globals():
    BOXPLOTIN = CHANGEHEADOUT
if not 'BOXPLOTOUT' in globals():
    BOXPLOTOUT = OUTDIR + 'plots/'

#GENELIST = FILTERGENESOUT + '{sample}.final_rna_table_filtered_clinical_genes.txt'

#def getGenesOfInterest(genelist):
#    output = []
#    if not 'GENELIST' in globals():
#        return ['NOGENELIST']
#    try:
#        open(genelist,'r')
#    except IOError:
#        return ['NOGENELIST']
#    with open(genelist, "r") as infile:
#        headerLine = infile.readline()
#        for line in infile:
#            lineSplit = line.strip().split("\t")
#            geneID = lineSplit[0]
#            output.append(geneID)
#        infile.close()
#    return output

# This rule generates the variable GENELIST containing a list with all genes of interest




# This rule plots a boxplot of the expression values of all cohort samples and marks the position of the patient sample in the boxplot. It does so for the subset of genes in the list that it is given.
rule boxplot_expression_value:
    input:
        patient = expand(BOXPLOTIN + '{sample}.genes.results.pruned.normalized.header.txt', sample = getSampleNames()),
        filtered = expand(FILTERGENESOUT + '{sample}.final_rna_table_dgidb_results_filtered_clinical_genes.txt', sample = getSampleNames())
    output:
        out = BOXPLOTOUT + '{sample}_boxplots_success.txt'
    params:
        mem = config['tools']['boxplot_expression_value']['mem'],
        time = config['tools']['boxplot_expression_value']['time'],
        scratch = config['tools']['boxplot_expression_value']['scratch'],
        lsferrfile = BOXPLOTOUT + "log/{sample}.boxplot_expression_value.lsferr.log",
        lsfoutfile = BOXPLOTOUT + "log/{sample}.boxplot_expression_value.lsfout.log",
        cohort = config['tools']['boxplot_expression_value']['cohort'],
        sampleName = '{sample}',
        outdir = BOXPLOTOUT
    benchmark:
        BOXPLOTOUT + 'log/{sample}.boxplot_expression_value.benchmark'
    threads:
        config['tools']['boxplot_expression_value']['threads']
    shell:
        'Rscript {config[tools][boxplot_expression_value][call]} {params.outdir} {params.cohort} {input.patient} {input.filtered} {params.sampleName}; touch {params.outdir}{params.sampleName}_boxplots_success.txt'


if not 'DGIDBIN' in globals():
    DGIDBIN = OUTDIR + 'databaseQuery/'
if not 'DGIDBOUT' in globals():
    DGIDBOUT = OUTDIR + 'databaseQuery/'

# query identified variants at dgidb
rule dgidbQuery:
    input:
        infile = DGIDBIN + '{sample}.final_rna_table.txt'
    output:
        outfile = DGIDBOUT + '{sample}.dgidb.txt',
        outfileCompleteTable = DGIDBOUT + '{sample}.dgidb.txt.CompleteTable.txt',
        outfileGeneCategory = DGIDBOUT + '{sample}.dgidb.txt.GeneCategories.txt'
    params:
        lsfoutfile = DGIDBOUT + '{sample}.dgidbQuery.lsfout.log',
        lsferrfile = DGIDBOUT + '{sample}.dgidbQuery.lsferr.log',
        scratch = config['tools']['queryDGIDB']['scratch'],
        mem = config['tools']['queryDGIDB']['mem'],
        time = config['tools']['queryDGIDB']['time'],
        minDatabaseNum = config['tools']['queryDGIDB']['minDatabaseNum'],
        colName_genes = config['tools']['queryDGIDB']['colName_genes']
    threads:
        config['tools']['queryDGIDB']['threads']
    benchmark:
        DGIDBOUT + '{sample}.dgidbQuery.benchmark'
    shell:
         '{config[tools][queryDGIDB][call]} {input.infile} {output.outfile} {params.minDatabaseNum} {params.colName_genes}'


if not 'ADDDGIDBIN' in globals():
    ADDDGIDBIN = CLEANTABLEOUT
if not 'ADDDGIDBOUT' in globals():
    ADDDGIDBOUT = OUTDIR + 'final_table/'

# This rule adds the result of the rDGIdb query to the final output table of the tcga workflow
rule add_dgidb_column:
    input:
        inTable = ADDDGIDBIN + '{sample}.final_rna_table.txt',
        dgidbTable = DGIDBOUT + '{sample}.dgidb.txt'
    output:
        out = ADDDGIDBOUT + '{sample}.final_rna_table_dgidb_results.txt'
    params:
        mem = config['tools']['add_dgidb_column']['mem'],
        time = config['tools']['add_dgidb_column']['time'],
        scratch = config['tools']['add_dgidb_column']['scratch'],
        lsferrfile = ADDDGIDBOUT + "/log/{sample}.add_dgidb_column.lsferr.log",
        lsfoutfile = ADDDGIDBOUT + "/log/{sample}.add_dgidb_column.lsfout.log",
        hgncCol = config['tools']['add_dgidb_column']['hgncCol']
    benchmark:
        ADDDGIDBOUT + '/log/{sample}.add_dgidb_column.benchmark'
    threads:
        config['tools']['add_dgidb_column']['threads']
    shell:
        '{config[tools][add_dgidb_column][call]} --inFile {input.inTable} --outFile {output.out} --hgncCol {params.hgncCol} --dgidbTable {input.dgidbTable}'

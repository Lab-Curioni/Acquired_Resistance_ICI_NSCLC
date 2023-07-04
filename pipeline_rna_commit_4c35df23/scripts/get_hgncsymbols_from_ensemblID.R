####################################
## File name: get_hgncsymbols_from_ensemblID.R
## Anne Richter, May 2018
## R Version: 3.4.3 (2017-11-30)
## In command line type:
## Rscript get_hgncsymbols_from_ensemblID.R --inFile /path/to/inputTable --outFile /path/to/outputTable
####################################

library(optparse)
library(biomaRt)

# argument parser
option_list = list(
  make_option("--inFile", type = "character", help = "Input tablen that might be the output of htseq-count ({sample}.htseq_counts.txt).
              Tab separated. Without header. In the first column are the ensembl gene IDs with version, in the second are the gene counts."),
  make_option("--outFile", type = "character", help = "/path/to/outputPlot enter full absolute path")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# load test data
#inTable <- read.table("/Users/annricht/Documents/icgc_test/dev_filter_icgc_output_table/Patient_6_RNA.htseq_counts.txt", header =F)

# load files from command line
inTable <- read.table(opt$inFile, header = FALSE)
names(inTable) <- c("ensembl_ID", "gene_count")

# remove transcript ID of ensembl ID, e.g. ENSG00000000003.10 to ENSG00000000003
# inTable$ensembl_ID_short <- gsub("\\..*","", inTable$ensembl_ID)

# generate mart object
mart_obj = useEnsembl(biomart="ensembl", host="www.ensembl.org", dataset="hsapiens_gene_ensembl", GRCh=37)

#if (genomeVersion != 'hg19'){
#  print('Use genome version GRCh38.')
#  mart_obj = useEnsembl(biomart="ensembl", host="www.ensembl.org", dataset="hsapiens_gene_ensembl")
#}

# get table with different gene IDs from biomaRt
hgncGeneMapping = getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "hgnc_symbol", "description"),  mart = mart_obj, uniqueRows=T)
hgncGeneMapping_noNa = na.omit(hgncGeneMapping)

# table with column "ensembl_gene_id_version" unique
hgncGeneMapping_noNa_unique = hgncGeneMapping_noNa[!duplicated(hgncGeneMapping_noNa$ensembl_gene_id_version),]

# checking number of hgnc symbols returned by biomaRt
hgnc_symbols <- hgncGeneMapping$hgnc_symbol
hgnc_symbols_uni <- unique(hgnc_symbols)
print("Number of hgnc symbols retrieved by BiomaRt:")
print(length(hgnc_symbols))
print("Number of unique hgnc symbols retrieved by BiomaRt:")
print(length(hgnc_symbols_uni))

# subset of inTable that is also in list by mart
inTable_subset <- subset(inTable, inTable$ensembl_ID %in% hgncGeneMapping_noNa_unique$ensembl_gene_id_version)

# merge of inTable and table with different IDs returned by biomaRt
mergedTable <- merge(inTable_subset, hgncGeneMapping_noNa_unique, by.x = "ensembl_ID", by.y = "ensembl_gene_id_version")

print("Number of unique hgnc symbols assigned to ensembl IDs in inTable:")
length(unique(mergedTable$hgnc_symbol))

# change order of columns for further processing, hgnc_symbol is now first column
mergedTable_order <- mergedTable[c(4,1,2,5)]

write.table(mergedTable_order, file = opt$outFile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

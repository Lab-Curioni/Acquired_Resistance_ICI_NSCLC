library(gplots)
library(DESeq2)
library("AnnotationDbi")
library("org.Hs.eg.db")


# Assuming that project_dir is your project directory
project_dir = "/project/dir/"

output_dir = paste0(project_dir, 'out_rna/log2fc/')

gene_counts_dir = paste0(project_dir, 'out_rna/gene_counts/')
gene_counts_suffix = '_RNA_non_rRNA.htseq_counts.txt' 
samples = c('BB_before', 'BB_resistance', 'VA_before', 'VA_resistance')

# read in gene counts
K = read.table(paste0(gene_counts_dir,samples[1],gene_counts_suffix), sep="\t", header=F, quote="", row.names=1)
for (i in 2:length(samples)){
    D = read.table(paste0(gene_counts_dir,samples[i],gene_counts_suffix), sep="\t", header=F, quote="", row.names=1)
    K = cbind(K, D[,1])
}

colnames(K) = samples

# remove __
x = grep("^__", rownames(K))
K = K[-x,]

# remove everything after period in ENSG ids
K_rownames = rownames(K)
rownames(K) = sub("\\..*", "", K_rownames)


# add gene symbols
gene_symbols = mapIds(org.Hs.eg.db, keys=rownames(K), column="SYMBOL", keytype="ENSEMBL", multiVals="first")


inds <- which(!is.na(gene_symbols))
found_genes <- gene_symbols[inds]
 
# subset your data frame based on the found_genes
df <- K[names(found_genes), ]
df$SYMBOL = found_genes

K3 = aggregate(list(BB_before=df$BB_before, BB_resistance=df$BB_resistance, VA_before=df$VA_before, VA_resistance=df$VA_resistance), by = list(df$SYMBOL), FUN=sum)
M = as.matrix(K3[,c(2:5)])
rownames(M) = K3$Group.1


countdata = M
coldata = as.data.frame(c('before', 'resistance', 'resistance', 'before'))
rownames(coldata) = c('BB_before', 'BB_resistance', 'VA_before', 'VA_resistance')
colnames(coldata) = c('label')


BB = c(1,2)
VA = c(4,3)

ddsFullCountTable_BB <- DESeqDataSetFromMatrix(
countData = countdata[,BB],
colData = coldata[BB,, drop=FALSE],
design = ~ label)

dds_BB = DESeq(ddsFullCountTable_BB)
res_BB = results(dds_BB)
order_BB = order(res_BB[,2], na.last=NA, decreasing=TRUE)

write.table(res_BB[order_BB, 2, drop=FALSE], file=paste0(output_dir, 'BB_log2fc.txt'), quote=FALSE, sep="\t", col.names=NA)


ddsFullCountTable_VA <- DESeqDataSetFromMatrix(
countData = countdata[,VA],
colData = coldata[VA,, drop=FALSE],
design = ~ label)

dds_VA = DESeq(ddsFullCountTable_VA)
res_VA = results(dds_VA)
order_VA = order(res_VA[,2], na.last=NA, decreasing=TRUE)

write.table(res_VA[order_VA, 2, drop=FALSE], file=paste0(output_dir, 'VA_log2fc.txt'), quote=FALSE, sep="\t", col.names=NA)

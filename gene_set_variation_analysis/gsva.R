library(GSVA)
library(GSEABase)
library(gplots)
library(limma)
library("AnnotationDbi")
library("org.Hs.eg.db")


# Assuming that project_dir is your project directory
project_dir = "/project/dir/"

output_dir = paste0(project_dir, 'out_rna/gsva/')
output_prefix = 'myprefix_'
gene_sets_gmt = paste0(project_dir, 'data/h.all.v6.2.symbols.gmt') #hallmark
gene_sets=getGmt(gene_sets_gmt)

selected_gene_sets = c('HALLMARK_ANGIOGENESIS','HALLMARK_APOPTOSIS', 'HALLMARK_DNA_REPAIR', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'HALLMARK_GLYCOLYSIS', 'HALLMARK_IL2_STAT5_SIGNALING', 'HALLMARK_IL6_JAK_STAT3_SIGNALING', 'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_KRAS_SIGNALING_DN', 'HALLMARK_KRAS_SIGNALING_UP', 'HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB');

gene_counts_dir = '/cluster/project/nexus/curioni/hiltbrunner_resistancePilotPart2_2018/out_rna/gene_counts/'
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

X = gsva(M, gene_sets)$es.obs;

pdf(paste0(output_dir, output_prefix, "gsva_pathways.pdf"))
par(mar=c(7,4,4,2)+0.1)
heatmap.2(X[selected_gene_sets,c(1,2,4,3)], dendrogram="none", main="", scale="none", density.info="none", trace="none", cexRow=0.5, cexCol=0.5, margins=c(12,12), Colv=FALSE, col=bluered(100))
dev.off()

write.table(X[selected_gene_sets,], paste0(output_dir, output_prefix, "gsva_pathways.txt"), quote=F, sep='\t', col.names=NA)



subtype = data.frame(row.names=samples)
subtype[,1] = c('before', 'resistant', 'resistant', 'before')

adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)
design = model.matrix(~factor(subtype$V1))
colnames(design) = c("before", "beforeVsResistant")
fit = lmFit(X[selected_gene_sets,], design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef="beforeVsResistant", number=Inf)
write.table(allGeneSets[,c('logFC', 'P.Value', 'adj.P.Val')], paste0(output_dir, output_prefix, "gsva_pathways.DE.txt"), quote=F, sep='\t', col.names=NA)



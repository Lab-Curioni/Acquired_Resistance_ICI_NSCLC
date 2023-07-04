#source("https://bioconductor.org/biocLite.R")
#biocLite("GenVisR")
require(RColorBrewer)
library(ggplot2)
library(reshape2)
library(GenVisR)
args = commandArgs(trailingOnly = T);

## Read command line arguments.
input_file = args[1]
output_file = args[2]
w = as.integer(args[3])
h = as.integer(args[4])

Mutations_SNV <- read.table(input_file, header=T, stringsAsFactors = FALSE)
mutation_class <- c("frameshift_variant", "stop_gained", "missense_variant", "3_prime_UTR_variant", "synonymous_variant",
                      "intron_variant", "splice_acceptor_variant&intron_variant", "splice_region_variant&intron_variant",
                      "5_prime_UTR_premature_start_codon_gain_variant", "downstream_gene_variant", "upstream_gene_variant",
                      "intergenic_region", "protein_protein_contact", "no_mutation_found")
mutation_class_snvs <-c("frameshift_variant", "missense_variant", "start_lost", "stop_gained", "stop_lost")
Mutations_SNV <- Mutations_SNV[Mutations_SNV$variant_class %in% mutation_class_snvs,]
col_snvs<- c(brewer.pal(8, "Set1"),brewer.pal(12, "Set3"))

pdf(file=output_file, width=w, height=h)
waterfall(Mutations_SNV, fileType = "Custom", variant_class_order = mutation_class_snvs, mainGrid=TRUE, mainDropMut= TRUE, plotMutBurden = FALSE, mainXlabel = TRUE,  mainPalette = col_snvs, main_geneLabSize=8, mainLabelSize=4)
dev.off()

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

## Read in the mutation table.
Mutations <- read.table(input_file, header=T,stringsAsFactors = FALSE)

## Name Mutation Type in an extra columns in the cnv table.
Mutations$CNV <- ifelse(as.integer(Mutations$variant_class)>=6, "strong Amplification", "-")
Mutations$CNV <- ifelse(as.integer(Mutations$variant_class)<6&as.integer(Mutations$variant_class)>2, "Amplification", Mutations$CNV)
Mutations$CNV <- ifelse(as.integer(Mutations$variant_class)==1, "hemizygous Deletion", Mutations$CNV)
Mutations$CNV <- ifelse(as.integer(Mutations$variant_class)==0, "homozygous Deletion", Mutations$CNV)
names(Mutations) <- c("sample","gene","CNV", "variant_class")

mutation_class <-c("strong Amplification", "Amplification", "hemizygous Deletion", "homozygous Deletion")
Mutations<- Mutations[Mutations$variant_class %in% mutation_class,]
col <- c("red", "orange", "chartreuse3", "chartreuse4")

pdf(file=output_file, width=w, height=h)
waterfall(Mutations, fileType = "Custom", variant_class_order = mutation_class, mainGrid=TRUE, mainDropMut= TRUE, plotMutBurden = FALSE, mainXlabel = TRUE,  mainPalette = col)
dev.off()

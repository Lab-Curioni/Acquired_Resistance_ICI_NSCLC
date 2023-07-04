#!/opt/local/bin/RScript

# version adapted from script part of ngs-pipe, authored by Jochen Singer

#library(ape,lib="/nfs/nas22.ethz.ch/fs2201/eth_nexus_pht_2/utilities/Rlibs/")
#library(gplots,lib="/nfs/nas22.ethz.ch/fs2201/eth_nexus_pht_2/utilities/Rlibs/")

library(ape)
library(gplots)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

df <- read.table(args[1], header = TRUE)
print("1")
print(names(df))

sampleNames = unlist(lapply(names(df), function(x) strsplit(x,"[.]")[[1]][2]))
print("2")
names(df) = sampleNames
names(df)
gt_dist <- dist.gene(t(df), method = "pairwise")
m_gt_dist <- as.matrix(gt_dist)
dim(m_gt_dist)
png(args[2],width = 2200, height = 1600, res = 300)
#heatmap.2(m_gt_dist, trace="none", margins = c(10, 10), main = "Sample similarity", cexRow=0.45, cexCol=0.45, row.names = names(m_gt_dist), col.names = names(m_gt_dist))
heatmap.2(m_gt_dist, trace="none", margins = c(4, 8), main = "Sample similarity", cexRow=0.45, cexCol=0.45)
dev.off()

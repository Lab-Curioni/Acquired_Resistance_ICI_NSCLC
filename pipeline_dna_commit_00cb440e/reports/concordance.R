library(gplots)
library(gsubfn)
library(gplots)

args = commandArgs(trailingOnly = T);
names = strsplit(args[1], ",")
names = unlist(names)
results_folder = args[2]
output_file = paste(results_folder, "/concordance.pdf", sep="")

getNumberPart <- function(x) {
  pat <- "(-?(\\d*\\.*\\d+|\\d+\\.))"
  strapply(x, pattern=pat, FUN=as.numeric, simplify=TRUE, empty=NA)
}

N = length(names)

concordance = matrix(0L, nrow=N, ncol=N)
rownames(concordance) <- names
colnames(concordance) <- names
for (i in 1:N) {
  for (j in 1:N) {
    D = file(paste0(results_folder, "/", names[i],"_", names[j], "_concordance.txt"), "r")
    percentage <- round(getNumberPart(readLines(D, n=1)))
    print(paste(names[i], names[j]))
    close(D)
    concordance[i,j]=percentage
  }
}

pdf(file=output_file, width=10, height=10)
heatmap.2(concordance, cellnote=concordance, notecex=0.2, Rowv=F, Colv=F, scale="none", trace="none", margins=c(12,9), cexRow=0.3, cexCol=0.3, notecol="black")
dev.off()

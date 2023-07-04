args <- commandArgs(trailingOnly = TRUE)
base_dir=args[1]
fnames=Sys.glob(file.path(base_dir, "*.rseqc.geneBodyCoverage.txt"))
samples = gsub(".rseqc.geneBodyCoverage.txt","",basename(fnames))

N=length(samples)
x=c(1:100)
gbcov = list()
colors=rainbow(N)

#max_dens=0

for (i in seq_along(samples)){
  input = fnames[i]
  print(input)
  x <- read.table(input, row.names=1)
  #x[2,] = x[2,]/rowSums(x)[2]
  x[2,] = x[2,]/max(x[2,])
  gbcov[[i]] <- as.matrix(x)
  #max_dens = max(max_dens, max(x[2,]))
}



pdf(paste0(base_dir, "/genebodycoverage.pdf"))
#plot(gbcov[[1]][1,], gbcov[[1]][2,], xlim=c(0,100), type='l', col=colors[1], main="Gene body coverage", ylim=c(0,max_dens), xlab="5' - 3'", ylab="Density")
plot(gbcov[[1]][1,], gbcov[[1]][2,], xlim=c(0,100), type='l', col=colors[1], main="Gene body coverage", ylim=c(0,1), xlab="5' -> 3'", ylab="Coverage")
print(N)
for (i in 2:N){
  lines(gbcov[[i]][1,], gbcov[[i]][2,], col=colors[i])
    print(i)
    print(gbcov[[i]])
}

legend("bottomright",samples,cex=0.7, lty=c(1,1), lwd=c(2,2),col=colors)

dev.off()



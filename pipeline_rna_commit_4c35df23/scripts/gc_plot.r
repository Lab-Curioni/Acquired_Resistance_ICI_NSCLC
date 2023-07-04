args <- commandArgs(trailingOnly = TRUE)
base_dir=args[1]
fnames=Sys.glob(file.path(base_dir,"*rseqc.read_gc.GC.xls"))
samples = gsub(".rseqc.read_gc.GC.xls","",basename(fnames))

N=length(samples)
x=c(1:100)
gc = list()
colors=rainbow(N)

max_cnt=0
for (i in seq_along(samples)){
  input = fnames[i]
  x <- read.table(input, header=T)
  x1 <- x[order(x$GC.),]
  integral = sum(x1$read_count[-1]*diff(x1$GC.))
  x1$read_count <- x1$read_count/integral
  gc[[i]] <- x1
  max_cnt = max(max_cnt, max(x1$read_count))
}

pdf(paste0(base_dir, "/gc_content.pdf"))
plot(gc[[1]]$GC., gc[[1]]$read_count, xlim=c(0,100), type='l', col=colors[1], main="GC content", ylim=c(0,max_cnt), xlab="%GC content", ylab="Density")
if (N>=2){
  for (i in 2:N){
    lines(gc[[i]]$GC., gc[[i]]$read_count, col=colors[i])
  }
} else {
 lines(gc[[1]]$GC., gc[[1]]$read_count, col=colors[1])
}
legend("topright",samples,cex=0.7, lty=c(1,1), lwd=c(2,2),col=colors)

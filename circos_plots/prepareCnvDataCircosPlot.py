

import sys

incnvs=sys.argv[1]
out = sys.argv[2]


in_cnvs=open(incnvs,'r')
out_file = open(out, 'w')


for line in in_cnvs:
	if line.startswith("Chromosome"):
		out_file.write("Chromosmome\tStart\tEnd\tCopyNumber\tOption\n")
		continue
	split=line.split("\t")
	chrom=split[0].replace("chr","hs")
	start=split[1]
	end=split[2]
	copynum=split[5]
	print(copynum)
	if float(copynum) == 0:
		option="color=yellow"
	if float(copynum)>0 and float(copynum) <2:
		option="color=vlblue"
	if float(copynum)>2 and float(copynum)<=3.5:	
		option="color=lred"
	if float(copynum)>3.5:
		option="color=dred"
	out_file.write(chrom+"\t"+start+"\t"+end+"\t"+copynum+"\t"+option+"\n")






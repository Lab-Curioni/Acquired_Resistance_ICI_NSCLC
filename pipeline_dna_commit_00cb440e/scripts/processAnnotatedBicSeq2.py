#!/usr/bin/env python

'''
Process annotated bicSeq2 results
Franziska Singer, January 2016
'''

import sys
import numpy as np
import os
import re


if "-h" in sys.argv[1]:
	print "Process annotated bicSeq2 results."
	print "Usage: python processAnnotatedBicSeq2.py [infileBicSeq2] [outfile]"
	sys.exit(1)

infileBicSeq2 = sys.argv[1]
outfileName = sys.argv[2]

infile = open(infileBicSeq2,'r')
outfile = open(outfileName,'w')
outfileCbioportal = open(outfileName + "_listForCbioportalQuery.txt" ,'w')
outfileCbioportal.write("Gene\tCopynumber\n")

outfile.write("Chromosome\tstart\tend\tlog2.copyRatio\tpvalue\testimated copy number\tGenes\n")

currentGenes = []
currentCNV = []

allGenes = []

allCNVs = 1
allCNVGenes = 0
#lines are of the form: chr7	10013	159128640	0.408493	0.00661941	3.4	chr7	149596	154868	uc003sio.1	AK024243

for line in infile:
	lineSplit = line.strip().split("\t")
	chrom = lineSplit[0]
	start = lineSplit[1]
	end = lineSplit[2]
	
	if len(currentCNV) == 0:
		# first line, simply extract and initialize
		currentCNV = [chrom, start, end, lineSplit[3],lineSplit[4],lineSplit[5]]
		
	if (chrom != currentCNV[0]) or (start != currentCNV[1]) or (end != currentCNV[2]):
		# arrived at a new cnv
		allCNVs += 1
		outfile.write("\t".join(currentCNV) + "\t" + ";".join(currentGenes) + "\n")
		currentCNV = [chrom, start, end, lineSplit[3],lineSplit[4],lineSplit[5]]
		currentGenes = []
		
	# also include lines without any annotation
	
	if len(lineSplit) <= 8:
		# no genes annotated
		continue

	if lineSplit[10] not in currentGenes:
		currentGenes.append(lineSplit[10])
		allCNVGenes += 1
	if lineSplit[10] not in allGenes:
		allGenes.append(lineSplit[10])
		outfileCbioportal.write("%s\t%s\n"%(lineSplit[10],lineSplit[5]))

outfile.write("\t".join(currentCNV) + "\t" + ";".join(currentGenes) + "\n")  # last cnv

infile.close()
outfile.close()
outfileCbioportal.close()

outfileDGIDB = open(outfileName + "_distinctGeneNames.txt", 'w')
outfileDGIDB.write("gene_names\n")
outString = "\n".join(allGenes)
outfileDGIDB.write(outString)
outfileDGIDB.close()

print "Parsed %s CNVs overlapping %s genes." %(allCNVs,allCNVGenes)

#!/usr/bin/env python

'''
Process annotated excavator results
Franziska Singer, January 2016
'''

import sys
import numpy as np
import os
import re


if "-h" in sys.argv[1]:
	print "Process annotated excavator results."
	print "Usage: python processAndFilterAnnotatedExcavator.py [infileExcavator] [outfile] [probabilityThreshold]"
	sys.exit(1)

infileExcavator = sys.argv[1]
outfileName = sys.argv[2]
probThreshold = float(sys.argv[3])

infile = open(infileExcavator,'r')
outfile = open(outfileName,'w')
outfileCbioportal = open(outfileName + "_listForCbioportalQuery.txt" ,'w')
outfileCbioportal.write("Gene\tCopynumber\n")

outfile.write("Chromosome\tStart\tEnd\tSegment\tCNF\tCN\tCall\tProbCall\tGenes\n")

currentGenes = []
currentCNV = []

allGenes = []

allCNVs = 1
allCNVGenes = 0
#filteredCNVs = 0
#lines are of the form: chr1	1871754	3747839	0.666372717427075	3.17415530362767	3	1	0.950569590047249	chr1	1884751	1935276	uc001aim.1	KIAA1751
probIndex = 7
for line in infile:
	lineSplit = line.strip().split("\t")
	chrom = lineSplit[0]
	start = lineSplit[1]
	end = lineSplit[2]
	probCall = lineSplit[probIndex]
	
	if float(probCall) < probThreshold:
		#filteredCNVs += 1
		continue
		
	if len(currentCNV) == 0:
		# first line, simply extract and initialize
		currentCNV = [chrom, start, end, lineSplit[3],lineSplit[4],lineSplit[5],lineSplit[6],lineSplit[7]]
		
	if (chrom != currentCNV[0]) or (start != currentCNV[1]) or (end != currentCNV[2]):
		# arrived at a new cnv
		allCNVs += 1
		outfile.write("\t".join(currentCNV) + "\t" + ";".join(currentGenes) + "\n")
		currentCNV = [chrom, start, end, lineSplit[3],lineSplit[4],lineSplit[5],lineSplit[6],lineSplit[7]]
		currentGenes = []

	# also include lines without any annotation
	
	if len(lineSplit) <= 8:
		# no genes annotated
		continue
		
	if lineSplit[12] not in currentGenes:
		currentGenes.append(lineSplit[12])
		allCNVGenes += 1
	if lineSplit[12] not in allGenes:
		allGenes.append(lineSplit[12])
		outfileCbioportal.write("%s\t%s\n"%(lineSplit[12],lineSplit[5]))

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

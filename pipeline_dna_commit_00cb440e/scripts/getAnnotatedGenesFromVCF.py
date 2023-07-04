#!/usr/bin/env python

'''
Get annotated genes from a vcf file
Franziska Singer, August 2015
'''

import sys
import numpy as np
import os
import re


if len(sys.argv) <= 1:
	print "Get annotated genes from a vcf file."
	print "Usage: python getAnnotatedGenesFromVCF.py [vcfFile] [outfile]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Get annotated genes from a vcf file."
	print "Usage: python getAnnotatedGenesFromVCF.py [vcfFile] [outfile]"
	sys.exit(1)

vcfFile = sys.argv[1]
outfileName = sys.argv[2]

infile = open(vcfFile,'r')
outfile = open(outfileName,'w')
outfileSift = open(outfileName + "_listForSift.txt",'w')
outfile.write("Chromosome\tPosition\tChange\tFrequency\tGenes\tAnnotation\tVariantID\tFunctionalAnnotation\n")
geneArr = []
allVariants = 0

for line in infile:
	if line.startswith("#"):
		continue
	allVariants += 1
	lineSplit = line.strip().split("\t")
	if "ANN=" in line:
		annoTag = lineSplit[7].split("ANN=")[1].split(";")[0]	
		geneName = annoTag.split("|")[3]
		geneAnno = annoTag.split("|")[9] + "|" + annoTag.split("|")[10]
		geneNameTempArr = [geneName]
		geneAnnoTempArr = [geneAnno]
		annoTagArr = annoTag.split(",")
		for i in range(0,len(annoTagArr)):
			annoTemp = annoTagArr[i]
			geneNameTemp = annoTemp.split("|")[3]
			if geneNameTemp not in geneNameTempArr:
				geneNameTempArr.append(geneNameTemp)
				geneAnnoTemp = annoTag.split("|")[9] + "|" + annoTag.split("|")[10]
				geneAnnoTempArr.append(geneAnnoTemp)
			if geneNameTemp not in geneArr:
				geneArr.append(geneNameTemp)
		allTempGeneNames = ";".join(geneNameTempArr)
		allTempGeneAnnos = ";".join(geneAnnoTempArr)
		# note: currently the frequency extraction is adapted to the combination of varscan, mutect, and a variable other caller! TODO!
		frequency = lineSplit[10].split(":")[-2]
		if "%" in frequency:
			frequency = frequency.split("%")[0]
		else:
			frequency = float(frequency) * 100.0
			
		functionalAnnotate = "."
		if "dbNSFP" in line:  # extract functional annotation
			functionalAnnotate = "dbNSFP_CADD" + lineSplit[7].split(";set=")[0].split("dbNSFP_CADD")[1]  # add dbNSFP_CADD to use dbNSFP_CADD for the split but include the annotation
			
		outfile.write(lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[3] + ">" + lineSplit[4] + "\t" + str(frequency) + "\t" + allTempGeneNames + "\t" + allTempGeneAnnos + "\t" + lineSplit[2] + "\t" + functionalAnnotate + "\n")
		chromTag = re.sub("[^0-9]", "", lineSplit[0])
		if len(chromTag) == 0:
			chromTag = lineSplit[0][-1]
		outfileSift.write(chromTag + "," + lineSplit[1] + "," + lineSplit[3] + "," + lineSplit[4] + "\n")
				
	else:
		print "Error! no annotation tag in line %s." %(line)
		
		
infile.close()
outfile.close()
outfileSift.close()

outfile1 = open(outfileName + "_distinctGeneNames.txt",'w')
outString = "\n".join(geneArr)
outfile1.write(outString)
outfile1.close()
print "Extracted %s genes associated to %s variants." %(len(geneArr),allVariants)

#!/usr/bin/env python

'''
Combine the snv overview table with results cnv calling, create overview for appendix creation
snv input files are the results of the following script extractProteinCodingMutations.py
Franziska Singer, March 2016
'''

import sys
import numpy as np
import os
import re
from natsort import natsort

def cnvTypeFunction(cnv):
	if cnv < 2:
		if cnv <= 0.5:
			return(0,"homozygous deletion")
		else:
			return (1,"hemizygous deletion")
	if cnv == 2:
		return (2,"neutral")
	if cnv > 2:
		if cnv >= 3.5:
			return(round(cnv),"high-level amplification")
		else:
			return(3,"low-level amplification")

if len(sys.argv) <= 1:
	print "Combine snv and cnv and create appendix table."
	print "Usage: python combine_snv_cnv_forAppendix.py [variantInfofile] [cnvFile] [damagingSNVs_file] [outputName]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Combine snv and cnv and create appendix table."
	print "Usage: python combine_snv_cnv_forAppendix.py [variantInfofile] [cnvFile] [damagingSNVs_file] [outputName]"
	sys.exit(1)

	
variantFile = sys.argv[1]  # this file is the result of extractProteinCodingMutations.py  (use the dgidb independent file)
cnvFile = sys.argv[2]      # either wgs (bicseq2) oder exome (excavator) based  (use file sample.CNVannotated.overview.txt or sample.CNVannotated.overviewFiltered.txt, respectively)
variantfile_damaging = sys.argv[3] # this file is the result of extractProteinCodingMutations.py  (use the dgidb independent _damaging file)
outfileName = sys.argv[4]

outfileAll = open(outfileName+"_allMutations.txt", 'w')
outfileAll.write("Gene\tMutation\tVariant frequency(%%)\\newline or Copy number\tDatabase identifier\tTechnical confidence\n")
outfileDamage = open(outfileName+"_damagingMutations.txt", 'w')
outfileDamage.write("Gene\tMutation\tVariant frequency(%%)\tDatabase identifier\tTechnical confidence\n")
outfileCNV = open(outfileName+"_cnvsAppendix.txt", 'w')
outfileCNV.write("Chromosome\tMutation\tStart\tEnd\tEstimated copy number\tPvalue\n")

infileDamage = open(variantfile_damaging,'r')
infileDamage.readline()  # skip first line

damagingGenes = 0
for lineDamage in infileDamage:
	lineSplit = lineDamage.strip().split("\t")
	frequency = lineSplit[4]
	genes = lineSplit[0]
	mutationInfo = lineSplit[5]
	databaseID = lineSplit[7]
	
	if len(databaseID) == 1:
		databaseID = "-"
	else:
		if ";" in databaseID:
			dataSplit = databaseID.split(";")
			databaseID = "; ".join(dataSplit)
	
	genesSplit = genes.split(";")
	for i in range(0,len(genesSplit)):
		gene = genesSplit[i]
		if "_" in gene:
			splitTemp = gene.split("_")
			gene = "\\_".join(splitTemp)
		damagingGenes += 1
		mutInfo = mutationInfo.split(";")[i]
		if "|p." in mutInfo:
			mutInfo = mutInfo.replace("|","\\newline ")  # tex \\newline
		else:
			mutInfo = mutInfo.replace("|","")
		if "_" in mutInfo:
			splitTemp = mutInfo.split("_")
			mutInfo = "\\_".join(splitTemp)
		
		outfileDamage.write(gene + "\t" + mutInfo + "\t" + frequency + "\t" + databaseID + "\t1\n") #%(gene,mutInfo,frequency,databaseID)  # note: technical confidence for now always 1 in lack of a suitable calculation. Change! TODO!
	
print "Found %s damaging genes (this includes multiple gene annotations for one position)." %(damagingGenes)
outfileDamage.close()

mySortedChromNames = []
dict_myCNVs = {} # chromosome to list of cnvs

infileCNV = open(cnvFile,'r')
cnvNum = 0
cnvGenesNum = 0

line1 = infileCNV.readline()
isWGS = False

if "pvalue" in line1: # then WGS (bicSeq2)
	isWGS = True
	
for lineCNV in infileCNV:
	lineSplit = lineCNV.strip().split("\t")
	
	
	chrom = lineSplit[0]
	if chrom not in mySortedChromNames:
		mySortedChromNames.append(chrom)
		dict_myCNVs[chrom] = []
		
	genes = []
	if isWGS:
		if len(lineSplit) > 6: # in case no annotation found for this CNV
			genes = lineSplit[6]
	else:
		if len(lineSplit) > 8:
			genes = lineSplit[8]
	
	copyNum = float(lineSplit[5])  # same column for wgs and exome
		
	(copyNum,type) = cnvTypeFunction(copyNum)
	
	myTempCNV = [chrom,type,lineSplit[1],lineSplit[2],str(copyNum),lineSplit[4]]
	dict_myCNVs[chrom].append(myTempCNV)
	cnvNum += 1
	
	if len(genes) == 0: # no annotated genes found
		continue
	
	geneSplit = genes.split(";")
	allCopyNumGenes = []
	for gene in geneSplit:
		if gene not in allCopyNumGenes:
			if "_" in gene:
				splitTemp = gene.split("_")
				gene = "\\_".join(splitTemp)
		
			outfileAll.write("%s\t%s\t%s\t-\t1\n" %(gene,type,copyNum))  # note: technical confidence for now always 1 in lack of a suitable calculation. Change! TODO!
		cnvGenesNum += 1
	
print "Found %s cnvs, affecting %s genes." %(cnvNum,cnvGenesNum)
infileCNV.close()

mySortedChromNames = natsort(mySortedChromNames)
for chrom in mySortedChromNames:
	for cnv in dict_myCNVs[chrom]:
		outfileCNV.write("%s\n"%("\t".join(cnv)))
		
outfileCNV.close()


infile = open(variantFile,'r')

allVariants = 0
infile.readline()

for line in infile:
	lineSplit = line.strip().split("\t")
	frequency = lineSplit[4]
	genes = lineSplit[0]
	mutationInfo = lineSplit[5]
	databaseID = lineSplit[7]
	
	if len(databaseID) == 1:
		databaseID = "-"
	else:
		if ";" in databaseID:
			dataSplit = databaseID.split(";")
			databaseID = "; ".join(dataSplit)
	
	genesSplit = genes.split(";")
	for i in range(0,len(genesSplit)):
		gene = genesSplit[i]
		if "_" in gene:
			splitTemp = gene.split("_")
			gene = "\\_".join(splitTemp)
		allVariants += 1
		mutInfo = mutationInfo.split(";")[i]
		if "|p." in mutInfo:
			mutInfo = mutInfo.replace("|","\\newline ")  # tex
		else:
			mutInfo = mutInfo.replace("|","")
			
		outfileAll.write("%s\t%s\t%s\t%s\t1\n" %(gene,mutInfo,frequency,databaseID))  # note: technical confidence for now always 1 in lack of a suitable calculation. Change! TODO!

outfileAll.close()
infile.close()
print "Found %s variants (this includes multiple gene annotations for one position)." %(allVariants)

callSortDamage = "head -n 1 %s > %s ; grep -v -e \"Mutation\" %s | awk 'NF!=1' | sort >> %s" %(outfileName+"_damagingMutations.txt",outfileName+"_damagingMutations_sorted.txt",outfileName+"_damagingMutations.txt",outfileName+"_damagingMutations_sorted.txt")
callSortAll = "head -n 1 %s > %s ; grep -v -e \"Mutation\" %s | awk 'NF!=1' | sort >> %s" %(outfileName+"_allMutations.txt",outfileName+"_allMutations_sorted.txt",outfileName+"_allMutations.txt",outfileName+"_allMutations_sorted.txt")

os.system(callSortDamage)
os.system(callSortAll)

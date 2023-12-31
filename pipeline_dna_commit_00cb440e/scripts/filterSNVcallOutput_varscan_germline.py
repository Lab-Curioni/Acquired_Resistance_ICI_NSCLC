#!/usr/bin/env python

'''
Filter the output of SNV calling by VarScan (nonsomatic mdoe)
current filter criteria: SS=1, SS=5, indelError, str10, for LOH normal variant frequency is below threshold
as well as homopolymers and p-value threshold (adjusted p-value) and strand bias check
deprecated: 
freqCompareFilter = float(sys.argv[5])    # if set to 1, no filter is applied. Otherwise we try to filter out variants which frequencies are rather close in normal and tumor tissue. E.g. 50% and 60% -> with high coverage, p-value will be small; however, we might not want to regard this variant
threshVarFreqDiff = float(sys.argv[6])	  # if set to 0, no filtering is applied. Otherwise, remove all variants which frequencies only differ by at most threshVarFreqDiff percent -> NOTE: this will penalize low frequency mutations!
Franziska Singer, January 2016
'''

import sys
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

if len(sys.argv) <= 1:
	print "Filter the output of SNV calling by VarScan (germline mode)."
	print "Usage: python filterSNVcallOutput_varscan_germline_mochBreast.py [inputVCF] [outfile] [minimumVariantSupport] [pvalue] [minimumNucleotideCoverage] [checkStrands y/n] [tumorFreqThreshold] [filterHomoPolymer y/n] [filterSilent y/n] [filterCommon y/n]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Filter the output of SNV calling by VarScan (germline mode)."
	print "Usage: python filterSNVcallOutput_varscan_germline_mochBreast.py [inputVCF] [outfile] [minimumVariantSupport] [pvalue] [minimumNucleotideCoverage] [checkStrands y/n] [tumorFreqThreshold] [filterHomoPolymer y/n] [filterSilent y/n] [filterCommon y/n]"
	sys.exit(1)

########
#functions
########	
def checkBases(allAllels):
		baseArr = [0,0,0,0]  # number of A, C, T, G that is found in the alleles
		for thisAllel in allAllels:
			if ("A" in thisAllel) or ("a" in thisAllel):
				baseArr[0] += 1
			if ("C" in thisAllel) or ("c" in thisAllel):
				baseArr[1] += 1
			if ("T" in thisAllel) or ("t" in thisAllel):
				baseArr[2] += 1
			if ("G" in thisAllel) or ("g" in thisAllel):
				baseArr[3] += 1
		return baseArr

inputVCF = sys.argv[1]
outputVCF = sys.argv[2]
minimumVariantSupport = int(sys.argv[3])		# minimum number of reads supporting variant in tumor
pvalueThreshold = float(sys.argv[4]) 		   # threshold for adjusted pvalue
minimumNucleotideCoverage = int(sys.argv[5])				# minimum nucleotide coverage in tumor 
checkStrands = sys.argv[6]
tumorFreqThreshold = float(sys.argv[7])	  # if set to 0, no filtering is applied. Otherwise remove variants with tumor frequency less than threshold
filterHomopolymer = sys.argv[8]   		  # if set to y, we check whether ref and alt allel show an indel which only consists of equal bases, e.g. AA -> A ; for ion torrent, this is likely to be a sequencing error
filterSilent = sys.argv[9]				  # if set to y, all mutations annotated with "synonymous" or "silent" are filtered
removeCommon = sys.argv[10]				  # filter variants annotated as COMMON or with MAF>5%


pValueArr = []    # p-value list

infileTemp = open(inputVCF,'r')   # only for p-value adjustment
for lineTemp in infileTemp:
	if lineTemp.startswith("#"):
		continue
	lineSplitTemp = lineTemp.strip().split("\t")
	thisPValTemp = float(float(lineSplitTemp[9].split(":")[7]))
	pValueArr.append(thisPValTemp)
infileTemp.close()

pValueArr.sort()
pValueArr = pValueArr[::-1]   # descending order
numTests = len(pValueArr)
dictPvalueToRank = {}

rankP = 0 
currentP = 0.0
for thisP in pValueArr:   # assign rank
	if thisP != currentP:
		rankP += 1
		currentP = thisP
		dictPvalueToRank[thisP] = rankP
	else:
		rankP += 1

infile = open(inputVCF,'r')
outfile = open(outputVCF,'w')

# different filter option counts, for statistics
allVariants = 0
filteredVariants = 0
filteredVariantsVarScanDefault = 0
filteredVariantsPVal = 0
filteredVariantsCoverageNucMin = 0
filteredVariantsCoverageVarMin = 0
filteredVariantsHP = 0
filteredVariantsSyn = 0
filteredVariantsStrand = 0
filteredVariantsFreq = 0
filteredVariantsCommon = 0

countCommon = 0
countHeterozygous = 0
countHomozygous = 0
countIndel = 0
countSNP = 0

for line in infile:
	if line.startswith("#"):
		outfile.write(line)
		continue
		
	allVariants += 1
	lineSplit = line.strip().split("\t")
	
	filter = False
	
	#VarScan filters
	if "indelError" in line:
		filteredVariantsVarScanDefault += 1
		filter = True
	if "str10" in line:
		filteredVariantsVarScanDefault += 1
		filter = True
	
	# first pvalue
	thisPVal = float(lineSplit[9].split(":")[7])
	pVal_adjusted = (float(numTests)/float(dictPvalueToRank[thisPVal])) * thisPVal
	if pVal_adjusted > pvalueThreshold:
		filteredVariantsPVal += 1
		filter = True
	
	# homopolymers
	if "y" in filterHomopolymer:
		refAllel = lineSplit[3]
		altAllel = lineSplit[4]
		
		baseArr = checkBases([refAllel,altAllel])
		sumDiffBases = 0
		for baseTemp in baseArr:
			if baseTemp > 0:
				sumDiffBases += 1
		if sumDiffBases == 1:
			filteredVariantsHP += 1
			filter = True
	
	# filter synonomous
	if "y" in filterSilent:
		# check if annotated as silent or synonymous
		if ("synonymous" in line) or ("silent" in line):
			filteredVariantsSyn += 1
			filter = True
	
	# strand bias
	if "y" in checkStrands:
		# check if one strand shows zero support
		strandFo = int(lineSplit[9].split(":")[12])
		strandRe = int(lineSplit[9].split(":")[13])
		
		if (strandFo == 0) or (strandRe == 0):
			filteredVariantsStrand += 1
			filter = True
			
	# minimum variant frequency
	if (minimumVariantSupport > 0):
		# filter all sites that do not achieve the specified minimum coverage of variant supporting reads
		thisVariantCoverage = int(lineSplit[9].split(":")[5])
		if thisVariantCoverage < minimumVariantSupport:
			filteredVariantsCoverageVarMin += 1
			filter = True
	
	# minimum nucleotide coverage
	if minimumNucleotideCoverage > 0:
		# filter all sites that do not achieve a minimum coverage with the total number of reads
		thisCoverage = int(lineSplit[9].split(":")[4]) + int(lineSplit[9].split(":")[5])  # calculate the sum to only respect the reads that exceeded the quality threshold
		if thisCoverage < minimumNucleotideCoverage:
			filteredVariantsCoverageNucMin += 1
			filter = True
		
	# frequency threshold
	if tumorFreqThreshold > 0:
		tumorFreq = float(lineSplit[9].split(":")[6].split("%")[0])
		if tumorFreq <= tumorFreqThreshold:
			filteredVariantsFreq += 1
			filter = True
	
	#remove common variants
	if "y" in removeCommon:
		if ("COMMON=1" in line) or ("G5" in line):
			filteredVariantsCommon += 1
			filter = True
							
	if filter:
		filteredVariants += 1
		continue
	outfile.write(line)
	
	if ("COMMON=1" in line) or ("G5" in line):
		countCommon += 1
	if "HET=1" in line:
		countHeterozygous += 1
	if "HOM=1" in line:
		countHomozygous += 1
		
	refAllel = lineSplit[3]
	altAllel = lineSplit[4]
	
	if (len(refAllel) > 1) or (len(altAllel) > 1):
		countIndel += 1
	else:
		countSNP += 1
	
infile.close()
outfile.close()

print "%s variants left, %s of them heterozygous, %s homozygous. %s are annotated as common. %s SNPs and %s InDels. Filtered %s of %s variants (%.3f%%)." %((allVariants-filteredVariants),countHeterozygous,countHomozygous,countCommon,countSNP,countIndel,filteredVariants,allVariants,(float(filteredVariants)/float(max(1,allVariants)))*100.0)
print "VarScan default: %s, pvalue: %s, minVarCoverage: %s, minNucCoverage: %s, homopolymer: %s, strand bias: %s, synonomous: %s, tumorVarFreq: %s, common: %s." %(filteredVariantsVarScanDefault,filteredVariantsPVal,filteredVariantsCoverageVarMin,filteredVariantsCoverageNucMin,filteredVariantsHP,filteredVariantsStrand,filteredVariantsSyn,filteredVariantsFreq,filteredVariantsCommon)
		


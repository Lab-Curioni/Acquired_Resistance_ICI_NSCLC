#!/usr/bin/env python

'''
Combines dgidb query and variant or cnv annotation
Further, include pathway information, cBioportal information, and clinicalTrial information.
(Thus, script queryClinicalTrials.py must have been applied before)
Franziska Singer, March 2016
'''

import sys
import numpy as np
import os
import re

#convert copy number to copy number type (Discrete copy-number calls from the RAE algorithm. -2 = putative homozygous deletion; -1 = putative hemizygous deletion; 0 = copy-neutral; 1 = low-level gain; 2 = high-level amplification.)
def getCallType(thisCopyNumber):
	if thisCopyNumber < 2:
		if thisCopyNumber <= 0.5:
			return "-2"
		else:
			return "-1"
	if thisCopyNumber == 2:
		return "0"
	if thisCopyNumber > 2:
		if thisCopyNumber >= 3.5:
			return "2"
		else:
			return "1"


if len(sys.argv) <= 1:
	print "Combine dgidb and variant information."
	print "Usage: python combine_dgidb_cnvs_includeDetails_woCbioportal.py [variantInfofile] [dgidb_geneList] [clinicalTrials] [pathwayDB] [outfile]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Combine dgidb and variant information."
	print "Usage: python combine_dgidb_cnvs_includeDetails_woCbioportal.py [variantInfofile] [dgidb_geneList] [clinicalTrials] [pathwayDB] [outfile]"
	sys.exit(1)

variantInfoFile = sys.argv[1]  # this file is the result of processBICSeq2.py (or processExcavator, for whole exome)
dgidbFile_categ = sys.argv[2]			# this file contains gene category information
clinicalTrialsFile = sys.argv[3]		# dgidb complete table including clinicalTrials information, result of queryClinicalTrials.py
pathwayDB = sys.argv[4]
outfileName = sys.argv[5]

infileDGIDB_categories = open(dgidbFile_categ,'r')
infileDGIDB_categories.readline()  # skip first line

# gene category information
dictDGIDB_categ = {}  # gene names to information

dgidbGenes = 0
for lineDGIDB in infileDGIDB_categories: # format: Gene	GeneName	Category
	lineSplit = lineDGIDB.strip().split("\t")
	gene = lineSplit[0]
	
	if gene in dictDGIDB_categ.keys():
		print "Error! %s already contained!" %(gene)
		continue
		
	dictDGIDB_categ[gene] = [lineSplit[1],lineSplit[2]]
	dgidbGenes += 1
	
infileDGIDB_categories.close()

# gene drug and trials information
infileDGIDB = open(clinicalTrialsFile,'r')
dictDGIDB = {}  # gene names to information on drugs and clinical trials, including support of drug gene interaction
				# [score,type,clinicalTrials]

drugNum = 0
firstSplit = infileDGIDB.readline().strip().split("\t")
indexScore = 17
indexType = 18
indexTrialInfo = 19

for i in range (0,len(firstSplit)):
	if "Score" in firstSplit[i]:
		indexScore = i
	if "Type" in firstSplit[i]:
		indexType = i
	if "ClinicalTrials" in firstSplit[i]:
		indexTrialInfo = i

for lineDGIDB in infileDGIDB:
	lineSplit = lineDGIDB.strip().split("\t")
	gene = lineSplit[0]
	drug = lineSplit[1]
	score = lineSplit[indexScore]
	drugType = lineSplit[indexType]
	clinicalTrialsInfo = lineSplit[indexTrialInfo]	
	
	if gene not in dictDGIDB.keys():
		dictDGIDB[gene] = {}
	if drug not in dictDGIDB[gene].keys():
		dictDGIDB[gene][drug] = []
	else:
		print "Warning! Frug %s already contained for gene %s!" %(drug,gene)
		continue
		
	dictDGIDB[gene][drug].append(score)
	dictDGIDB[gene][drug].append(drugType)
	dictDGIDB[gene][drug].append(clinicalTrialsInfo)
	drugNum += 1
	
infileDGIDB.close()

# sanity check: category and clinicalTrials file should include the same number of genes

if len(dictDGIDB.keys()) != len(dictDGIDB_categ.keys()):
	print "Error! Different number of genes!"
else:
	print "Number of different dgidb genes: %s." %(len(dictDGIDB_categ.keys())) 
	
print "Found %s drugs for dgidb genes.\n" %(drugNum)

# pathway information
infilePathways = open(pathwayDB,'r')
dictPathways = {}

for linePW in infilePathways:
	lineSplit = linePW.strip().split()
	geneName = lineSplit[0].strip("\"")  # strip to remove the " symbols
	allPWs = "".join(lineSplit[1:])
	#print allPWs
	dictPathways[geneName] = allPWs.strip("\"")

infilePathways.close()

# finally parse existing variant table and create new file with more detailed information per gene

infileVariant = open(variantInfoFile,'r')
outfile = open(outfileName,'w')

line1 = infileVariant.readline()
isWGS = False

if "pvalue" in line1: # then WGS (bicSeq2)
	isWGS = True

if isWGS:
	outfile.write("Gene\tChromosome\tStart\tEnd\tCopyNumber\tCallType\tpvalue\tlog2.copyRatio\tClinicalTrials_cancerType (phase,isRecruiting (y/n))\tClinicalTrials_notCancerType (phase,isRecruiting (y/n))\tPathway\tDGIDB-drugs(Score,Type)\tDGIDB-geneName\tDGIDB-categories\n")
else:
	outfile.write("Gene\tChromosome\tStart\tEnd\tCopyNumber\tCallType\tProbability\tClinicalTrials_cancerType (phase,isRecruiting (y/n))\tClinicalTrials_notCancerType (phase,isRecruiting (y/n))\tPathway\tDGIDB-drugs(Score,Type)\tDGIDB-geneName\tDGIDB-categories\n")

matchedDGIDB = []
pathwayInfo = 0
for line in infileVariant:
	lineSplit = line.strip().split("\t")
	genes = lineSplit[len(lineSplit)-1]
	
	for dgidbGene in dictDGIDB_categ.keys():
		if dgidbGene in genes:
			foundMatch = False
			if ";" in genes:  # split and check if indeed the correct gene
				semiSplit = genes.split(";")
				for temp in semiSplit:
					if str(temp.strip()) == str(dgidbGene.strip()):  # equal, no additional signs
						foundMatch = True
			else:
				if len(genes.strip()) == len(dgidbGene.strip()):  # equal, no additional signs
					foundMatch = True
			
			#if not foundMatch:
			#	print "NoMatch: dgidb: %s ; genes: %s" %(dgidbGene,genes)

			if foundMatch:
				if dgidbGene not in matchedDGIDB:
					matchedDGIDB.append(dgidbGene)
					
				# pathway info
				pathwayString = ""
				if dgidbGene in dictPathways.keys():
					pathwayInfo += 1
					pathwayString = dictPathways[dgidbGene]
				else:
					pathwayString = "."
					
				# drug and clincical trials info
				drugNames = ""
				clinicalTrialCancerType = ""
				clincialTrialNotCancerType = ""
				if dgidbGene not in dictDGIDB.keys():
					print "Error. Gene %s not conatined in dgidb clincial trials dict!" %(dgidbGene)
					continue
				else:
					for drugInteractions in dictDGIDB[dgidbGene].keys():
						drugNames += drugInteractions + "(" + dictDGIDB[dgidbGene][drugInteractions][0] + "," + dictDGIDB[dgidbGene][drugInteractions][1] + ");"
						clinTrialInfo = dictDGIDB[dgidbGene][drugInteractions][2]
						
						foundCT = False
						foundNotCT = False
						numTrialsCT = 0
						numTrialsNotCT = 0
						clinicalTrialCancerType_temp = ""
						clincialTrialNotCancerType_temp_1 = "" # for the different phases. Only report trial IDs for the highest phase found for this drug
						clincialTrialNotCancerType_temp_2 = ""
						clincialTrialNotCancerType_temp_3 = ""
						if ";" in clinTrialInfo:
							clinSplit = clinTrialInfo.strip().split(";")
							for clinEntry in clinSplit:
								trialInfo = clinEntry.split(",")
								if "yes" in trialInfo[3]: 
									clinicalTrialCancerType_temp += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
									foundCT = True
									numTrialsCT += 1
								else:
									if "3" in trialInfo[2]:
										clincialTrialNotCancerType_temp_3 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
									elif "2" in trialInfo[2]:
										clincialTrialNotCancerType_temp_2 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
									else:
										clincialTrialNotCancerType_temp_1 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
									foundNotCT = True
									numTrialsNotCT += 1
							if not foundCT:
								clinicalTrialCancerType_temp += ".;"
							if not foundNotCT:
								clincialTrialNotCancerType_temp_1 += ".;"
								clincialTrialNotCancerType_temp_2 += ".;"
								clincialTrialNotCancerType_temp_3 += ".;"
							
						else:
							if "," in clinTrialInfo:
								trialInfo = clinTrialInfo.split(",")
								if "yes" in trialInfo[3]:
									clinicalTrialCancerType_temp += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
									foundCT = True
									numTrialsCT += 1
								else:
									if "3" in trialInfo[2]:
										clincialTrialNotCancerType_temp_3 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
									elif "2" in trialInfo[2]:
										clincialTrialNotCancerType_temp_2 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
									else:
										clincialTrialNotCancerType_temp_1 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
									foundNotCT = True
									numTrialsNotCT += 1
						if not foundCT:
							clinicalTrialCancerType_temp += ".;"
						if not foundNotCT:
							clincialTrialNotCancerType_temp_1 += ".;"
							clincialTrialNotCancerType_temp_2 += ".;"
							clincialTrialNotCancerType_temp_3 += ".;"
						# now append everything to string for outfile
						clinicalTrialCancerType += drugInteractions + "(" + str(numTrialsCT) + "):" + clinicalTrialCancerType_temp
						if len(clincialTrialNotCancerType_temp_3) > 0:
							clincialTrialNotCancerType += drugInteractions + "(" + str(numTrialsNotCT) + "):" + clincialTrialNotCancerType_temp_3
						elif len(clincialTrialNotCancerType_temp_2) > 0:
							clincialTrialNotCancerType += drugInteractions + "(" + str(numTrialsNotCT) + "):" + clincialTrialNotCancerType_temp_2
						else:
							clincialTrialNotCancerType += drugInteractions + "(" + str(numTrialsNotCT) + "):" + clincialTrialNotCancerType_temp_1
					
				# write out and inlcude category info
				cnvType = ""
				if isWGS:
					thisCopyNumber = float(lineSplit[5]) # without cBioportal, we have to include the type conversion for the WGS CNVs here
					cnvType = getCallType(thisCopyNumber)
				
				if isWGS:
					outfile.write(dgidbGene + "\t" + lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[2] + "\t" + lineSplit[5] + "\t" + cnvType + "\t" + lineSplit[4] + "\t" + lineSplit[3] + "\t" + clinicalTrialCancerType + "\t" + clincialTrialNotCancerType + "\t" + pathwayString + "\t" + drugNames + "\t" + dictDGIDB_categ[dgidbGene][0] + "\t" + dictDGIDB_categ[dgidbGene][1] + "\n")
				else:
					outfile.write(dgidbGene + "\t" + lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[2] + "\t" + lineSplit[len(lineSplit)-4] + "\t" + lineSplit[len(lineSplit)-3] + "\t" + lineSplit[len(lineSplit)-2] + "\t" + clinicalTrialCancerType + "\t" + clincialTrialNotCancerType + "\t" + pathwayString + "\t" + drugNames + "\t" + dictDGIDB_categ[dgidbGene][0] + "\t" + dictDGIDB_categ[dgidbGene][1] + "\n")

infileVariant.close()
outfile.close()

print "Matched %s of %s genes from dgidb. Found pathway information for %s genes." %(len(matchedDGIDB),dgidbGenes,pathwayInfo)

#!/usr/bin/env python

'''
Combine the snv overview table with results from clinical trials and dgidb query
Franziska Singer, March 2016
'''

import sys
import numpy as np
import os
import re


if len(sys.argv) <= 1:
	print "Get annotated genes from a vcf file."
	print "Usage: python combine_dgidb_snvs_includeDetails_woCbioportal.py [variantInfofile] [dgidb_geneCategoryList] [clinicalTrials] [outfile]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Get annotated genes from a vcf file."
	print "Usage: python combine_dgidb_snvs_includeDetails_woCbioportal.py [variantInfofile] [dgidb_geneCategoryList] [clinicalTrials] [outfile]"
	sys.exit(1)

	
variantInfoFile = sys.argv[1]  # this file is the result of processAnnotatedVCF_includeDetails.py
dgidbFile_categ = sys.argv[2]			# this file contains gene category information
clinicalTrialsFile = sys.argv[3]		# dgidb complete table including clinicalTrials information, result of queryClinicalTrials.py
outfileName = sys.argv[4]

infileDGIDB_categories = open(dgidbFile_categ,'r')
infileDGIDB_categories.readline()  # skip first line

# gene category information
dictDGIDB_categ = {}  # gene names to information

dgidbGenes = 0
for lineDGIDB in infileDGIDB_categories:  # format: Gene	GeneName	Category
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
		print "Warning! Drug %s already contained for gene %s!" %(drug,gene)
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


# parse existing variant table and create new file with more detailed information per gene

infile = open(variantInfoFile,'r')
outfile = open(outfileName,'w')

outfile.write("Gene\tChromosome\tPosition\tChange\tFrequency\tAnnotation\tAnnotatedImpact\tVariantID\tFunctionalAnnotation\tClinicalTrials_cancerType (phase,isRecruiting (y/n))\tClinicalTrials_notCancerType (phase,isRecruiting (y/n))\tPathway\tDGIDB-drugs(Score,Type)\tDGIDB-geneName\tDGIDB-categories\n")

allVariants = 0

infile.readline()
matchedDGIDB = []

for line in infile:
	lineSplit = line.strip().split("\t")
	genes = lineSplit[4]
	position = int(lineSplit[1])
	
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
				#print "NoMatch: dgidb: %s ; genes: %s" %(dgidbGene,genes)

			if foundMatch:
				if dgidbGene not in matchedDGIDB:
					matchedDGIDB.append(dgidbGene)
					
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
							for clinEntry in clinSplit: # format: (ID,isRecruiting y/n,Phase,isCancerTypeSpecific yes/no)
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
				# already in file: Chromosome	Position	Change	Frequency	Genes	Annotation	VariantID	FunctionalAnnotation	Pathway VariantImpact
				# table header: Gene\tChromosome\tPosition\tChange\tFrequency\tAnnotation\tVariantID\tAnnotatedImpact\tFunctionalAnnotation\tClinicalTrials_recruiting\tClinicalTrials_completed\tCBioportal_cancerType_GENE(mutated/allCases/frequency)\tCBioportal_cancerType_VARIANT(mutated/allCases/frequency)\tCBioportal_studyMax_GENE(mutated/allCases/frequency)\tCBioportal_studyMax_VARIANT(mutated/allCases/frequency)\tPathway\tDGIDB-drugs(Score,Type)\tDGIDB-geneName\tDGIDB-categories\n"
				outfile.write(dgidbGene + "\t" + lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[2] + "\t" + lineSplit[3] + "\t" + lineSplit[5] + "\t" + lineSplit[9] + "\t" + lineSplit[6] + "\t" + lineSplit[7] + "\t" + clinicalTrialCancerType + "\t" + clincialTrialNotCancerType + "\t" + lineSplit[8] + "\t" + drugNames + "\t" + dictDGIDB_categ[dgidbGene][0] + "\t" + dictDGIDB_categ[dgidbGene][1] + "\n")
				
infile.close()
outfile.close()

print "Matched %s of %s genes from dgidb." %(len(matchedDGIDB),dgidbGenes)

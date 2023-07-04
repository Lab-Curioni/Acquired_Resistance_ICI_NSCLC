#!/usr/bin/env python

'''
Combine the snv overview table with results from cBioportal query, clinical trials information, and dgidb
Franziska Singer, March 2016
'''

import sys
import numpy as np
import os
import re


if len(sys.argv) <= 1:
	print "Get annotated genes from a vcf file."
	print "Usage: python combine_dgidb_snvs_includeDetails.py [variantInfofile] [dgidb_geneList] [cBioportal_max] [cBioportal_cancerTypeSpecific] [clinicalTrials] [outfile]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Get annotated genes from a vcf file."
	print "Usage: python combine_dgidb_snvs_includeDetails.py [variantInfofile] [dgidb_geneList] [cBioportal_max] [cBioportal_cancerTypeSpecific] [clinicalTrials] [outfile]"
	sys.exit(1)

	
variantInfoFile = sys.argv[1]  # this file is the result of processAnnotatedVCF_includeDetails.py
dgidbFile_categ = sys.argv[2]			# this file contains gene category information
inputCBioportal_max = sys.argv[3]
inputCBioportal_Cancertype = sys.argv[4]
clinicalTrialsFile = sys.argv[5]		# dgidb complete table including clinicalTrials information, result of queryClinicalTrials.py
outfileName = sys.argv[6]

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
headerLine = infileDGIDB.readline()
scoreColumn = 0
clinTrialColumn = 0
drugTypeColumn = 0

headerSplit = headerLine.strip().split("\t")
for headPos in range(0,len(headerSplit)):
	if headerSplit[headPos] == "Score":
		scoreColumn = headPos
		continue
	if headerSplit[headPos] == "Type":
		drugTypeColumn = headPos
		continue
	if headerSplit[headPos] == "ClinicalTrials (ID,isRecruiting y/n,Phase,isCancerTypeSpecific yes/no)":
		clinTrialColumn = headPos
		continue
	
for lineDGIDB in infileDGIDB:
	lineSplit = lineDGIDB.strip().split("\t")
	gene = lineSplit[0]
	drug = lineSplit[1]
	score = lineSplit[scoreColumn]
	drugType = lineSplit[drugTypeColumn]
	clinicalTrialsInfo = lineSplit[clinTrialColumn]	
	
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

# cBioportal information

infileCBP_ct= open(inputCBioportal_Cancertype,'r')
infileCBP_max = open(inputCBioportal_max,'r')

dictCbioportal_ct = {}   # gene to position to [studyID,#casesStudyCancerType,#casesMutCancerType,#casesWithVariantcancerType,frequencyMut,frequencyVariant]
					  # one entry for each position

dictCbioportal_max = {}   # gene to position to [studyID_max,#casesStudy,#casesMut,#casesWithVariant,frequencyMut,frequencyVariant]
					  # one entry for each position
					  
infileCBP_ct.readline()
for lineCBP_ct in infileCBP_ct: #study	allCases	gene	position	casesWithMutatedGene	casesWithVariant	frequencyMutatedGene	frequencyVariant
	lineSplit = lineCBP_ct.strip().split("\t")
	geneName = lineSplit[2]
	if geneName not in dictCbioportal_ct.keys():
		dictCbioportal_ct[geneName] = {}
	position = int(lineSplit[3])
	if position not in dictCbioportal_ct[geneName].keys():
		dictCbioportal_ct[geneName][position] = []
	dictCbioportal_ct[geneName][position].append(lineSplit[0])
	dictCbioportal_ct[geneName][position].append(lineSplit[1])
	for i in range(4,8):
		dictCbioportal_ct[geneName][position].append(lineSplit[i])

infileCBP_ct.close()

infileCBP_max.readline()
for lineCBP_mx in infileCBP_max: #study	allCases	gene	position	casesWithMutatedGene	casesWithVariant	frequencyMutatedGene	frequencyVariant
	lineSplit = lineCBP_mx.strip().split("\t")
	geneName = lineSplit[2] # note: always two lines per variant position of this gene! One for geneMax (multiple times same entry becuase equal for all positions of this gene), one for variantMax
	if geneName not in dictCbioportal_max.keys():
		dictCbioportal_max[geneName] = {}
	position = int(lineSplit[3])
	if position not in dictCbioportal_max[geneName].keys():
		dictCbioportal_max[geneName][position] = []
	tempArr = [lineSplit[0],lineSplit[1],lineSplit[4],lineSplit[5],lineSplit[6],lineSplit[7]]
	dictCbioportal_max[geneName][position].append(tempArr)

infileCBP_max.close()

# finally parse existing variant table and create new file with more detailed information per gene

infile = open(variantInfoFile,'r')
outfile = open(outfileName,'w')
outfileIndependent = open(outfileName + "_dgidbIndependent.txt",'w')

outfile.write("Gene\tChromosome\tPosition\tChange\tFrequency\tAnnotation\tAnnotatedImpact\tVariantID\tFunctionalAnnotation\tClinicalTrials_cancerType (phase,isRecruiting (y/n))\tClinicalTrials_notCancerType (phase,isRecruiting (y/n))\tCBioportal_cancerType_GENE(mutated/allCases/frequency)\tCBioportal_cancerType_VARIANT(mutated/allCases/frequency)\tCBioportal_studyMax_GENE(mutated/allCases/frequency)\tCBioportal_studyMax_VARIANT(mutated/allCases/frequency)\tPathway\tDGIDB-drugs(Score,Type)\tDGIDB-geneName\tDGIDB-categories\n")
outfileIndependent.write("Gene\tChromosome\tPosition\tChange\tFrequency\tAnnotation\tAnnotatedImpact\tVariantID\tFunctionalAnnotation\tClinicalTrials_cancerType (phase,isRecruiting (y/n))\tClinicalTrials_notCancerType (phase,isRecruiting (y/n))\tCBioportal_cancerType_GENE(mutated/allCases/frequency)\tCBioportal_cancerType_VARIANT(mutated/allCases/frequency)\tCBioportal_studyMax_GENE(mutated/allCases/frequency)\tCBioportal_studyMax_VARIANT(mutated/allCases/frequency)\tPathway\tDGIDB-drugs(Score,Type)\tDGIDB-geneName\tDGIDB-categories\n")

allVariants = 0

infile.readline()
matchedDGIDB = []

for line in infile:
	lineSplit = line.strip().split("\t")
	genes = lineSplit[4]
	position = int(lineSplit[1])
	
	includeInIndependentFile = True
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
				includeInIndependentFile = False
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
									# not cancer-type specific, check which phase
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
					
				# cbioportal info
				
				cbioCTString_gene = ""
				cbioMaxString_gene = ""
				cbioCTString_var = ""
				cbioMaxString_var = ""
				
				if dgidbGene in dictCbioportal_ct.keys():
					if position in dictCbioportal_ct[dgidbGene].keys():
						cbioportal_ct = dictCbioportal_ct[dgidbGene][position] #[studyID,#casesStudyCancerType,#casesMutCancerType,#casesWithVariantcancerType,frequencyMut,frequencyVariant]
						 
						if float(cbioportal_ct[4]) == 0.0:
							cbioCTString_gene += "0"
							cbioCTString_var += "0"
						else:
							cbioCTString_gene += "%s(%s/%s/%.1f%%)"%(cbioportal_ct[0],cbioportal_ct[2],cbioportal_ct[1],(float(cbioportal_ct[4])*100.0))
							cbioCTString_var += "%s(%s/%s/%.1f%%)"%(cbioportal_ct[0],cbioportal_ct[3],cbioportal_ct[1],(float(cbioportal_ct[5])*100.0))
					
				else:	
					print "Warning! Gene %s not contained in cBioportal cancer type-specific dictionary!" %(dgidbGene)
					cbioCTString_gene += "NA"
					cbioCTString_var += "NA"
					
				if dgidbGene in dictCbioportal_max.keys():
					cbioportal_max = dictCbioportal_max[dgidbGene]   #[studyID,#casesStudyCancerType,#casesMutCancerType,#casesWithVariantcancerType,frequencyMut,frequencyVariant]
					
					maxIndexGene = 0
					maxPosGene = 0
					maxFreqGene = 0.0
					
					maxIndexVar = 0
					maxFreqVar = 0.0
					foundPosi = False
					for posiDict in dictCbioportal_max[dgidbGene].keys():  # due to multiple entries in maximum table, find the maximum over all positions
						cbioportal_max = dictCbioportal_max[dgidbGene][posiDict]
						if posiDict == position:
							foundPosi = True
							for index in range(0,len(cbioportal_max)):
								if float(cbioportal_max[index][5]) > maxFreqVar:
									maxIndexVar = index
									maxFreqVar = float(cbioportal_max[index][5])
							
						for index in range(0,len(cbioportal_max)):
							if float(cbioportal_max[index][4]) > maxFreqGene:
								maxIndexGene = index
								maxPosGene = posiDict
								maxFreqGene = float(cbioportal_max[index][4])
								
								
					if maxFreqGene == 0.0:
						cbioMaxString_gene += "0"
						cbioMaxString_var += "0"
					else:
						cbioMaxString_gene += "%s(%s/%s/%.1f%%)"%(dictCbioportal_max[dgidbGene][maxPosGene][maxIndexGene][0],dictCbioportal_max[dgidbGene][maxPosGene][maxIndexGene][2],dictCbioportal_max[dgidbGene][maxPosGene][maxIndexGene][1],(float(dictCbioportal_max[dgidbGene][maxPosGene][maxIndexGene][4])*100.0))
						if foundPosi:
							cbioMaxString_var += "%s(%s/%s/%.1f%%)"%(dictCbioportal_max[dgidbGene][position][maxIndexVar][0],dictCbioportal_max[dgidbGene][position][maxIndexVar][3],dictCbioportal_max[dgidbGene][position][maxIndexVar][1],(float(dictCbioportal_max[dgidbGene][position][maxIndexVar][5])*100.0))
						else:
							print "Did not find position %s for gene %s!" %(position, dgidbGene)
				else:
					cbioMaxString_gene += "NA"
					cbioMaxString_var += "NA"
					
				# write out and inlcude category info
				# already in file: Chromosome	Position	Change	Frequency	Genes	Annotation	VariantID	FunctionalAnnotation	Pathway VariantImpact
				# table header: Gene\tChromosome\tPosition\tChange\tFrequency\tAnnotation\tAnnotatedImpact\tVariantID\tFunctionalAnnotation\tClinicalTrials_recruiting\tClinicalTrials_completed\tCBioportal_cancerType_GENE(mutated/allCases/frequency)\tCBioportal_cancerType_VARIANT(mutated/allCases/frequency)\tCBioportal_studyMax_GENE(mutated/allCases/frequency)\tCBioportal_studyMax_VARIANT(mutated/allCases/frequency)\tPathway\tDGIDB-drugs(Score,Type)\tDGIDB-categories\n"
				outfile.write(dgidbGene + "\t" + lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[2] + "\t" + lineSplit[3] + "\t" + lineSplit[5] + "\t" + lineSplit[9] + "\t" + lineSplit[6] + "\t" + lineSplit[7] + "\t" + clinicalTrialCancerType + "\t" + clincialTrialNotCancerType + "\t" + cbioCTString_gene + "\t" + cbioCTString_var + "\t" + cbioMaxString_gene + "\t" + cbioMaxString_var + "\t" + lineSplit[8] + "\t" + drugNames + "\t" + dictDGIDB_categ[dgidbGene][0] + "\t" + dictDGIDB_categ[dgidbGene][1] + "\n")
				outfileIndependent.write(dgidbGene + "\t" + lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[2] + "\t" + lineSplit[3] + "\t" + lineSplit[5] + "\t" + lineSplit[9] + "\t" + lineSplit[6] + "\t" + lineSplit[7] + "\t" + clinicalTrialCancerType + "\t" + clincialTrialNotCancerType + "\t" + cbioCTString_gene + "\t" + cbioCTString_var + "\t" + cbioMaxString_gene + "\t" + cbioMaxString_var + "\t" + lineSplit[8] + "\t" + drugNames + "\t" + dictDGIDB_categ[dgidbGene][0] + "\t" + dictDGIDB_categ[dgidbGene][1] + "\n")
				
	if includeInIndependentFile:
		# include in independent outfile
		
		#cbioportal:
			
		cbioCTString_gene = ""
		cbioMaxString_gene = ""
		cbioCTString_var = ""
		cbioMaxString_var = ""
		
		genesSplit = genes.split(";")
		for myGene in genesSplit:
		
			if myGene in dictCbioportal_ct.keys():
				if len(genesSplit) > 1:
					cbioCTString_gene += "%s:" %(myGene)
					cbioCTString_var += "%s:" %(myGene)
				
				if position in dictCbioportal_ct[myGene].keys():
					cbioportal_ct = dictCbioportal_ct[myGene][position] #[studyID,#casesStudyCancerType,#casesMutCancerType,#casesWithVariantcancerType,frequencyMut,frequencyVariant]
					 
					if float(cbioportal_ct[4]) == 0.0:
						cbioCTString_gene += "0"
						cbioCTString_var += "0"
					else:
						cbioCTString_gene += "%s(%s/%s/%.1f%%)"%(cbioportal_ct[0],cbioportal_ct[2],cbioportal_ct[1],(float(cbioportal_ct[4])*100.0))
						cbioCTString_var += "%s(%s/%s/%.1f%%)"%(cbioportal_ct[0],cbioportal_ct[3],cbioportal_ct[1],(float(cbioportal_ct[5])*100.0))		
				
				if len(genesSplit) > 1:
					cbioCTString_gene += ";" 
					cbioCTString_var += ";" 
				
			else:
				#print "Warning! Gene %s not contained in cBioportal cancer type-specific dictionary!" %(myGene)
				if len(genesSplit) > 1:
					cbioCTString_gene += "%s:NA;" %(myGene)
					cbioCTString_var += "%s:NA;" %(myGene)
				else:
					cbioCTString_gene += "NA"
					cbioCTString_var += "NA"
				
				
			if myGene in dictCbioportal_max.keys():
				if len(genesSplit) > 1:
					cbioMaxString_gene += "%s:" %(myGene)
					cbioMaxString_var += "%s:" %(myGene)
					
				cbioportal_max = dictCbioportal_max[myGene]   #[studyID,#casesStudyCancerType,#casesMutCancerType,#casesWithVariantcancerType,frequencyMut,frequencyVariant]
				
				maxIndexGene = 0
				maxPosGene = 0
				maxFreqGene = 0.0
				
				maxIndexVar = 0
				maxFreqVar = 0.0
				foundPosi = False
				for posiDict in dictCbioportal_max[myGene].keys():  # due to multiple entries in maximum table, find the maximum over all positions
					cbioportal_max = dictCbioportal_max[myGene][posiDict]
					if posiDict == position:
						foundPosi = True
						for index in range(0,len(cbioportal_max)):
							if float(cbioportal_max[index][5]) > maxFreqVar:
								maxIndexVar = index
								maxFreqVar = float(cbioportal_max[index][5])
						
					for index in range(0,len(cbioportal_max)):
						if float(cbioportal_max[index][4]) > maxFreqGene:
							maxIndexGene = index
							maxPosGene = posiDict
							maxFreqGene = float(cbioportal_max[index][4])
							
							
				if maxFreqGene == 0.0:
					cbioMaxString_gene += "0"
					cbioMaxString_var += "0"
				else:
					cbioMaxString_gene += "%s(%s/%s/%.1f%%)"%(dictCbioportal_max[myGene][maxPosGene][maxIndexGene][0],dictCbioportal_max[myGene][maxPosGene][maxIndexGene][2],dictCbioportal_max[myGene][maxPosGene][maxIndexGene][1],(float(dictCbioportal_max[myGene][maxPosGene][maxIndexGene][4])*100.0))
					if foundPosi:
						cbioMaxString_var += "%s(%s/%s/%.1f%%)"%(dictCbioportal_max[myGene][position][maxIndexVar][0],dictCbioportal_max[myGene][position][maxIndexVar][3],dictCbioportal_max[myGene][position][maxIndexVar][1],(float(dictCbioportal_max[myGene][position][maxIndexVar][5])*100.0))
					else:
						print "Did not find position %s for gene %s!" %(position, myGene)
				
				if len(genesSplit) > 1:
					cbioMaxString_gene += ";" 
					cbioMaxString_var += ";" 
					
			else:
				if len(genesSplit) > 1:
					cbioMaxString_gene += "%s:NA;" %(myGene)
					cbioMaxString_var += "%s:NA;" %(myGene)
				else:
					cbioMaxString_gene += "NA"
					cbioMaxString_var += "NA"
			
		outfileIndependent.write(genes + "\t" + lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[2] + "\t" + lineSplit[3] + "\t" + lineSplit[5] + "\t" + lineSplit[9] + "\t" + lineSplit[6] + "\t" + lineSplit[7] + "\tNA\tNA\t" + cbioCTString_gene + "\t" + cbioCTString_var + "\t" + cbioMaxString_gene + "\t" + cbioMaxString_var + "\t" + lineSplit[8] + "\tNA\tNA\tNA\n")
					
				
infile.close()
outfile.close()
outfileIndependent.close()

print "Matched %s of %s genes from dgidb." %(len(matchedDGIDB),dgidbGenes)

#!/usr/bin/env python

'''
seiler_tcga_2017
get drug information from clinicalAnnotation.txt
Franziska Singer, July 2018
'''

import sys
import numpy as np
import os
import re
import argparse
from collections import OrderedDict
from natsort import natsorted

'''
function definitions
'''


def parseInputFile(sample,name,dictSamples,dictDrugs,drugsToVariants):
	
	print("Sample %s. File: %s" %(name,sample))

	infile = open(sample,'r')
	header = infile.readline()
	(index_drugs,index_ct,index_nct,index_gene) = getColumnIndex(header)

	drugsInSample = []
	ct_drugs = 0
	nct_drugs = 0
	other_drugs = 0

	for line in infile:
		lineSplit = line.strip().split("\t")
		drugs = lineSplit[index_drugs].split(";") #OXAZEPAM(5,positive allosteric modulator,potentiator);CLOTIAZEPAM(5,potentiator)

		for drugEntry in drugs:
			drugSplit = drugEntry.split("(")
			drug = drugSplit[0]
			score = drugSplit[1].split(",")[0].strip()
			if len(drugSplit) > 2: # more brackets included , drug name containes brackets
				drug = "(".join(drugSplit[0:-1])
				score = drugSplit[-1].split(",")[0].strip()
			
			gene = lineSplit[index_gene]
			if drug not in drugsToVariants.keys():
				drugsToVariants[drug] = {}
				drugsToVariants[drug][name] = {}
			if name not in drugsToVariants[drug].keys():
				drugsToVariants[drug][name] = {}
			if gene not in drugsToVariants[drug][name].keys():
				drugsToVariants[drug][name][gene] = [1,score]
			else:
				drugsToVariants[drug][name][gene][0] += 1 # score stays the same for this drug gene combination

			trials_ct = lineSplit[index_ct]
			trials_nct = lineSplit[index_nct]
			
			alreadySeen = False
			for seenDrug in drugsInSample:  # slow if lots of drugs, but we need to ensure that drug names contained in other drug names are handled correctly
				if drug == seenDrug:
					alreadySeen = True
					break

			if alreadySeen:  # only if we have not yet seen this drug before, we need to look at it. Otherwise we have already included it where necessary
				continue

			drugsInSample.append(drug)

			isCT = False
			isNCT = False
			if drug in trials_ct: # cancer_type specific trials; note that drugs that are name parts of other drugs are marked, too!
				ct_drugs += 1
				isCT = True
			if drug in trials_nct:
				nct_drugs += 1
				isNCT = True
			if (isCT == False) and (isNCT == False):
				other_drugs += 1
			
			if drug not in dictDrugs.keys():
				dictDrugs[drug] = []  # array[0]= dgidbScore, array[1] = type (ct/nct/other) , arry[2-n] = sampleNames
				dictDrugs[drug].append(score)
				if isCT == True:
					dictDrugs[drug].append("ct")
				elif isNCT == True:
					dictDrugs[drug].append("nct")
				else:
					dictDrugs[drug].append("other")

			dictDrugs[drug].append(name)
		

	infile.close()
	dictSamples[name] = [len(drugsInSample),ct_drugs,nct_drugs,other_drugs]

	return (dictSamples,dictDrugs,drugsToVariants)

def getColumnIndex(header):
	headerSplit = header.strip().split("\t")

	index_drugs = -1
	index_ct = -1
	index_nct = -1
	index_gene = -1
	
	for pos in range(0,len(headerSplit)):
		if args.colName_drugs == headerSplit[pos]:
			index_drugs = pos
		if args.colName_trials_ct == headerSplit[pos]:
			index_ct = pos
		if args.colName_trials_nct == headerSplit[pos]:
			index_nct = pos
		if args.colName_gene == headerSplit[pos]:
			index_gene = pos

	if (index_drugs == -1) or (index_ct == -1) or (index_nct == -1) or (index_gene == -1):
		print("Error! Could not match all input columns in header %s." %(header))
		sys.exit(1)
	
	return (index_drugs,index_ct,index_nct,index_gene)



parser = argparse.ArgumentParser(description='Get drug information from files of form [sample].clinicalAnnotation.txt.')
parser.add_argument('--inputDir', dest='inputDir', required=True, help='Input directory with files of format [sample].clinicalAnnotation.txt.')
parser.add_argument('--outFileTag', dest='outFileTag', required=True, help='Name tag of the output files.')
parser.add_argument('--colName_drugs', dest='colName_drugs', required=True, help='Column name of column containing drug names.')
parser.add_argument('--colName_gene', dest='colName_gene', required=True, help='Column name of column containing gene names.')
parser.add_argument('--colName_trials_ct', dest='colName_trials_ct', required=True, help='Column name of column containing info on cancer type specific trials.')
parser.add_argument('--colName_trials_nct', dest='colName_trials_nct', required=True, help='Column name of column containing info on non cancer type specific trials.')
parser.add_argument('--fileEnding', dest='fileEnding', required=True, help='To get correct input files, specify the file ending of desired input files.')

args = parser.parse_args()

print("Parameters:\nInputDir: %s\noutFileTag: %s\ncolName_drugs: %s\ncolName_gene: %s\ncolName_trials_ct: %s\ncolName_trials_nct: %s\nfileEnding: %s\n" %(args.inputDir,args.outFileTag,args.colName_drugs,args.colName_gene,args.colName_trials_ct,args.colName_trials_nct,args.fileEnding))
# drug information
dictDrugs = {}

# sample information
dictSamples = {}
sampleNames = [] # for sorting
drugsToVariants = {}  # dictionary to match drugs to the underlying target genes and dgidb scores in each sample

for file in os.listdir(args.inputDir):
	sampleFile = "%s%s" %(args.inputDir,os.path.basename(file))
	if os.path.isfile(sampleFile) and sampleFile.endswith(args.fileEnding):
		name = os.path.basename(file).split("_")[0].split("-")[1]
		if name not in sampleNames:
			sampleNames.append(name)
			#dictSample[name] = [0,0,0] # all, ct, nct
		else:
			print("Error! Multiple files for %s!" %(name))
		(dictSamples,dictDrugs,drugsToVariants) = parseInputFile(sampleFile,name,dictSamples,dictDrugs,drugsToVariants)


sampleNames = natsorted(sampleNames)
outfile_samples = open(args.outFileTag + ".samples_drugNumbers.txt",'w')
outfile_samples.write("Sample\tAll_drugs\tDrugs_cancerType\tDrugs_not_cancerType\tDrugs_other\n")

outfile_drugs = open(args.outFileTag + ".drugs_sampleInfos.txt", 'w')
drug_headerLine = "Drug\tNumberSamples\tDGIdb_score\tCancer_related"

outfile_drugs_info = open(args.outFileTag + ".drugs_variantInfos.txt", 'w')
drug_info_headerLine = "Drug\tNumberSamples\tCancer_related"

#dictSampleIndex = {}
#index_temp = 4
for sample in sampleNames:
	if sample not in dictSamples.keys():
		print("Error! Sample %s not contained in sample names list " %(sample))
		continue
	drugNums = dictSamples[sample]
	outfile_samples.write("%s\t%s\t%s\t%s\t%s\n" %(sample,drugNums[0],drugNums[1],drugNums[2],drugNums[3]))

	#dictSampleIndex[sample] = index_temp
	#index_temp += 1
	drug_headerLine = drug_headerLine + "\t" + sample
	drug_info_headerLine = drug_info_headerLine + "\t" + sample

outfile_drugs.write(drug_headerLine + "\n")
outfile_drugs_info.write(drug_info_headerLine + "\n")

sortedDrugDict = OrderedDict(sorted(dictDrugs.items(), key=lambda t: len(t[1]), reverse=True))

allDrugs = 0
for drug in sortedDrugDict.keys():
	allDrugs += 1
	drugInfos = sortedDrugDict[drug]
	sampleNum = len(drugInfos) - 2
	drugLine = "%s\t%s\t%s\t%s" %(drug,sampleNum,drugInfos[0],drugInfos[1])
	drug_info_line = "%s\t%s\t%s" %(drug,sampleNum,drugInfos[1])
	
	for sample in sampleNames:
		if sample in drugInfos:  # if sample names contain each other, this will create problems!
			drugLine = drugLine + "\t" + "1"
		else:
			drugLine = drugLine + "\t" + "0"
	
		if drug not in drugsToVariants.keys():
			print("Warning! Drug %s not contained in dictionary!" %(drug))
			continue
		if sample in drugsToVariants[drug].keys():
			targetedGenes = []
			for gene in drugsToVariants[drug][sample].keys():
				if drugsToVariants[drug][sample][gene][0] == 1:
					targetedGenes.append("%s (score: %s)" %(gene,drugsToVariants[drug][sample][gene][1]))
				else:
					targetedGenes.append("%s (multiple variants, score: %s)" %(gene,drugsToVariants[drug][sample][gene][1]))
			drug_info_line = drug_info_line + "\t" + " | ".join(targetedGenes)
		else:
			drug_info_line = drug_info_line + "\t" + "0"

	outfile_drugs.write(drugLine + "\n")
	outfile_drugs_info.write(drug_info_line + "\n")

outfile_samples.close()
outfile_drugs.close()
outfile_drugs_info.close()

print("Found %s different drugs in %s samples." %(allDrugs,len(sampleNames)))

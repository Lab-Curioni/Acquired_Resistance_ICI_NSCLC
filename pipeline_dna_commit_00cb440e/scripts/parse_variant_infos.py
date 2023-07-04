#!/usr/bin/env python

'''
seiler_tcga_2017
get variant information for oncoprint heatmaps
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


def parseInputFile(sample,sampleName,dictImpact,dictGenes):
	
	print("Sample %s. File: %s" %(sampleName,sample))

	infile = open(sample,'r')
	header = infile.readline()
	(index_genes,index_impact) = getColumnIndex(header)

	impactTooLow = 0

	for line in infile:
		lineSplit = line.strip().split("\t")
		genes = lineSplit[index_genes].split(";")
		impacts = lineSplit[index_impact].split(";")

		worstImpact = 1000    # Initialize
		worstImpact_gene = "" # If several genes have the same ranking, take the first annotated one
		worstImpact_name = ""

		for gene in genes:

			tempImpact = 1000
			temp_impact_name = ""

			for impact in impacts:
				impactSplit = impact.split(":")
				impactGene = impactSplit[0]  # mouse-based

				if impactGene != gene:  # the gene annotation currently looked at belongs to another gene
					continue
				
				impactAnno = impactSplit[1]
				if "&" in impactAnno:  # e.g. splice_region_variant&intron_variant
					annoSplit = impactAnno.split("&")
					for anno in annoSplit:
						if anno not in dictImpact.keys():
							print("Warning! Impact %s not contained in impact dictionary!" %(anno))
							continue
						impact_rank = dictImpact[anno]
						if impact_rank < tempImpact:
							tempImpact = impact_rank
							temp_impact_name = anno
				else:
					if impactAnno not in dictImpact.keys():
						print("Warning! Impact %s not contained in impact dictionary!" %(impactAnno))
						continue
					impact_rank = dictImpact[impactAnno]
					if impact_rank < tempImpact:
						tempImpact = impact_rank
						temp_impact_name = impactAnno

			if tempImpact == 1000:
				print("Warning! No impact found for annotated gene. Please check impact string and genes for line %s." %(line))
				continue

			if tempImpact < worstImpact:
				worstImpact = tempImpact
				worstImpact_gene = gene
				worstImpact_name = temp_impact_name

		# for this variant line gene and impact have been extracted, now add to dict
		# also respect multiple variants per sample per gene

		# threshold
		if worstImpact >= int(args.impact_threshold):
			impactTooLow +=1
			continue

		if worstImpact_gene not in dictGenes.keys():
			dictGenes[worstImpact_gene] = {}
		
		if sampleName not in dictGenes[worstImpact_gene].keys():
			dictGenes[worstImpact_gene][sampleName] = worstImpact_name
			continue

		impactString = dictGenes[worstImpact_gene][sampleName]
		if ";" in impactString:  # already multiple annotations
			currentImpact_name = impactString.split(";")[0]
			if worstImpact < dictImpact[currentImpact_name]:  # change the worst impact annotation to current one
				dictGenes[worstImpact_gene][sampleName] = worstImpact_name + ";MULTIPLE"
		else:
			if worstImpact < dictImpact[impactString]:  # change the worst impact annotation to current one
				dictGenes[worstImpact_gene][sampleName] = worstImpact_name + ";MULTIPLE"
			else:
				dictGenes[worstImpact_gene][sampleName] = impactString + ";MULTIPLE"

	infile.close()

	print("Ignored %s variants because too low impact."%(impactTooLow))

	return (dictGenes)

def getColumnIndex(header):
	headerSplit = header.strip().split("\t")

	index_genes = -1
	index_impact = -1
	
	for pos in range(0,len(headerSplit)):
		if args.colName_genes == headerSplit[pos]:
			index_genes = pos
		if args.colName_impact == headerSplit[pos]:
			index_impact = pos

	if (index_genes == -1) or (index_impact == -1):
		print("Error! Could not match all input columns in header %s." %(header))
		sys.exit(1)
	
	return (index_genes,index_impact)

def readIn_impactDict(dictImpact):
	infile_impacts = open(args.impactDict,'r')

	infile_impacts.readline()  # header line

	for line in infile_impacts:
		lineSplit = line.strip().split("\t")
		impact = lineSplit[0]
		ordering = int(lineSplit[1])

		if impact not in dictImpact.keys():
			dictImpact[impact] = ordering
		else:
			print("Warning! Impact %s is already included in impact dictionary!" %(impact))

	infile_impacts.close()
	return dictImpact



parser = argparse.ArgumentParser(description='Get variant information from files of form [sample].clinicalAnnotation.txt.')
parser.add_argument('--inputDir', dest='inputDir', required=True, help='Input directory with files of format [sample].clinicalAnnotation.txt.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name of the output file.')
parser.add_argument('--colName_genes', dest='colName_genes', required=True, help='Column name of column containing gene names.')
parser.add_argument('--colName_impact', dest='colName_impact', required=True, help='Column name of column containing info on the variant impact.')
parser.add_argument('--fileEnding', dest='fileEnding', required=True, help='To get correct input files, specify the file ending of desired input files.')
parser.add_argument('--impactDict', dest='impactDict', required=True, help='File with priortization of variant impacts.')
parser.add_argument('--impact_threshold', dest='impact_threshold', required=True, help='Threshold for disregarding impacts.')

args = parser.parse_args()

print("Parameters:\nInputDir: %s\noutFile: %s\ncolName_genes: %s\ncolName_impact: %s\nfileEnding: %s\nImpact dictionary: %s\nImpact threshold: %s\n" %(args.inputDir,args.outFile,args.colName_genes,args.colName_impact,args.fileEnding,args.impactDict,args.impact_threshold))

dictImpact = {}
dictImpact = readIn_impactDict(dictImpact)

# gene information
dictGenes = {} # gene to No. samples, (isDamaging), sampleEntry

# sample information
sampleNames = [] # for sorting

for file in os.listdir(args.inputDir):
	sampleFile = "%s%s" %(args.inputDir,os.path.basename(file))
	if os.path.isfile(sampleFile) and sampleFile.endswith(args.fileEnding):
		name = os.path.basename(file).split("_")[0]
		if name not in sampleNames:
			sampleNames.append(name)
		else:
			print("Error! Multiple files for %s!" %(name))
		(dictGenes) = parseInputFile(sampleFile,name,dictImpact,dictGenes)

sampleNames = natsorted(sampleNames)
outfileHeaderString = "\t".join(sampleNames)
outfile = open(args.outFile,'w')
outfile.write("Gene\tNum_mutated_samples\t" + outfileHeaderString + "\n")

sortedGeneDict = OrderedDict(sorted(dictGenes.items(), key=lambda t: len(t[1]), reverse=True))

allGenes = 0
for gene in sortedGeneDict.keys():
	allGenes += 1
	geneInfos = sortedGeneDict[gene]
	numSamples = len(geneInfos)

	geneLine = "%s\t%s"%(gene,numSamples)
	
	for sample in sampleNames:
		if sample in geneInfos:
			geneLine = geneLine + "\t" + geneInfos[sample]
		else:
			geneLine = geneLine + "\t" + "0"

	outfile.write(geneLine + "\n")

outfile.close()

print("Found %s different genes in %s samples." %(allGenes,len(sampleNames)))

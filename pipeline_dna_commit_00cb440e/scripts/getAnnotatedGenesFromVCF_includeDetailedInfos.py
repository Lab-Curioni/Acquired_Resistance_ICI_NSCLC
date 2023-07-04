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
	print("Get annotated genes from a vcf file.")
	print("Usage: python getAnnotatedGenesFromVCF_includeDetailedInfos.py [vcfFile] [outfile] [pathwayDB] [tumorName]")
	sys.exit(1)
if "-h" in sys.argv[1]:
	print("Get annotated genes from a vcf file.")
	print("Usage: python getAnnotatedGenesFromVCF_includeDetailedInfos.py [vcfFile] [outfile] [pathwayDB] [tumorName]")
	sys.exit(1)

vcfFile = sys.argv[1]
outfileName = sys.argv[2]
pathwayDB = sys.argv[3]
tumorName = sys.argv[4]

infilePathways = open(pathwayDB,'r')
dictPathways = {}

for linePW in infilePathways:
	lineSplit = linePW.strip().split()
	geneName = lineSplit[0].strip("\"")  # strip to remove the " symbols
	allPWs = "".join(lineSplit[1:])
	dictPathways[geneName] = allPWs.strip("\"")

infilePathways.close()


infile = open(vcfFile,'r')
outfile = open(outfileName,'w')
outfileCbioportal = open(outfileName + "_listForCbioportalQuery.txt",'w')
outfileCbioportal.write("Gene\tPosition\n")

outfile.write("Chromosome\tPosition\tChange\tFrequency\tGenes\tAnnotation\tVariantID\tFunctionalAnnotation\tPathway\tVariantImpact\tClinicalSignificance\n")
geneArr = []
allVariants = 0

pathwayInfo = 0
seenGenes = []  # only to provide correct number of genes for which we found  pathway

tumorIndex = 10
for line in infile:
	if line.startswith("##"):
		continue
	if line.startswith("#CHROM"):
		#header line, get the correct column of the tumor sample
		lineSplit = line.strip().split("\t")
		for i in range(0,len(lineSplit)):
			if lineSplit[i] == tumorName:
				tumorIndex = i
				print(tumorIndex)
				break
		continue

	lineSplit = line.strip().split("\t")
    # GATK CombineVariants nMin=2 also keeps variants called only by one caller if another variant was called by the same or another caller at the same site ... exclude those variants
	caller_set = lineSplit[7].split(";set=")[1]
	if (caller_set != "Intersection") and ("-" not in caller_set):
		continue

	allVariants += 1
	if "ANN=" in line:
		annoTag = lineSplit[7].split("ANN=")[1].split(";")[0]
		geneName = annoTag.split("|")[3]
		geneAnno = annoTag.split("|")[9] + "|" + annoTag.split("|")[10]
		mutationImpact = annoTag.split("|")[1]
		geneNameTempArr = [geneName]
		geneAnnoTempArr = [geneName+":"+geneAnno]
		mutaImpactTempArr = [geneName+":"+mutationImpact]

		dictTemp_genes = {}
		dictTemp_genes[geneName] = [geneName+":"+geneAnno]
		annoTagArr = annoTag.split(",")

		for i in range(0,len(annoTagArr)):
			annoTemp = annoTagArr[i]
			geneNameTemp = annoTemp.split("|")[3]
			mutationImpactTemp = annoTemp.split("|")[1]
			geneAnnoTemp = annoTemp.split("|")[9] + "|" + annoTemp.split("|")[10]
			if geneNameTemp not in geneNameTempArr:
				geneNameTempArr.append(geneNameTemp)
				geneAnnoTempArr.append(geneNameTemp+":"+geneAnnoTemp)
				mutaImpactTempArr.append(geneNameTemp+":"+mutationImpactTemp)
			if geneNameTemp not in dictTemp_genes:
				dictTemp_genes[geneNameTemp] = [geneNameTemp+":"+geneAnnoTemp]
			else:
				# in case of multiple transcripts, check whether a new annotation was found
				foundAnnoTemp = False
				for annoEntry in dictTemp_genes[geneNameTemp]:
					geneAnno_and_gene = geneNameTemp+":"+geneAnnoTemp
					if annoEntry == geneAnno_and_gene:
						foundAnnoTemp = True
						break
				if not foundAnnoTemp:
					# add this annotation to the annotation column
					geneAnnoTempArr.append(geneNameTemp+":"+geneAnnoTemp)
					mutaImpactTempArr.append(geneNameTemp+":"+mutationImpactTemp)
					dictTemp_genes[geneNameTemp].append(geneNameTemp+":"+geneAnnoTemp)
			if geneNameTemp not in geneArr:
				geneArr.append(geneNameTemp)
		allTempGeneNames = ";".join(geneNameTempArr)
		allTempGeneAnnos = ";".join(geneAnnoTempArr)
		allTempImpacts = ";".join(mutaImpactTempArr)

		### Reads out the CLNSIG_Tag for the corresponding mutation.
		### CLNVAR reports different mutations (e.g C>T, C>A) at the same position separated by a comma. 5|5|5,2|2
		### If there are multiple conditions reported for a particular mutations, CLNSIG is reported for each separated by a pipe.
		CLNSIG_Tag = "." # NOTE: 2018-01-18(FS) bugfix, because previously CLNSIG_Tag has not been initialized
		if "CLNSIG=" in line:
			CLNSIG_Tag = lineSplit[7].split("CLNSIG=")[1].split(";")[0]
			if len(CLNSIG_Tag)>1 and not CLNSIG_Tag=="255":
				Change = lineSplit[3]+">"+lineSplit[4]
				CLNHGVS_Tag = lineSplit[7].split("CLNHGVS=")[1].split(";")[0]
				Block = CLNHGVS_Tag.split(",")
				for i in range(1,len(Block)):
					if Change in Block[i]:
						CLNSIG_Tag=CLNSIG_Tag.split(",")[i]
						print(CLNSIG_Tag)
		#else: NOTE: 2018-01-17 (FS): this overrides the actual CLNSIG_Tag and thus would lead to a crash in case of multiple CLNHGVS_Tag blocks 
		#	CLNSIG_Tag = "."


		# note: currently the frequency extraction is adapted to the combination of varscan, mutect, and a variable other caller! TODO!
		frequency = lineSplit[tumorIndex].split(":")[-2]
		if "%" in frequency:
			frequency = frequency.split("%")[0]
		else:
			frequency = float(frequency) * 100.0

		functionalAnnotate = "."
		if "dbNSFP" in line:  # extract functional annotation
			functionalAnnotate = "dbNSFP_CADD" + lineSplit[7].split(";set=")[0].split("dbNSFP_CADD")[1]  # add dbNSFP_CADD to use dbNSFP_CADD for the split but include the annotation

		pathway = ""

		for tempGene in geneNameTempArr:
			outfileCbioportal.write(tempGene + "\t" + lineSplit[1] + "\n")
			if tempGene in dictPathways.keys():
				if tempGene not in seenGenes:
					pathwayInfo += 1
					seenGenes.append(tempGene)
				pathway += tempGene + ":" + dictPathways[tempGene] + ";"
		if len(pathway) == 0: # found none of the genes in pathway file
			pathway = "."

		outfile.write(lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[3] + ">" + lineSplit[4] + "\t" + str(frequency) + "\t" + allTempGeneNames + "\t" + allTempGeneAnnos + "\t" + lineSplit[2] + "\t" + functionalAnnotate + "\t" + pathway + "\t" + allTempImpacts + "\t" + CLNSIG_Tag + "\n")


	else:
		print("Error! no annotation tag in line %s." %(line))


infile.close()
outfile.close()
outfileCbioportal.close()

outfile1 = open(outfileName + "_distinctGeneNames.txt",'w')
outfile1.write("gene_names\n")
outString = "\n".join(geneArr)
outfile1.write(outString)
outfile1.close()
print("Extracted %s genes associated to %s variants, found pathway information for %s genes." %(len(geneArr),allVariants,pathwayInfo))

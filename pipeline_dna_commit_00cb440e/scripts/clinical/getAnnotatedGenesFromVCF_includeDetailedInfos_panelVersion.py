#!/usr/bin/env python

'''
Get annotated genes from a vcf file
assumes annotation at least with snpEff (searches the ANN tag in info column)
note: frequency format assumes tumor only panels, converted with script "convert_panelFile_to_vcfFile.py"
Unlike original getAnnotatedGenesFromVCF script, here the tumor name does not need to be provided, as only one column "SAMPLE" is included in the IonTorrent-based VCFs.
It might be necessary to change this for other panels.
Franziska Singer, August 2015
'''

import sys
import numpy as np

if len(sys.argv) <= 1:
	print "Get annotated genes from a vcf file."
	print "Usage: python getAnnotatedGenesFromVCF_includeDetailedInfos_panelVersion.py [vcfFile] [outfile] [pathwayDB]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Get annotated genes from a vcf file."
	print "Usage: python getAnnotatedGenesFromVCF_includeDetailedInfos_panelVersion.py [vcfFile] [outfile] [pathwayDB]"
	sys.exit(1)

vcfFile = sys.argv[1]
outfileName = sys.argv[2]
pathwayDB = sys.argv[3]

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

#outfile.write("Chromosome\tPosition\tChange\tFrequency\tGenes\tAnnotation\tVariantID\tFunctionalAnnotation\tPathway\n")
outfile.write("Chromosome\tPosition\tChange\tFrequency\tGenes\tAnnotation\tVariantID\tFunctionalAnnotation\tPathway\tVariantImpact\tClinicalSignificance\n")
geneArr = []
allVariants = 0

pathwayInfo = 0
seenGenes = []  # only to provide correct number of genes for which we found a pathway

for line in infile:
	if line.startswith("#"):
		continue

	lineSplit = line.strip().split("\t")
    
	allVariants += 1
	if "ANN=" in line:
		annoTag = lineSplit[7].split("ANN=")[1].split(";")[0]	
		geneName = annoTag.split("|")[3]
		geneAnno = annoTag.split("|")[9] + "|" + annoTag.split("|")[10]
		mutationImpact = annoTag.split("|")[1]
		geneNameTempArr = [geneName]
		geneAnnoTempArr = [geneAnno]
		mutaImpactTempArr = [mutationImpact]
		
		dictTemp_genes = {}
		dictTemp_genes[geneName] = [geneAnno]
		annoTagArr = annoTag.split(",")
		
		for i in range(0,len(annoTagArr)):
			annoTemp = annoTagArr[i]
			geneNameTemp = annoTemp.split("|")[3]
			mutationImpactTemp = annoTemp.split("|")[1]
			geneAnnoTemp = annoTemp.split("|")[9] + "|" + annoTemp.split("|")[10]
			if geneNameTemp not in geneNameTempArr:
				geneNameTempArr.append(geneNameTemp)
				geneAnnoTempArr.append(geneAnnoTemp)
				mutaImpactTempArr.append(mutationImpactTemp)
			if geneNameTemp not in dictTemp_genes:
				dictTemp_genes[geneNameTemp] = [geneAnnoTemp]
			else:
				# in case of multiple transcripts, check whether a new annotation was found
				geneAnnoTemp = annoTemp.split("|")[9] + "|" + annoTemp.split("|")[10]
				foundAnnoTemp = False
				for annoEntry in dictTemp_genes[geneNameTemp]:
					if annoEntry == geneAnnoTemp:
						foundAnnoTemp = True
						break
				if not foundAnnoTemp:
					# add this annotation to the annotation column
					geneAnnoTempArr.append(geneAnnoTemp)
					mutaImpactTempArr.append(mutationImpactTemp)
					dictTemp_genes[geneNameTemp].append(geneAnnoTemp)
				
			if geneNameTemp not in geneArr:
				geneArr.append(geneNameTemp)
		allTempGeneNames = ";".join(geneNameTempArr)
		allTempGeneAnnos = ";".join(geneAnnoTempArr)
		allTempImpacts = ";".join(mutaImpactTempArr)
		
		### Reads out the CLNSIG_Tag for the corresponding mutation.
		### CLNVAR reports different mutations (e.g C>T, C>A) at the same position separated by a comma. 5|5|5,2|2
		### If there are multiple conditions reported for a particular mutations, CLNSIG is reported for each separated by a pipe.
		CLNSIG_Tag = "."
		if "CLNSIG" in line:
			CLNSIG_Tag = lineSplit[7].split("CLNSIG=")[1].split(";")[0]
			if len(CLNSIG_Tag)>1 and not CLNSIG_Tag=="255":
				Change = lineSplit[3]+">"+lineSplit[4]
				CLNHGVS_Tag = lineSplit[7].split("CLNHGVS=")[1].split(";")[0]
				Block = CLNHGVS_Tag.split(",")
				for i in range(1,len(Block)):
					if Change in Block[i]:
						print Block[i]
						print CLNSIG_Tag
						CLNSIG_Tag=CLNSIG_Tag.split(",")[i]
						print(CLNSIG_Tag)
					#else:
					#	CLNSIG_Tag = "."

		# extract frequency, format: COV:AC:FREQ:PVAL:DP4 or GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
		formatTagArr = lineSplit[8].split(":")
		freqIndex = -1
		for i in range (0,len(formatTagArr)):
			if formatTagArr[i] == "FREQ":
				freqIndex = i
				break
		if freqIndex == -1:
			print "Error! Found no frequency tag in format column. Expected: COV:AC:FREQ:PVAL:DP4 or GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR"
			continue
		
		frequency = lineSplit[9].split(":")[freqIndex]  
		if "%" in frequency:
			frequency = frequency.split("%")[0]
			
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
			
		#outfile.write(lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[3] + ">" + lineSplit[4] + "\t" + str(frequency) + "\t" + allTempGeneNames + "\t" + allTempGeneAnnos + "\t" + lineSplit[2] + "\t" + functionalAnnotate + "\t" + pathway + "\n")
		outfile.write(lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[3] + ">" + lineSplit[4] + "\t" + str(frequency) + "\t" + allTempGeneNames + "\t" + allTempGeneAnnos + "\t" + lineSplit[2] + "\t" + functionalAnnotate + "\t" + pathway + "\t" + allTempImpacts + "\t" + CLNSIG_Tag + "\n")
				
	else:
		print "Error! no annotation tag in line %s." %(line)
		
		
infile.close()
outfile.close()
outfileCbioportal.close()

outfile1 = open(outfileName + "_distinctGeneNames.txt",'w')
outString = "\n".join(geneArr)
outfile1.write(outString)
outfile1.close()
print "Extracted %s genes associated to %s variants, found pathway information for %s genes." %(len(geneArr),allVariants,pathwayInfo)

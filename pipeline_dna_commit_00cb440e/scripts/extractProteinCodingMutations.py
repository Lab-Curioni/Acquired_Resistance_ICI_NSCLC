#!/usr/bin/env python

'''
Extract all protein coding, nonsynonmous, and damaging mutations from a given mutation table
Franziska Singer, June 2016
'''

import sys
import subprocess
import re

if len(sys.argv) <= 1:
	print("Extract protein coding mutations and compute mutational load.")
	print("Usage: python extractProteinCodingMutations.py [infile] [infileRegionLength] [outfile_mutationalLoad]")
	sys.exit(1)
if "-h" in sys.argv[1]:
	print("Extract protein coding mutations and compute mutational load.")
	print("Usage: python extractProteinCodingMutations.py [infile] [infileRegionLength] [outfile_mutationalLoad]")
	sys.exit(1)

infile = sys.argv[1]
infileRegionLength = sys.argv[2]
outfileNameTag = sys.argv[3]
outfileMutationalLoad = sys.argv[4]

allProteinCoding = 0
allNonsyno = 0
allDamaging = 0
intronVariants = 0
allVariants = 0

inF = open(infile,"r")
outfile = open(outfileNameTag+".nonsynonymous.txt","w")
outfileDam = open(outfileNameTag + ".damaging.txt",'w')
outfileIntron = open(outfileNameTag + ".intronVariants.txt",'w')

headerline = inF.readline()
outfile.write(headerline.strip() + "\n")
outfileDam.write(headerline.strip() + "\n")
outfileIntron.write(headerline.strip() + "\n")
headerSplit = headerline.strip().split("\t")

indexVI = -1
indexAN = -1
# extract the relevant column indices
for posTemp in range(0,len(headerSplit)):
	if "VariantImpact" in headerSplit[posTemp]:
		indexVI = posTemp
	if "Annotation" == headerSplit[posTemp]:
		indexAN = posTemp

if (indexVI == -1):
	print("Error! Did not find column VariantImpact!")
	sys.exit(0)
if (indexAN == -1):
	print("Error! Did not find column Annotation!")
	sys.exit(0)

for line in inF:
	lineSplit = line.strip().split("\t")
	
	if ("intron_variant" in lineSplit[indexVI]) or ("splice" in lineSplit[indexVI]):   # currently supports snpEff annotations
		outfileIntron.write(line.strip() + "\n")
		intronVariants += 1
	
	tempSplit = lineSplit[indexAN].split("|p.")
	allVariants += 1
	isSynonymous = False
	if len(tempSplit) > 1:
		allProteinCoding += 1
		#tempSplit = lineSplit[5].split("|p.")
		myString = tempSplit[1].split("|")[0].split(";")[0]  # split twice in case there are several annotations
		myNonnumberString = re.sub("[0-9]", ",", myString)
		finalSplit = myNonnumberString.split(",")
		if finalSplit[0] == finalSplit[len(finalSplit)-1]:
			# synonymous mutation
			isSynonymous = True
			continue
		else:
			allNonsyno += 1
			outfile.write(line.strip() + "\n")
	# check if annotated as damaging, include all databases
	if (isSynonymous == False) and ("_pred=D" in line):
    #if (isSynonymous == False) and ("_pred=D" in lineSplit[8]):
		allDamaging += 1
		outfileDam.write(line.strip() + "\n")
		
print("Parsed %s variants. Found %s protein coding mutations. Of them %s are nonsynonymous. Found %s damaging mutations. Found %s intron variants." %(allVariants,allProteinCoding,allNonsyno,allDamaging,intronVariants))

infileRL = open(infileRegionLength,'r')
regionLength = float(infileRL.readline().strip())
infileRL.close()

mutLoad = (float(allNonsyno)/regionLength) * (1000000.0)
outfileMutBurden = open(outfileMutationalLoad,'w')
outfileMutBurden.write("%s"%(mutLoad))
outfileMutBurden.close()

outfileStats = open(outfileNameTag+".stats.txt",'w')
outfileStats.write("all_variants\tprotein_coding\tnon_synonymous_protein_coding\tdamaging\tintronic\n%s\t%s\t%s\t%s\t%s\n"%(allVariants,allProteinCoding,allNonsyno,allDamaging,intronVariants))
outfileStats.close()

inF.close()
outfile.close()
outfileDam.close()
outfileIntron.close()

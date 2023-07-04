#!/usr/bin/env python
"""

Convert a panel "vcf" file to a vcf file ready for annotation and other post-processing steps
Franziska Singer, August 2016

"""

usage = """
%prog [args]

Help on python script convert_panelFile_to_vcfFile.py

Script to convert the table like variant output format of different 
panel analysis methods to vcf file format that can be post-processed e.g. with 
annotation methods etc. 

Additional panel types can be added as separate arguments, guiding which coverting method is chosen.

"""


import sys
import argparse
import os
from natsort import natsort
import re

''' 
--------------------------------------------------------------
		method definitions 
--------------------------------------------------------------
'''

'''
extract the indices of all interesting columns in variant table for DNA-only panels from Basel (tumor only)
stores each index in a dictionary
Additionally, extract format description and info tag dexcription for vcf header
format line of form: ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
info line of form: ##INFO=<ID=AA,Number=1,Type=String,Description="Peptide annotation">
'''
def getColumnIndices_Basel(lineSplit):
	dict_featureIndices = {}
	formatLines = ""
	formatLines += "##FORMAT=<ID=FREQ,Number=.,Type=String,Description=\"Frequency of variant allele\">\n"
	infoLines = ""
	for i in range(0,len(lineSplit)):
		if "Locus" == lineSplit[i]:
			dict_featureIndices["ChromPos"] = i
			continue
		if "Genotype" == lineSplit[i]:
			dict_featureIndices["ALT"] = i
			continue
		if "Ref" == lineSplit[i]:
			dict_featureIndices["REF"] = i
			continue
		if "Genes" == lineSplit[i]:
			dict_featureIndices["Genes"] = i
			infoLines += "##INFO=<ID=Gene,Number=.,Type=String,Description=\"Gene annotation\">\n" 
			continue
		if "Variant ID" == lineSplit[i]:
			dict_featureIndices["VariantID"] = i
			continue
		if "Transcript" == lineSplit[i]:
			dict_featureIndices["Transcript"] = i
			infoLines += "##INFO=<ID=Transcript,Number=.,Type=String,Description=\"Transcript of gene annotation\">\n"
			continue
		if "Coding" == lineSplit[i]:
			dict_featureIndices["Coding"] = i
			infoLines += "##INFO=<ID=Coding,Number=.,Type=String,Description=\"Change on genomic level\">\n"
			continue
		if "Amino Acid Change" == lineSplit[i]:
			dict_featureIndices["AAchange"] = i
			infoLines += "##INFO=<ID=AA,Number=.,Type=String,Description=\"Amino acid change\">\n"
			continue
		if "Variant Effect" == lineSplit[i]:
			dict_featureIndices["VarEff"] = i
			infoLines += "##INFO=<ID=VarEff,Number=.,Type=String,Description=\"Variant effect\">\n"
			continue
		if "SIFT" == lineSplit[i]:
			dict_featureIndices["SIFT"] = i
			infoLines += "##INFO=<ID=SIFT,Number=.,Type=float,Description=\"SIFT score\">\n"
			continue
		if "PolyPhen" == lineSplit[i]:
			dict_featureIndices["PolyPhen"] = i
			infoLines += "##INFO=<ID=PolyPhen,Number=.,Type=float,Description=\"PolyPhen score\">\n"
			continue
		if "PFAM" == lineSplit[i]:
			dict_featureIndices["PfamDomain"] = i
			infoLines += "##INFO=<ID=Pfam,Number=.,Type=String,Description=\"Domain annotation from Pfam\">\n"
			continue
		if "dbSNP" == lineSplit[i]:
			dict_featureIndices["dbSNP"] = i
			continue
		if "MAF" == lineSplit[i]:
			dict_featureIndices["MAF"] = i
			infoLines += "##INFO=<ID=MAF,Number=.,Type=String,Description=\"MAF score\">\n"
			continue
		if "UCSC Common SNPs" == lineSplit[i]:
			dict_featureIndices["COMMON"] = i
			infoLines += "##INFO=<ID=UCSCcommon,Number=.,Type=String,Description=\"Annotated as common variant in UCSC\">\n"
			continue
		if "COSMIC" == lineSplit[i]:
			dict_featureIndices["CosmicEntry"] = i
			infoLines += "##INFO=<ID=CosmicDescr,Number=.,Type=String,Description=\"Disease description from COSMIC\">\n"
			continue
		if "OMIM" == lineSplit[i]:
			dict_featureIndices["GeneName"] = i
			infoLines += "##INFO=<ID=OMIM,Number=.,Type=String,Description=\"Full gene name description\">\n"
			continue
		if "DrugBank" == lineSplit[i]:
			dict_featureIndices["DrugBank"] = i
			infoLines += "##INFO=<ID=DrugBank,Number=.,Type=String,Description=\"Entries in DrugBank associated to variant\">\n"
			continue
		if "ClinVar" == lineSplit[i]:
			dict_featureIndices["ClinVar"] = i
			infoLines += "##INFO=<ID=ClinVar,Number=.,Type=String,Description=\"ClinVar entry associated to variant\">\n"
			continue
		if "Allele Coverage" == lineSplit[i]:
			dict_featureIndices["AlleleCov"] = i
			formatLines += "##FORMAT=<ID=AC,Number=.,Type=String,Description=\"Allele coverage of observed alleles\">\n"
			continue
		if "p-value" == lineSplit[i]:
			dict_featureIndices["pVal"] = i
			formatLines += "##FORMAT=<ID=PVAL,Number=.,Type=String,Description=\"Pvalue associated to variant\">\n"
			continue
		if "Coverage" == lineSplit[i]:
			dict_featureIndices["Coverage"] = i
			formatLines += "##FORMAT=<ID=COV,Number=.,Type=String,Description=\"Nucleotide coverage at variant position\">\n"
			continue
		if "Ref+/Ref-/Var+/Var-" == lineSplit[i]:
			dict_featureIndices["DP4"] = i
			formatLines += "##FORMAT=<ID=DP4,Number=.,Type=String,Description=\"Coverage of reference forward, reference reverse, alternative forward, alternative reverse\">\n"
			continue
		if "Homopolymer Length" == lineSplit[i]:
			dict_featureIndices["HPlength"] = i
			infoLines += "##INFO=<ID=HPlength,Number=.,Type=String,Description=\"Length of neighboring homopolymers\">\n"
			continue
	return (dict_featureIndices,formatLines,infoLines)
		

'''
fill the INFO column with information from DNA-only panels from Basel (tumor only)
'''
def fillInfoString_Basel(lineSplit,dict_featureIndices):
	infoString = ""
	infoString += "Gene=%s;" %(lineSplit[dict_featureIndices["Genes"]] if len(lineSplit[dict_featureIndices["Genes"]]) > 0 else "N/A")
	infoString += "Transcript=%s;" %(lineSplit[dict_featureIndices["Transcript"]] if len(lineSplit[dict_featureIndices["Transcript"]]) > 0 else "N/A")
	infoString += "Coding=%s;" %(lineSplit[dict_featureIndices["Coding"]] if len(lineSplit[dict_featureIndices["Coding"]]) > 0 else "N/A")
	infoString += "AA=%s;" %(lineSplit[dict_featureIndices["AAchange"]] if len(lineSplit[dict_featureIndices["AAchange"]]) > 0 else "N/A")
	infoString += "VarEff=%s;" %(lineSplit[dict_featureIndices["VarEff"]] if len(lineSplit[dict_featureIndices["VarEff"]]) > 0 else "N/A")
	infoString += "SIFT=%s;" %(lineSplit[dict_featureIndices["SIFT"]] if len(lineSplit[dict_featureIndices["SIFT"]]) > 0 else "N/A")
	infoString += "PolyPhen=%s;" %(lineSplit[dict_featureIndices["PolyPhen"]] if len(lineSplit[dict_featureIndices["PolyPhen"]]) > 0 else "N/A")
	infoString += "Pfam=%s;" %(lineSplit[dict_featureIndices["PfamDomain"]] if len(lineSplit[dict_featureIndices["PfamDomain"]]) > 0 else "N/A")
	infoString += "MAF=%s;" %(lineSplit[dict_featureIndices["MAF"]] if len(lineSplit[dict_featureIndices["MAF"]]) > 0 else "N/A")
	infoString += "UCSCcommon=%s;" %(lineSplit[dict_featureIndices["COMMON"]] if len(lineSplit[dict_featureIndices["COMMON"]]) > 0 else "N/A")
	infoString += "CosmicDescr=%s;" %(lineSplit[dict_featureIndices["CosmicEntry"]] if len(lineSplit[dict_featureIndices["CosmicEntry"]]) > 0 else "N/A")
	infoString += "OMIM=%s;" %(lineSplit[dict_featureIndices["GeneName"]] if len(lineSplit[dict_featureIndices["GeneName"]]) > 0 else "N/A")
	infoString += "DrugBank=%s;" %(lineSplit[dict_featureIndices["DrugBank"]] if len(lineSplit[dict_featureIndices["DrugBank"]]) > 0 else "N/A")
	infoString += "ClinVar=%s;" %(lineSplit[dict_featureIndices["ClinVar"]] if len(lineSplit[dict_featureIndices["ClinVar"]]) > 0 else "N/A")
	infoString += "HPlength=%s;" %(lineSplit[dict_featureIndices["HPlength"]] if len(lineSplit[dict_featureIndices["HPlength"]]) > 0 else "N/A")
	
	return infoString

'''
given the allele coverage column, use the genotype column to extract the correct reference and variant allele
and calculate frequency; also calculate dp4
'''
def resolveAlleleCoverages_Basel(lineSplit,dict_featureIndices):
	genotype = lineSplit[dict_featureIndices["ALT"]]
	refString = genotype.split("/")[0]
	altString = genotype.split("/")[1]
	
	coverage = float(lineSplit[dict_featureIndices["Coverage"]])
	if (coverage == 0):
		print "Warning! Found zero nucleotide coverage for entry: %s in locus: %s" %(lineSplit[dict_featureIndices["Coverage"]],lineSplit[dict_featureIndices["ChromPos"]])
		
	refCov = 0
	altCov = 0
	
	alleleCovs = lineSplit[dict_featureIndices["AlleleCov"]].split(",")  # format example: G=268, A=1732, T=0
	
	for allele in alleleCovs:
		alleleSplit = allele.strip().split("=")
		
		if refString == alleleSplit[0].strip():  # found reference coverage
			refCov = int(alleleSplit[1])
			continue
		if altString == alleleSplit[0].strip():	# found alternative coverage 
			altCov = int(alleleSplit[1])
			continue
	if (refCov == 0) or (altCov == 0):
		print "Warning! Found zero allele coverage for entry: %s in locus: %s" %(lineSplit[dict_featureIndices["AlleleCov"]],lineSplit[dict_featureIndices["ChromPos"]])
		
	acString = str(refCov) + "," + str(altCov)
	frequency = float(float(altCov)/coverage) * 100.0
	
	dp4Column = lineSplit[dict_featureIndices["DP4"]].split(",")
	
	refFo = ""
	refRev = ""
	altFo = ""
	altRev = ""
	
	for dp in dp4Column:
		dpSplit = dp.strip().split("=")
		if refString == dpSplit[0].strip():  # found reference entry
			refFo = dpSplit[1].split("/")[0]
			refRev = dpSplit[1].split("/")[1]
			continue
		if altString == dpSplit[0].strip():	# found alternative entry 
			altFo = dpSplit[1].split("/")[0]
			altRev = dpSplit[1].split("/")[1]
			continue
	
	dp4 = refFo + "," + refRev + "," + altFo + "," + altRev
	freqString = "%.2f%%" %(frequency)
	return (acString,freqString,dp4)
	
	
	
'''
Converts variant tables of format type of IonReporter result based on DNA-only panels from Basel (tumor only).
Columns in table:
Locus	Genotype	Ref	Type	Genes	Location	Length	Info	Variant ID	Variant Name	
Strand	Transcript	Coding	Amino Acid Change	Variant Effect	PhyloP	SIFT	Grantham	
PolyPhen	PFAM	dbSNP	DGV	MAF	EMAF	AMAF	GMAF	UCSC Common SNPs	COSMIC	OMIM	
Gene Ontology	DrugBank	ClinVar	Allele Coverage	Allele Ratio	p-value	Coverage	Ref+/Ref-/Var+/Var-	Homopolymer Length 
'''
def convertToVCF_Basel(args):
	infile = open(args.infilePanel,'r')
	outfile = open(args.outfileVCF,'w')
	
	outfile.write("##fileformat=VCFv4.2\n")
	
	dict_chromVariants = {}  # maps each chromosome to variant positions to the actual description. Necessary because tables are not sorted.
	counter_variants = 0
	
	for line in infile:
		if line.startswith("##"):  # simply add comment lines to file
			outfile.write(line.strip() + "\n")
			continue
		
		lineSplit = line.strip().split("\t")
		if line.startswith("Locus"): # get the index of all interesting columns, additionally write correct header line and info lines
			(dict_featureIndices,formatLines,infoLines) = getColumnIndices_Basel(lineSplit)
			outfile.write(infoLines)
			outfile.write(formatLines)
			outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")				
			continue
				
		chromPos = lineSplit[dict_featureIndices["ChromPos"]]
		
		if len(chromPos) == 0:
			continue
		
		counter_variants += 1
		chromosome = chromPos.split(":")[0]
		position = chromPos.split(":")[1]
		
		if chromosome not in dict_chromVariants:
			dict_chromVariants[chromosome] = {}
		if position not in dict_chromVariants[chromosome]:
			dict_chromVariants[chromosome][position] = ""
			
		variantString = chromosome + "\t" + str(position)  # add columns CHROM and POS
		
		idCosmic = lineSplit[dict_featureIndices["VariantID"]].strip("\"")
		idDBSNP = lineSplit[dict_featureIndices["dbSNP"]].strip("\"")
		
		ids = ""
		if len(idCosmic) > 0:
			ids += idCosmic
			if len(idDBSNP) > 0:
				ids += ";" + idDBSNP
		else:
			if len(idDBSNP) > 0:
				ids += idDBSNP
		
		variantString += "\t" + ids + "\t" + lineSplit[dict_featureIndices["REF"]] + "\t" + lineSplit[dict_featureIndices["ALT"]].split("/")[1] # add columns ID, REF, and ALT	
		variantString += "\t.\tPASS"  # add columns QUAL and FILTER (here always same values)
		
		infoString = fillInfoString_Basel(lineSplit,dict_featureIndices)		
		formatString = "COV:AC:FREQ:PVAL:DP4"
		
		variantString += "\t" + infoString + "\t" + formatString # add info tag and format tag
		
		(acString,frequency,dp4) = resolveAlleleCoverages_Basel(lineSplit,dict_featureIndices)
		sampleString = lineSplit[dict_featureIndices["Coverage"]] + ":" + acString + ":" + frequency + ":" + re.sub(",", ".", lineSplit[dict_featureIndices["pVal"]]) + ":" + dp4
		
		variantString += "\t" + sampleString  # add the sample features
		
		# now everything is added to the variant line, so include in dictionary for later sorting	
		dict_chromVariants[chromosome][position] = variantString
		
	infile.close()
	
	print "Parsed %s variants." %(counter_variants)
	
	# sort dictionary
	
	for chrom in natsort(dict_chromVariants):
		for pos in natsort(dict_chromVariants[chrom]):
			outfile.write(dict_chromVariants[chrom][pos] + "\n")
			
	outfile.close()


''' 
--------------------------------------------------------------
		main program 
--------------------------------------------------------------
'''
	
# usage = """%prog [args]""" 

if __name__ == "__main__":
    
	parser = argparse.ArgumentParser()   # Default: %default , required=True
	parser.add_argument("-i", "--infile", dest="infilePanel", default=None, help="Specify the input panel variant file.",required = True)
	parser.add_argument("-o", "--outfile", dest="outfileVCF", default="panelVariants_convertedToVCF.vcf", help="Specify the output vcf file. Is overwritten if already existing.")
	parser.add_argument("-t", "--panelType", dest="panelType", default=None, help="Specify the type of the panel. Choose one of the following: Basel_IonTorrent_DNA.",required=True)
	parser.add_argument("-f", "--filterBasel", dest="filterBasel", action="store_true", default=False, help="If specified, apply homopolymer filter to Basel panel data.")
	#parser.add_argument("-t", "--threshold", dest="thresh", type=float, default=0.8, help="Threshold for required overlap of different predictions (percent of length). Default: 0.8.")
    
	args = parser.parse_args()


# first check if input file is actually there

if not os.path.exists(args.infilePanel):
	raise Exception("Could not find input file with name '%s'."%args.infilePanel)
			
# check which panel type and choose converting method accordingly	

if "Basel_IonTorrent_DNA" in args.panelType:
	convertToVCF_Basel(args)
	
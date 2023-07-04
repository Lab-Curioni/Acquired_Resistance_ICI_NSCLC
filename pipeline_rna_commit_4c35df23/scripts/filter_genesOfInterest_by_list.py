#!/usr/bin/env python

'''
Script for filtering the output table of the tcga cohort comparison RNA-Seq data, which is {sample}.final_rna_table.txt, for clinically interesting genes that are given to script in separate gene list

Anne Richter, Mai 2018
'''

import argparse

# retrieve command line arguments
parser = argparse.ArgumentParser(description='Script filters rows of tcga cohort comparison final output for clinically interesting genes')
parser.add_argument('--inFile', dest='inFile', required=True, help='Input file of type {sample}.final_rna_table.txt containing expression data of RNA-Seq experiment of patient compared to tcga cohort')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name/path of output file.')
parser.add_argument('--geneList', dest='geneList', required=True, help='List of genes that are clinically interesting for the given patient, tab delimited file with header and list of genes in first column')

args = parser.parse_args()

# dictionary will be filled with genes from clinical list as keys and rest of lines of table as values
dictGenes = {}

# get list of genes of interest into dictGenes
genelist = open(args.geneList, 'r')
headerLineGenes = genelist.readline()

for line in genelist:
    lineSplit = line.strip().split('\t')
    geneName = lineSplit[0]
    geneName = geneName.upper()
    dictGenes.setdefault(geneName,[])
genelist.close()

# get infile which is table with final output of tcga RNA-Seq cohort comparison
infile = open(args.inFile, 'r')
outfile = open(args.outFile, 'w')

# retrieve header from infile and write to outfile
headerLine = infile.readline()
outfile.write(headerLine)

# go through infile
for line in infile:
    lineSplit = line.strip().split('\t')
    geneID = lineSplit[0]
    geneID = geneID.upper()
# write rows of filtered genes to outfile
    if geneID in dictGenes:
        outfile.write(line)
#        print(line)
        continue

infile.close()
outfile.close()

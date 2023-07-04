#!/usr/bin/env python

'''
Given a tab separated output table including one column with gene IDs (hgnc symbols) and an output table of rDGIdb, this script adds one column to the tab separated table with the results of a rDGIdb query.

Anne Richter, June 2018
'''

import argparse

# import input table, header of gene column, dgidb result table and the path to the output file
parser = argparse.ArgumentParser(description='Script adds one column with interacting drugs to an existing table including one column with gene IDs (hgnc symbols), using the output of rDGIdb query')
parser.add_argument('--inFile', dest='inFile', required=True, help='Input table, tab delimited, with one column with hgnc symbols.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Path and name of output table')
parser.add_argument('--hgncCol', dest='hgncCol', required=True, help='Header of column that contains hgnc symbols')
parser.add_argument('--dgidbTable', dest='dgidbTable', required=True, help='Output table of the rDGIdb query')

args = parser.parse_args()

# dictionary with gene hgnc symbol as key and drugs as value
dictGenesDGI = {}

# open dgidb output table and fill dictionary with genes (keys) and drugs (values)
dgidbtable = open(args.dgidbTable, 'r')
headerDGI = dgidbtable.readline()

for line in dgidbtable:
    lineSplit = line.strip().split('\t')
    geneName = lineSplit[0]
    drugs = lineSplit[1]
    dictGenesDGI[geneName] = drugs

dgidbtable.close()

# open input table for reading and output table for writing
infile = open(args.inFile,'r')
outfile = open(args.outFile,'w')

# save header line of input table and find column with hgnc symbols
headerLine = infile.readline()
headerLine = headerLine.strip()
headerSplit = headerLine.split('\t')
index_hgncCol = None
for pos in range(0,len(headerSplit)):
    if args.hgncCol in headerSplit[pos]:
        index_hgncCol = pos

# write header to outfile
outHeader = headerLine + "\tdgidb_drugs"
outfile.write(outHeader + "\n")

# match each line to dgidb result and write line to outfile
for line in infile:
    lineStrip = line.strip()
    lineSplit = lineStrip.split('\t')
    hgncSymbol = lineSplit[index_hgncCol]
    print(hgncSymbol)
    if hgncSymbol in dictGenesDGI:
        outLine = lineStrip + "\t" + dictGenesDGI[hgncSymbol]
        outfile.write(outLine + "\n")
    else:
        outLine = lineStrip + "\t" + "NA"
        outfile.write(outLine + "\n")

outfile.close()

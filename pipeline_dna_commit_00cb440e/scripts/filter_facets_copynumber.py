#!/usr/bin/env python

'''
Filter copy number calls from facets
Franziska Singer, October 2018
'''

import sys
import numpy as np
import os
import re
import argparse
import math

'''
function definitions
'''

def getColumnIndex(firstInputLine):
    firstLineSplit = firstInputLine.split('\t')
    index_totalCopyCol = -1
    for pos in range(0,len(firstLineSplit)):
        if args.colName_totalCopy in firstLineSplit[pos]:
            index_totalCopyCol = pos
            break
    if index_totalCopyCol == -1:
        print("Error! Did not find a column matching %s in input header %s." %(args.colName_totalCopy,firstInputLine))
        sys.exit(1)

    return index_totalCopyCol

'''
main body
'''

parser = argparse.ArgumentParser(description='Filter facets copy number calls. Remove normal (here: total copy t.cn = 2) events and provide categorization for remaining CNVs (based on TCGA ABSOLUTE categories).')
parser.add_argument('--infile', dest='inFile', required=True, help='Input table with copy number calls, tab separated.')
parser.add_argument('--outfile', dest='outFile', required=True, help='Name of the filtered output file.')
parser.add_argument('--colName_totalCopy', dest='colName_totalCopy', required=True, help='Column name of column containing the total copy number state.')
#parser.add_argument('--colName_minorCopy', dest='colName_minorCopy', required=True, help='Column name of column containing the minor alelle copy state.')

args = parser.parse_args()

infile = open(args.inFile,'r')
outfile = open(args.outFile,'w')

firstInputLine = infile.readline().strip()
index_totalCopyCol = getColumnIndex(firstInputLine)

outHeader = firstInputLine
outHeader += "\tCNV_ratio\tCNV_log2ratio\tCNV_category"
outfile.write(outHeader + "\n")

num_all = 0
num_outfilter = 0
num_del = 0
num_loss = 0
num_gain = 0
num_amp = 0

for line in infile:
	lineSplit = line.strip().split("\t")
	num_all += 1
	totalCopy = float(lineSplit[index_totalCopyCol])
	
	if totalCopy == 2:
		num_outfilter += 1
		continue

	outfileLine = line.strip()

	ratio = float(totalCopy / 2.0)
	log2ratio = -10.0 # to make sure that copy number = 0 is counted as a deletion
	if ratio > 0.0:
		log2ratio = float(math.log(ratio,2))
	category = "N.A"
	if log2ratio > 1.0:
		category = "AMP"
		num_amp += 1
	elif log2ratio >= 0.58:
		category = "GAIN"
		num_gain += 1
	elif log2ratio < -1.0:
		category = "DEL"
		num_del += 1
	elif log2ratio < -0.42:
		category = "LOSS"
		num_loss += 1

	if category == "N.A":
		print("Warning! Could not categorize copy number for line %s." %(line))

	#print(num_all,num_outfilter,num_del,num_loss,num_gain,num_amp)
	outfileLine = "%s\t%s\t%.3f\t%s" %(outfileLine,ratio,log2ratio,category)
	outfile.write(outfileLine + "\n")

infile.close()
outfile.close()

print("Parsed events: %s\nFiltered out (neutral): %s\nThe remaining CNVs include:\nHomozygous deletions (DEL): %s\nCopynumber losses (LOSS): %s\nCopynumber gains (GAIN): %s\nHigh amplifications (AMP): %s" %(num_all,num_outfilter,num_del,num_loss,num_gain,num_amp))

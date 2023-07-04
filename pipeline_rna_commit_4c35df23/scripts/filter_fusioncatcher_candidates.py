#!/usr/bin/env python

'''
Script for filtering the output table of fusioncatcher, {sample}_final-list_candidate-fusion-genes.GRCh37.txt. All lines of the output table that have one of the keywords that are specified in the file description_tags_false_positive_fusioncatcher.txt in the description column, are filtered out.

Anne Richter, June 2018
'''


import argparse

# retrieve command line arguments
parser = argparse.ArgumentParser(description='Script filters rows of output file of fusioncatcher. Rows that contain at least one of the keywords specified in the file description_tags_false_positive_fusioncatcher.txt are discarded.')
parser.add_argument('--inFile', dest='inFile', required=True, help='Input file of type {sample}_final-list_candidate-fusion-genes.GRCh37.txt containing candidate fusion genes predicted by fusioncatcher.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name/path of filtered output file.')
parser.add_argument('--blacklist', dest='blackList', required=True, help='List of keywords that mark a candidate fusion gene to have a higher probability of being a false positive.')

args = parser.parse_args()

# list will be filled with keywords that mark a fusion candidate as having a high probability of being false positive
keywords = []

blacklist = open(args.blackList, 'r')

for line in blacklist:
    line = line.strip()
    keywords.append(line)

blacklist.close()

# get input file that contains fusion candidate genes
infile = open(args.inFile, 'r')
outfile = open(args.outFile, 'w')

# retrieve header from infile and write to outfile
headerLine = infile.readline()
outfile.write(headerLine)

# go through infile
for line in infile:
    lineSplit = line.strip().split('\t')
    descriptions = lineSplit[2].strip().split(',')
    keyword_count = 0
    for element in descriptions:
        if element in keywords:
            keyword_count = keyword_count + 1;
    if keyword_count == 0:
        outfile.write(line)

infile.close()
outfile.close()

#!/usr/bin/env python

'''
Adapt the names of given panel result files (format does not matter)
and write sample map for snakemake
Franziska Singer, June 2016
'''

import sys
import os

if len(sys.argv) <= 1:
	print "Prepare samples for clinical snakemake pipleine."
	print "Usage: python prepareMutationFilesForSnakemakePipeline.py [inDrectory] [outDirectory] [snakemakeSampleMapFile]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Prepare samples for clinical snakemake pipleine."
	print "Usage: python prepareMutationFilesForSnakemakePipeline.py [inDrectory] [outDirectory] [snakemakeSampleMapFile]"
	sys.exit(1)

inDirec = sys.argv[1]
outDirec = sys.argv[2]
snakeMap = sys.argv[3]

# files are written to new directory. In addition the sampleMap for snakemake is created

if not os.path.exists(outDirec):		# if the output directory does not exist, create it
    os.makedirs(outDirec)

outfile = open(snakeMap,'w')

sampleCounter = 1
for file in os.listdir(inDirec):
	sample = "%s%s" %(inDirec,os.path.basename(file))
	if os.path.isfile(sample):
		sampleName = os.path.basename(file).split(".")[0]
		
		outfile.write("S%s\t%s\tT\t1\nS%s\tnormal\tN\t1\n" %(sampleCounter, sampleName, sampleCounter))  # two lines, one for tumor (e.g. S1\tsampleName\tT\t1), one for normal (e.g. S1\tnormal\tN\t1)
		sampleCounter += 1
		
		cpAndRename = "cp %s %s/ ; mv %s/%s %s/%s_vs_normal.%s" %(sample, outDirec, outDirec, os.path.basename(file), outDirec, sampleName, os.path.basename(file).split(".")[1])
		os.system(cpAndRename)
		

outfile.close()
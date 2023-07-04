#!/usr/bin/env python
'''
reads out the amplicon length from a gtf file
@author: zickmannf
'''

import sys
import numpy

pathToBed = sys.argv[1]
pathToOutput = sys.argv[2]

infile = open(pathToBed, 'r')

exonLength = 0

for line in infile:
	if line.startswith("#") or line.startswith("browser") or line.startswith("track"):
		continue
	arr = line.split("\t")
	exonLength = exonLength + (int(arr[2]) - int(arr[1])) + 1

infile.close()

outfile = open(pathToOutput,'w')
outfile.write("%s"%(exonLength))
outfile.close()

print("Amplicon length = %s" %(exonLength))

#!/usr/bin/env python

'''
Create the readInput text file required by excavator 
Franziska Singer, May 2016
'''

import sys
import os

if len(sys.argv) <= 1:
	print "Create the readInput text file required by excavator."
	print "Usage: python createExcavatorReadInputfile.py [outfilenameReadInput] [experimentTarget] [assembly] [tumorBam] [normalBam]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Create the readInput text file required by excavator."
	print "Usage: python createExcavatorReadInputfile.py [outfilenameReadInput] [experimentTarget] [assembly] [tumorBam] [normalBam]"
	sys.exit(1)

outfilenameReadInput = sys.argv[1]
experimentTarget = sys.argv[2]
assembly = sys.argv[3]
tumorBam = sys.argv[4]
normalBam = sys.argv[5]

outfileReadInput = open(outfilenameReadInput,'w')
tumorName = os.path.basename(tumorBam).split(".bam")[0]
normalName = os.path.basename(normalBam).split(".bam")[0]
	
outfileReadInput.write("%s %s %s %s C1\n" %(experimentTarget,assembly,normalBam,normalName))
outfileReadInput.write("%s %s %s %s T1\n" %(experimentTarget,assembly,tumorBam,tumorName))
outfileReadInput.close()
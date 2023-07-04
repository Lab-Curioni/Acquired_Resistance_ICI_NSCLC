#!/usr/bin/env python

'''
Create the target text file required by excavator 
Franziska Singer, May 2016
'''

import sys
import os

if len(sys.argv) <= 1:
	print "Create the target file required by excavator."
	print "Usage: python createExcavatorTargetfile.py [outfileNameTarget] [wiggle] [reference]"
	sys.exit(1)
if "-h" in sys.argv[1]:
	print "Create the target file required by excavator."
	print "Usage: python createExcavatorTargetfile.py [outfileNameTarget] [wiggle] [reference]"
	sys.exit(1)

outfileNameTarget = sys.argv[1]
wiggle = sys.argv[2]
reference = sys.argv[3]

outfileTarget = open(outfileNameTarget,'w')

outfileTarget.write("%s %s\n" %(wiggle,reference))
outfileTarget.close()
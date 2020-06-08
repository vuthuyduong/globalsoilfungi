#!/usr/bin/env python
import sys
import os
import encodings

fastafilename = sys.argv[1] # soil.old.fasta
outputfilename = sys.argv[2] #soil.fasta

fastafile=open(fastafilename)
outputfile=open(outputfilename,"w")

isempty=False
for line in fastafile:
	if line.startswith(">"):
		isempty=False
		if line.rstrip() == ">":
			isempty=True
		else:
			outputfile.write(line)
	else:
		if isempty==False:
			outputfile.write(line)

fastafile.close()
outputfile.close()

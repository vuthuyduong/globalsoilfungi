#!/usr/bin/env python
import sys
import os
import encodings

blastresultfilename = sys.argv[1] # soil.old.genus.blastresult
outputfilename = sys.argv[2] #soil.genus.blastresult

blastresultfile=open(blastresultfilename)
outputfile=open(outputfilename,"w")

isempty=False
for line in blastresultfile:
	isempty=False
	words=line.split("\t")
	if words[0] != "":
		outputfile.write(line)
blastresultfile.close()
outputfile.close()

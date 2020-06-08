#!/usr/bin/env python
import sys
import os
from Bio import SeqIO

otufilename = sys.argv[1]
maxlineNo= 1000
if len(sys.argv) > 2:
	maxlineNo =int(sys.argv[2])

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

resultfolder="funguild"
if os.path.isdir(resultfolder)==False:
	os.system("mkdir " + resultfolder)
inputfile = open(otufilename)
header=inputfile.readline()
i=0
n=0
outputfilename= resultfolder + "/" + GetBase(otufilename) + ".0.txt"
outputfile=open(outputfilename,"w")
outputfile.write(header)
for line in inputfile:
	if i==maxlineNo:
		outputfile.close()
		n=n+1
		i=0
		outputfilename = resultfolder + "/" + GetBase(otufilename) + "." + str(n) + ".txt"
		outputfile=open(outputfilename,"w")
		outputfile.write(header)
	outputfile.write(line)
	i=i+1
inputfile.close()
outputfile.close()

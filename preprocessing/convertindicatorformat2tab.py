#!/usr/bin/env python
import sys
import os
import encodings

indicatorresult = sys.argv[1] # species indicator result file
outputfilename = sys.argv[2] #soil.fasta

indicatorfile=open(indicatorresult)
outputfile=open(outputfilename,"w")
start=False
for line in indicatorfile:
	
	if "Group" in line:
		start=True		
		outputfile.write(line)
	elif start==True:
		line=line.rstrip()
		if (not "p.value" in line):
			newline=""
			words=line.split(" ")
			for word in words:
				word=word.replace("*","").replace(" ","")
				if word!="":
					if newline=="":
						newline=word
					else:
						newline=newline + "\t" + word
			outputfile.write(newline + "\n")
indicatorfile.close()
outputfile.close()

#!/usr/bin/env python
import sys
import os
from Bio import SeqIO

otufilename = sys.argv[1]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

inputfile = open(otufilename)
outputfilename= GetBase(otufilename) + ".funguild.txt"
outputfile=open(outputfilename,"w")
header=inputfile.readline()
words=header.split("\t")

p_s=-1
p_g=-1
p_f=-1
p_o=-1
p_c=-1
p_p=-1
p_k=-1
i=0
for word in words:
	if word.lower()=="kingdom":
		p_k=i
	if word.lower()=="phylum":
		p_p=i	
	if word.lower()=="class":
		p_c=i	
	if word.lower()=="order":
		p_o=i	
	if word.lower()=="family":
		p_f=i	
	if word.lower()=="genus":
		p_g=i	
	i=i+1
outputfile.write("OTU_ID\ttaxonomy\n")
for line in inputfile:
	line=line.rstrip()
	words=line.split("\t")
	taxa=str(round(float(words[2])*100,2)) + "%|"
	kingdom=""
	level =int(words[3])
	if p_k< len(words) and level >=1:
		kingdom=words[p_k]
	if kingdom=="":
		kingdom="unidentified"		   
	taxa=taxa + "k__" + kingdom + ";"
	phylum=""
	if p_p < len(words) and level >=2:
		phylum=words[p_p]
	if phylum=="":
		phylum="unidentified"
	taxa=taxa + "p__" + phylum + ";"
	bioclass=""
	if p_c<len(words) and level >=3:
		bioclass=words[p_c]
	if bioclass=="":
		bioclass="unidentified"
	taxa=taxa + "c__" + bioclass + ";"
	order=""
	if p_o<len(words) and level >=4:
		order=words[p_o]
	if order=="":
		order="unidentified"
	taxa=taxa + "o__" + order + ";"
	family=""
	if p_f<len(words) and level >=5:
		family=words[p_f]
	if family=="":
		family="unidentified"
	taxa=taxa + "f__" + family + ";"
	genus=""
	if p_g<len(words) and level >=6:
		genus=words[p_g]
	if genus=="":
		genus="unidentified"
	
	taxa=taxa + "g__" + genus + ";"
	seqid=words[0].split("|")[0]
	outputfile.write(seqid + "\t" + taxa + "\n")
inputfile.close()
outputfile.close()
print("The output is saved in file: " + outputfilename)

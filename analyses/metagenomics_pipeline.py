#!/usr/bin/env python
# FILE: addtexttosequenceheaders.py, requiring Python 2.7 if classification file contains non-ascii characters
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys, argparse
import numpy as np
import os
from Bio import SeqIO
import json
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='metagenomics_pipepline.py',  
							   usage="%(prog)s [options] -i inputfolder -o output", 
							   description='''Script that t . ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the folder containing all fastq files.')
parser.add_argument('-o','--output', required=True, help='the output folder.') 
parser.add_argument('-metadata','--metadata', required=True, help='the text file containing information about the samples in the input folder.') 
parser.add_argument('-sampleidpos','--sampleidposition', default=0, type=int, required=True, help='the sample id position in the metadata file.') 



args=parser.parse_args()
inputfolder= args.input
outputfolder=args.output 
metadata=args.metadata
sampleidpos=args.sampleidposition

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))] 

def mergepairs(fastq1,fastq2,fasta):
	print("Merging fastq files " + fastq1 + " and " + fastq2 + "...")
	command="time vsearch --fastq_mergepairs " + fastq1 + " --reverse " + fastq2 + " --fastaout " + fasta + " --fastq_allowmergestagger"
	print(command)
	os.system(command)

def unique(fastainput,fastaoutput):
	print("removing duplicates for " + fastainput + "...")
	command="time vsearch --derep_prefix " + fastainput + " --output " + fastaoutput
	print(command)
	os.system(command)
	
def sortfasta(fastainput,fastaoutput):
	print("sorting file " + fastainput + "...")
	command="vsearch --sortbysize " + fastainput + " --output " + fastaoutput
	print(command)
	os.system(command)
	
def cluster(fastainput,threshold,clusterfolder,otufilename):
	#cluster the sequencces
	print("cluster the sequences...")	
	if not os.path.exists(clusterfolder):
		os.system("mkdir " + clusterfolder)
	command="vsearch -cluster_smallmem " + fastainput + " -usersort -id " + str(threshold) + " -centroids " + otufilename + " -clusters " + clusterfolder
	print(command)
	os.system(command)
	return clusterfolder

def createnewsequenceids(code,fastainput,fastaoutput):
	fastafile=open(fastainput)
	newfile=open(fastaoutput,"w")
	count=0
	for line in fastafile:
		if line.startswith(">"):
			count=count+1
			newfile.write(">" + code + "|" + str(count) + "\n")
		else:
			newfile.write(line)
	fastafile.close()
	newfile.close()
	
def createOTUtable(otufastafile,clusterfolder,sampleids,otutablename):
	otutabledict={}
	otus=SeqIO.to_dict(SeqIO.parse(otufastafile,"fasta"))
	for seqid in otus.keys():
		otutabledict.setdefault(seqid,{})
		for sampleid in sampleids:
			otutabledict[seqid].setdefault(sampleid,0)
	for clustername in os.listdir(clusterfolder):
		clusterdict=SeqIO.to_dict(SeqIO.parse(clusterfolder+"/"+clustername,"fasta"))
		repid=list(clusterdict.keys())[0]
		for seqid in clusterdict.keys():
			sampleid=seqid.split("|")[0]
			if sampleid in sampleids:
				otutabledict[repid][sampleid]=otutabledict[repid][sampleid] + 1
	otutablefile=open(otutablename,"w") 
	header="OTUID"
	for sampleid in sampleids:
		header=header + "\t" + sampleid
	header=header + "\n"	
	otutablefile.write(header)
	for otuid in otutabledict.keys():
		line=otuid 
		for sampleid in sampleids:
			line = line  + "\t" + str(otutabledict[otuid][sampleid])
		line =line + "\n"
		otutablefile.write(line)
	return otutablename

def LoadSampleIds(metadata,sampleidpos,allsampleids):
	sampleids=[]
	metadatafile=open(metadata)
	next(metadatafile)
	for line in metadatafile:
		texts=line.split("\t")
		sampleid=""
		if sampleidpos >0 and sampleidpos < len(texts):
			sampleid=texts[sampleidpos].rstrip()
		if sampleid=="":
			continue 
		if sampleid in allsampleids:
			sampleids.append(sampleid)		
	metadatafile.close()
	return sampleids
	
def createfastafile(filedict,outputfolder,sampleids,tablename):	
	for code in filedict.keys():
		fasta=outputfolder + "/" + code + ".fasta"
		fastq1=filedict[code]["R1"]
		fastq2=filedict[code]["R2"]	
		fasta=filedict[code]["fasta"]
		uniquefasta=filedict[code]["uniquefasta"]
		cleanfasta=filedict[code]["cleanfasta"]
		if not os.path.exists(filedict[code]["fasta"]):
			#merge fastq files
			mergepairs(fastq1,fastq2,fasta)
		if not os.path.exists(filedict[code]["uniquefasta"]):	
			#removing duplicates
			unique(fasta,uniquefasta)
		if not os.path.exists(filedict[code]["cleanfasta"]):	
			createnewsequenceids(code,uniquefasta,cleanfasta)
	#concatenate all unique fasta files:
	allfilename=outputfolder + "/all.fasta"
	if not os.path.exists(allfilename):
		command="cat " + outputfolder + "/*.clean.fasta > " + allfilename
		os.system(command)
	#sort unique fasta
	sortedfilename=outputfolder + "/sorted.fasta"
	if not os.path.exists(sortedfilename):
		sortfasta(allfilename,sortedfilename)
	#cluster the sequences
	otufilename=outputfolder + "/otus.fasta"
	clusterfolder=outputfolder + "/clusters"
	if not os.path.exists(otufilename):
		cluster(sortedfilename,0.99,clusterfolder,otufilename)
	#save otutable
	otutablename=outputfolder + "/" + tablename
	if not os.path.exists(otutablename):
		createOTUtable(otufilename,clusterfolder,sampleids,otutablename)	
		print("The otu table is saved in file " + tablename + ".")
	return otufilename

#########MAIN
swabcodes=[]
filenames= os.listdir(inputfolder)
filedict={}
if not os.path.exists(outputfolder):
	os.system("mkdir " + outputfolder)
for filename in filenames:
	texts=filename.split("_")
	if texts[2]=="ITS" or ("56061" in filename) or ("56062" in filename):
		code=texts[0]+"_"+texts[1]
		if not code in filedict.keys():
			filedict.setdefault(code,{})
		if "_R1_" in filename:				
			filedict[code]["R1"]=inputfolder + "/" + filename
		else:
			filedict[code]["R2"]=inputfolder + "/" + filename
		filedict[code]["fasta"]=outputfolder + "/" + code + ".fasta"	
		filedict[code]["uniquefasta"]=outputfolder + "/" + code + ".unique.fasta"	
		filedict[code]["cleanfasta"]=outputfolder + "/" + code + ".clean.fasta"	
		
#load sampleids
sampleids=LoadSampleIds(metadata,sampleidpos,list(filedict.keys()))
tablename=os.path.basename(metadata) + ".table"
#create fasta file
otufilename=createfastafile(filedict,outputfolder,sampleids,tablename)


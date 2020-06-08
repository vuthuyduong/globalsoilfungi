#!/usr/bin/env python
import sys
import os
import encodings
from Bio import SeqIO

fastafilename = sys.argv[1] # soil_ITS2.fasta
tablefilename = sys.argv[2] #soil.table
classificationfilename =sys.argv[3] #soil_ITS2.classification
guildfilename =sys.argv[4] # soil_ITS2.funguilds.guild.txt
outputtable=sys.argv[5] #soil_ITS2.table
outputguild =sys.argv[6] #soil_ITS.guild

def LoadClassification(classificationfilename):
	taxa=[]
	classificationfile=open(classificationfilename)
	header=classificationfile.readline()
	seqids=[]
	p_s=-1
	p_g=-1
	p_f=-1
	p_f=-1
	p_f=-1
	p_c=-1
	p_p=-1
	p_k=-1
	p_l=-1
	p_score=-1
	i=0
	
	for word in header.split("\t"):
		if word.rstrip().lower()=="species":
			p_s=i
		if word.rstrip().lower()=="genus":
			p_g=i
		if word.rstrip().lower()=="family":
			p_f=i
		if word.rstrip().lower()=="order":
			p_o=i
		if word.rstrip().lower()=="class":
			p_c=i
		if word.rstrip().lower()=="phylum":
			p_p=i	
		if word.rstrip().lower()=="kingdom":
			p_k=i		
		if word.rstrip().lower()=="classification position":
			p_l=i	
		if word.rstrip().lower()=="score":
			p_score=i		
		i=i+1

	for line in classificationfile:
		words=line.rstrip().split("\t")
		seqid=words[0].split("|")[0]
		
		kingdom="unidentified"
		if p_k<len(words):
			kingdom=words[p_k]
		phylum="unidentified"	
		if p_p<len(words):
			phylum=words[p_p]
		bioclass="unidentified"	
		if p_c<len(words):	
			bioclass=words[p_c]
		order="unidentified"	
		if p_o<len(words):	
			order=words[p_o]
		family="unidentified"	
		if p_f<len(words):
			family=words[p_f]
		genus="unidentified"	
		if p_g<len(words):	
			genus=words[p_g]
		species="unidentified"
		if p_s<len(words):
			species=words[p_s]
		level=int(words[p_l])
		score=words[p_score]
		if level<7:
			species="unidentified"
		if level<6:
			genus="unidentified"	
		if level<5:
			family="unidentified"	
		if level<4:
			order="unidentified"	
		if level<3:
			bioclass="unidentified"
		if level<2:
			phylum="unidentified"	
		if level<1:
			kingdom="unidentified"	
		classification=kingdom +"\t" +phylum+"\t" + 	bioclass +"\t"+order+"\t"+ family + "\t" + genus + "\t" + species+"\t"+score
		taxa.append(classification)
		seqids.append(seqid)
	classificationfile.close()
	header="kingdom\tphylum\tbioclass\torder\tfamily\tgenus\tspecies\tscore"		
	return seqids,taxa,header

def LoadGuild(guildfilename):
	guildfile=open(guildfilename)
	header=guildfile.readline()
	header=header[header.index("\t")+1:]
	header=header[header.index("\t")+1:]
	header=header[header.index("\t")+1:]
	header=header[header.index("\t")+1:]
	seqids=[]
	guilds=[]
	for line in guildfile:
		seqid=line.split("\t")[0]
		line=line[line.index("\t")+1:]
		line=line[line.index("\t")+1:]
		line=line[line.index("\t")+1:]
		line=line[line.index("\t")+1:]
		guilds.append(line)
		seqids.append(seqid)
	return seqids,guilds,header
def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))] 

taxonids,taxa,header=LoadClassification(classificationfilename)
guildids,guilds,guildheader=LoadGuild(guildfilename)

featurenumber=len(guildheader.split("\t"))
emptyguild=""
for i in range(0,featurenumber-1):
	emptyguild=emptyguild + "\t"
emptyguild=emptyguild+"\n"	

emptyclassification=""
for i in range(0,len(header.split("\t"))-1):
	emptyclassification=emptyclassification + "\t"
	
header="OTU_ID\t" + header + "\t" + guildheader

outputguildfile=open(outputguild,"w")
outputguildfile.write(header)
seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
seqids=[]
for seqrec in seqrecords:
	seqid=seqrec.id.split("|")[0]
	seqids.append(seqid)

outputfile=open(outputtable,"w")
tablefile=open(tablefilename)
header=tablefile.readline()
outputfile.write(header)
for line in tablefile:
	seqid=line.split("\t")[0]
	if seqid in seqids:
		outputfile.write(line)
		classification=emptyclassification
		if seqid in taxonids:
			classification=taxa[taxonids.index(seqid)]
		guild=emptyguild
		if seqid in guildids:
			guild=guilds[guildids.index(seqid)]
		guild=seqid + "\t" + classification + "\t" + guild
		outputguildfile.write(guild)
		
outputguildfile.close()
outputfile.close()
tablefile.close()


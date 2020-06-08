#!/usr/bin/env python
import sys
import os
import encodings

sciencefilename = sys.argv[1] # Tedersoo-DataFileS1 file

fastafilename = "soil.fasta"
classificationfilename = "soil.classification"
otustable="soil.table"
metadata="soil.metadata"

sciencefile = open(sciencefilename,encoding="utf-8", errors='ignore')
fastafile=open(fastafilename,"w")
classificationfile=open(classificationfilename,"w")
otufile=open(otustable,"w")
metafile=open(metadata,"w")
count=0
metafile.write("SampleID\tCountry\tBiome\tRichness\tn")
sampleids=[]
countries=[]
biomes=[]
richness=[]
otufile.write("OTU_ID")
classificationfile.write("Sequence_ID\tOTU_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")
for line in sciencefile:
	#line=unicode(line,errors='ignore')
	words=line.split("\t")
	if count==0:
		for i in range(18, len(words)):
			sampleids.append(words[i].rstrip())	
	elif count==1:
		for i in range(18, len(words)):
			countries.append(words[i].rstrip())	
	elif count==2:	
		for i in range(18, len(words)):
			biomes.append(words[i].rstrip())
	elif count==3:	
		for i in range(18, len(words)):
			richness.append(words[i].rstrip())
		#write to otufile and metafile
		i=0
		for sampleid in sampleids:
			metafile.write(sampleid + "\t" + countries[i] + "\t" + biomes[i] + "\t" + richness[i] + "\n")
			otufile.write("\t" + sampleid)
			i =i+1
		otufile.write("\n")	
	if count <4:
		count=count+1
		continue
	if line.rstrip()=="":
		continue
	seqid=words[0]	
	seq=words[2]
	if seqid =="" or seq=="":
		continue 
	
	fastafile.write(">" + seqid + "\n")
	fastafile.write(seq + "\n")	
	otuid=words[17]
	#otufile.write(otuid)
	otufile.write(seqid)
	for i in range(18, len(words)):
		otufile.write("\t" + words[i].rstrip())
	otufile.write("\n")
	kingdom=words[9]
	phylum=words[10]
	bioclass=words[11]
	order=words[12]
	family=words[13]
	genus=words[14]
	species=words[16]
	classificationfile.write(seqid + "\t" + otuid + "\t" + kingdom + "\t" + phylum + "\t" + bioclass  + "\t" + order + "\t" + family + "\t" + genus + "\t" + species  + "\n")
	count=count+1
	
sciencefile.close()
fastafile.close()
classificationfile.close()
otufile.close()
metafile.close()

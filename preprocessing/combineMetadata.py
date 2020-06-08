#!/usr/bin/env python
import sys
import os
import encodings
from Bio import SeqIO
from  geopy.geocoders import Nominatim

metadatafilename = sys.argv[1] # metadata.txt
otu_metadatafilename = sys.argv[2] #soil.metadata
output =sys.argv[3] #soil_ITS2.metadata

def LoadMetadata(metadatafilename):
	samplenames=[]
	header="pH\tN\tC\tElevation\tLatitude\tLongtitude"
	metadata=[]
	metadatafile=open(metadatafilename)
	metadatafile.readline()
	for line in metadatafile:
		words=line.split("\t")
		samplename=words[1]
		pH=words[2]
		N=words[3]
		C=words[4]
		elevation=words[5]
		LL=words[9]
		latitude=""
		longtitude=""
		if " " in LL:
			latitude=LL.split(" ")[0]
			longtitude=LL.split(" ")[1]
		else:
			latitude=LL
			location=words[17]
			print(location+"\n")
			location=location.replace(" a Que?i","")
			location=location.replace("Lag. Mascardi","Lago Mascardi")
			location=location.replace("Lag. Los Moscos","Lago Los Moscos")
			location=location.replace("Arroyo Los Notros","Los Notros")
			location=location.replace("Colombia: Varillal: El Zafire","Peru:Varillal")
			country=location.split(":")[0]
			city=location.split(":")[1]
			geolocator = Nominatim()
			loc = geolocator.geocode(city+','+ country)			 
			latitude=str(loc.latitude)
			longtitude=str(loc.longitude)	
			print(LL + "\t" + latitude + "\t" + longtitude + "\n") 
		if samplename not in samplenames:
			samplenames.append(samplename)
			metadata.append(pH+ "\t" + N + "\t" + C + "\t" + elevation+"\t"+latitude+"\t"+longtitude)
	return samplenames,metadata,header


samplenames,metadata,header=LoadMetadata(metadatafilename)
otu_metadatafile=open(otu_metadatafilename)
outputfile=open(output,"w")
header=otu_metadatafile.readline().rstrip() + "\t" + header + "\n"
outputfile.write(header)
for line in otu_metadatafile:
	samplename=line.split("\t")[0]
	if samplename in samplenames:
		data=metadata[samplenames.index(samplename)]
		outputfile.write(line.replace("\n","") + "\t" + data + "\n")
	else:
		outputfile.write(line.rstrip() + "\t\t\t\t\t\t\n")		
outputfile.close()
otu_metadatafile.close()


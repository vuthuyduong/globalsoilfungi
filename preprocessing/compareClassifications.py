#!/usr/bin/env python
import sys
import os
import encodings
from Bio import SeqIO

UNITEclassificationfilename = sys.argv[1] # UNITE_classification.txt
cbsclassificationfilename = sys.argv[2] #soil_ITS2.classification
CBSaccessionfilename = sys.argv[3] # ITSacc.txt
remainedfastafilename=sys.argv[4] # as output
output="soil_ITS2.classification.compared"
newbycbs="newbycbs.txt"


#load remaining sequences
remainedsequences=list(SeqIO.parse(remainedfastafilename, "fasta"))
remainedseqids=[]
for seq in remainedsequences:
	seqid=seq.id
	if "|" in seqid:
		seqid=seqid.split("|")[0]
	remainedseqids.append(seqid)
	
#load cbs accession ids
def LoadCBSaccessionids(CBSaccessionfilename):
	ids=[]
	CBSaccessionfile=open(CBSaccessionfilename)
	for line in CBSaccessionfile:
		words=line.split("\t")
		ids.append(words[0].replace(".1",""))
	return ids

ids=LoadCBSaccessionids(CBSaccessionfilename)
#load CBS+UNITE classification
cbsseqids=[]
cbssims=[]
cbsspecies=[]
cbsgenera=[]
cbsfamilies=[]
cbsorders=[]
cbsclasses=[]
cbsphyla=[]
cbskingdoms=[]
cbsrefids=[]
newcbsspecies=[]
newcbsgenera=[]
newcbsfamilies=[]
newcbsorders=[]
newcbsclasses=[]
newcbsphyla=[]
newcbskingdoms=[]
cbsclassificationfile=open(cbsclassificationfilename)
outputfile=open(output,"w")

next(cbsclassificationfile)
c_s=0
c_g=0
c_f=0
c_o=0
c_c=0
c_p=0
c_k=0
uc_s=0
uc_g=0
uc_f=0
uc_o=0
uc_c=0
uc_p=0
uc_k=0
for line in cbsclassificationfile:
	words=line.split("\t")
	if words[0]=="":
		continue
	seqid=words[0].split("|")[0]
	cbsseqids.append(seqid)
	sim=0
	species=""
	genus=""
	family=""
	order=""
	bioclass=""
	phylum=""
	kingdom=""
	refid=words[1]
	sim=float(words[2])
	cbssims.append(sim)
	if sim >=0.99:
		species=words[12]	
		if species !="":
				uc_s=uc_s+1
				if refid in ids:
					c_s=c_s+1
					if not species in newcbsspecies:
						newcbsspecies.append(species)
	if sim >=0.949:
		genus=words[11]
		if genus !="":
				uc_g=uc_g+1
				if refid in ids:
					c_g=c_g+1
					if not genus in newcbsgenera:
						newcbsgenera.append(genus)					
	if sim >=0.906:
		family=words[10]
		if family !="":
				uc_f=uc_f+1
				if refid in ids:
					c_f=c_f+1
					if not family in newcbsfamilies:
						newcbsfamilies.append(family)
	if sim>=0.826:
		order=words[9]
		if order !="":
				uc_o=uc_o+1
				if refid in ids:					
					c_o=c_o+1
					if not order in newcbsorders:
						newcbsorders.append(order)
	if sim>=0.535:
		bioclass =words[8]
		if bioclass !="":
				uc_c=uc_c+1
				if refid in ids:
					c_c=c_c+1
					if not bioclass in newcbsclasses:
						newcbsclasses.append(bioclass)
	if sim>=0.4:
		phylum=words[7]
		if phylum !="":
				uc_p=uc_p+1
				if refid in ids:
					c_p=c_p+1
					if not phylum in newcbsphyla:
						newcbsphyla.append(phylum)
	if sim>=0.2:
		kingdom=words[6]
		if kingdom !="":
				uc_k=uc_k+1
				if refid in ids:
					c_k=c_k+1
					if not kingdom in cbskingdoms:
						newcbskingdoms.append(kingdom)	
	cbsspecies.append(species)
	cbsgenera.append(genus)
	cbsfamilies.append(family)
	cbsorders.append(order)
	cbsclasses.append(bioclass)
	cbsphyla.append(phylum)
	cbskingdoms.append(kingdom)
	cbsrefids.append(refid)
#write to newbycbs.txt
newbycbsfile=open(newbycbs,"w")
newbycbsfile.write("Number of sequences classified at species level by CBS: " + str(len(newcbsspecies)) + "\n")
for species in newcbsspecies:
	newbycbsfile.write(species + "\n")
newbycbsfile.write("Number of sequences classified at genus level by CBS: " + str(len(newcbsgenera)) + "\n")
for genus in newcbsgenera:
	newbycbsfile.write(genus + "\n")
newbycbsfile.write("Number of sequences classified at family level by CBS: " + str(len(newcbsfamilies)) + "\n")
for family in newcbsfamilies:
	newbycbsfile.write(family + "\n")
newbycbsfile.write("Number of sequences classified at order level by CBS: " + str(len(newcbsorders)) + "\n")
for order in newcbsorders:
	newbycbsfile.write(order + "\n")
newbycbsfile.write("Number of sequences classified at class level by CBS: " + str(len(newcbsclasses)) + "\n")
for bioclass in newcbsclasses:
	newbycbsfile.write(bioclass + "\n")
newbycbsfile.write("Number of sequences classified at phylum level by CBS: " + str(len(newcbsphyla)) + "\n")
for phylum in newcbsphyla:
	newbycbsfile.write(phylum + "\n")
newbycbsfile.write("Number of sequences classified at kingdom level by CBS: " + str(len(newcbskingdoms)) + "\n")
for kingdom in newcbskingdoms:
	newbycbsfile.write(kingdom + "\n")

	
newbycbsfile.close()
#load UNITE classification
UNITEclassificationfile=open(UNITEclassificationfilename)
uniteseqids=[]
unitesims=[]
unitespecies=[]
unitegenera=[]
unitefamilies=[]
uniteorders=[]
uniteclasses=[]
unitephyla=[]
unitekingdoms=[]
uniterefids=[]

next(UNITEclassificationfile)
u_s=0
u_g=0
u_f=0
u_o=0
u_c=0
u_p=0
u_k=0
u_c_s=0
u_c_g=0
u_c_f=0
u_c_o=0
u_c_c=0
u_c_p=0
u_c_k=0
u_c_all=0
g_s=0 #number of the genera classified the same
g_d=0 #number of the genera classified the same
g_c=0 #number of genera classified by the CBS barcodes only
f_s=0
f_d=0
o_s=0
o_d=0
c_s=0
c_d=0

for line in UNITEclassificationfile:
	words=line.split("\t")
	if words[0]=="":
		continue
	seqid=words[0]
	uniteseqids.append(seqid)
	sim=0
	species=""
	genus=""
	family=""
	order=""
	bioclass=""
	phylum=""
	kingdom=""
	refid=words[15]
	if "%" in words[5]:
		sim=float(words[5].replace("%",""))/100
	unitesims.append(sim)
	if sim >=0.99:
		species=words[16]	
		if species !="":
				u_s=u_s+1
				if (seqid not in cbsseqids) and (seqid not in remainedseqids):
					u_c_s=u_c_s+1
	if sim >=0.9:
		genus=words[14]
		if genus !="":
				u_g=u_g+1
				if seqid in cbsseqids:
					index=cbsseqids.index(seqid)
					genus_c=cbsgenera[index]
					if genus==genus_c:
						g_s=g_s+1
					else:
						g_d=g_d+1
				if (seqid not in cbsseqids) and (seqid not in remainedseqids):
					u_c_g=u_c_g+1					
	if sim >=0.85:
		family=words[13]
		if family !="":
				u_f=u_f+1
				if seqid in cbsseqids:
					index=cbsseqids.index(seqid)
					family_c=cbsfamilies[index]
					if family==family_c:
						f_s=f_s+1
					else:
						f_d=f_d+1
				if (seqid not in cbsseqids) and (seqid not in remainedseqids):
					u_c_f=u_c_f+1				
	if sim>=0.8:
		order=words[12]
		if order !="":
				u_o=u_o+1
				if seqid in cbsseqids:
					index=cbsseqids.index(seqid)
					order_c=cbsorders[index]
					if order==order_c:
						o_s=o_s+1
					else:
						o_d=o_d+1
				if (seqid not in cbsseqids) and (seqid not in remainedseqids):
					u_c_o=u_c_o+1				
	if sim>=0.75:
		bioclass =words[11]
		if bioclass !="":
				u_c=u_c+1
				if seqid in cbsseqids:
					index=cbsseqids.index(seqid)
					class_c=cbsclasses[index]
					if bioclass==class_c:
						c_s=c_s+1
					else:
						c_d=c_d+1
				if (seqid not in cbsseqids) and (seqid not in remainedseqids):
					u_c_c=u_c_c+1				
	if sim>=0.7:
		phylum=words[10]
		if phylum !="":
				u_p=u_p+1
				if (seqid not in cbsseqids) and (seqid not in remainedseqids):
					u_c_p=u_c_p+1
	if sim>=0.65:
		kingdom=words[9]
		if kingdom !="":
				u_k=u_k+1
	if (seqid not in cbsseqids) and (seqid not in remainedseqids):
		u_c_all=u_c_all +1 		
	unitespecies.append(species)
	unitegenera.append(genus)
	unitefamilies.append(family)
	uniteorders.append(order)
	uniteclasses.append(bioclass)
	unitephyla.append(phylum)
	unitekingdoms.append(kingdom)
	uniterefids.append(refid)	
#write to output:
outputfile.write("Number of sequences identified by UNITE but are not ITS2 sequences: "  + str(u_c_all) + "\n")
outputfile.write("Number of genera UNITE vs. UNITEbutNotCBS vs. CBS+UNITE: \t"+ str(u_g) + "\t" + str(u_c_g) + "\t" + str(uc_g)+ "\n")
outputfile.write("Number of families UNITE vs. UNITEbutNotCBS vs. CBS+UNITE: \t" + str(u_f) + "\t" + str(u_c_f)  + "\t" + str(uc_f)+ "\n")
outputfile.write("Number of orders UNITE vs. UNITEbutNotCBS vs. CBS+UNITE: \t" + str(u_o) + "\t" + str(u_c_o)  + "\t" + str(uc_o)+ "\n")
outputfile.write("Number of classes UNITE vs. UNITEbutNotCBS vs. CBS+UNITE: \t" + str(u_c) + "\t" + str(u_c_c)  + "\t" + str(uc_c)+ "\n")

outputfile.write("Number of the same genera: \t"+ str(g_s) + "\n")
outputfile.write("Number of different genera: \t"+ str(g_d) + "\n")
outputfile.write("Number of the CBS genera: \t"+ str(c_g) + "\n")
outputfile.write("Number of the UNITE_ITS2 genera: \t"+ str(u_g-u_c_g-g_s-g_d) + "\n")
outputfile.write("Number of the new UNITE_ITS2 genera: \t"+ str(uc_g-c_g-g_s-g_d) + "\n")

outputfile.write("Number of the same families: \t"+ str(f_s) + "\n")
outputfile.write("Number of different families: \t"+ str(f_d) + "\n")
outputfile.write("Number of the CBS families: \t"+ str(c_f) + "\n")
outputfile.write("Number of the UNITE_ITS2 families: \t"+ str(u_f-u_c_f-f_s-f_d) + "\n")
outputfile.write("Number of the new UNITE_ITS2 families: \t"+ str(uc_f-c_f-f_s-f_d) + "\n")

outputfile.write("Number of the same orders: \t"+ str(o_s) + "\n")
outputfile.write("Number of different orders: \t"+ str(o_d) + "\n")
outputfile.write("Number of the CBS orders: \t"+ str(c_o) + "\n")
outputfile.write("Number of the UNITE_ITS2 genera: \t"+ str(u_o-u_c_o-o_s-o_d) + "\n")
outputfile.write("Number of the new UNITE_ITS2 genera: \t"+ str(uc_o-c_o-o_s-o_d) + "\n")


outputfile.write("Number of the same classes: \t"+ str(c_s) + "\n")
outputfile.write("Number of different classes: \t"+ str(c_d) + "\n")
outputfile.write("Number of the CBS classes: \t"+ str(c_c) + "\n")
outputfile.write("Number of the UNITE_ITS2 classes: \t"+ str(u_c-u_c_c-c_s-c_d) + "\n")
outputfile.write("Number of the new UNITE_ITS2 classes: \t"+ str(uc_c-c_c-c_s-c_d) + "\n")

outputfile.write("seq_id\tunite_refid\tunite_score\tunite_kingdom\tunite_phylum\tunite_class\tunite_order\tunite_family\tunite_genus\tunite_species\tcbs_refid\tcbs_score\tcbs_kingdom\tcbs_phylum\tcbs_class\tcbs_order\tcbs_family\tcbs_genus\tcbs_species\n")

i=0
for seqid in uniteseqids:
	outputfile.write(seqid + "\t" + uniterefids[i] + "\t" + str(unitesims[i]) + "\t" + unitekingdoms[i] + "\t" + unitephyla[i] + "\t" + uniteclasses[i] + "\t" + uniteorders[i] + "\t" + unitefamilies[i] + "\t" + unitegenera[i] + "\t" + unitespecies[i] + "\t")
	if seqid in cbsseqids:
		j=cbsseqids.index(seqid)
		outputfile.write(cbsrefids[j] + "\t" + str(cbssims[j]) + "\t" + cbskingdoms[j] + "\t" + cbsphyla[j] + "\t" + cbsclasses[j] + "\t" + cbsorders[j] + "\t" + cbsfamilies[j] + "\t" + cbsgenera[j] + "\t" + cbsspecies[j] + "\n")
	else:
		outputfile.write("\t\t\t\t\t\t\t\t\n")
	i=i+1	
UNITEclassificationfile.close()
cbsclassificationfile.close()
outputfile.close()

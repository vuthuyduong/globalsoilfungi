#####################
#Data loading ----
#setwd("/home/duong/ProgLang/Python/dnabarcoder/alldata/soil_samplesUNITE")
#setwd("/home/duong/Data/Metagenomics/references/yeastCBSITS2/UNITEsoilITS2")
setwd("C:/Users/Duong Vu/Documents/CBSPapers/DnaBarcoder_2022/Data/soil_samplesUNITE")
#vegan
#install.packages("vegan")
#install.packages("ggsignif", type="win.binary") 
#install.packages("ggpubr", type="win.binary") 

#install.packages("ggsignif") 
#install.packages("ggpubr") 

library("ggpubr")
library("vegan")
library(grid)

#install.packages("extrafont")
library(extrafont)
font_import()
#load otu table
allotu_table <-read.delim("soil_ITS2.table",row.names=1,check.names=FALSE)
#meta_table <-read.delim("/home/duong/Data/Metagenomics/soil_samplesUNITE/sciencepaper_otus/soil.unite.metadata",row.names=1,check.names=FALSE)
meta_table <-read.delim("soil_ITS2.metadata",row.names=1,check.names=FALSE)
#classification_table <-read.delim("dnabarcoder/soil_ITS2.moldITS2_BLAST.classification.better",row.names=1,check.names=FALSE)
classification_table <-read.delim("dnabarcoder/soil_ITS2.moldITS2_BLAST.classification.worseclassificationremoved",row.names=1,check.names=FALSE)
otu_table <- allotu_table[rownames(classification_table), ]
guild_table <-read.delim("dnabarcoder/soil_ITS2.CBSITS_BLAST.classification.worseclassificationremoved.guilds.txt",row.names=1,check.names=FALSE)
#guild_table<-guild_table[rownames(classification_table), ]
  
#save class_tables
#write.table(classification_table,file = "results/classification_table.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(otu_table,file = "results_moldITS2/moldotutable.txt", sep = "\t", quote = FALSE, row.names = T)
#write.table(guild_table,file = "results/guild_table.txt", sep = "\t", quote = FALSE, row.names = T)
#write.table(meta_table,file = "results/meta_table.txt", sep = "\t", quote = FALSE, row.names = T)



otu_number <-dim(otu_table)[1]
sample_number <-dim(otu_table)[2]
biomes <-as.character(meta_table$Biome)
countries <- meta_table$Country
labels <-paste(biomes," ",countries)
r_sum <- rowSums(otu_table[,1:sample_number])
c_sum <-colSums(otu_table[1:otu_number,])

#reduce the biome annotation 
biomes[biomes=="Arctic tundra"] <- "AT"
biomes[biomes=="Boreal forests"] <- "BF"
biomes[biomes=="Grassland and shrubland"] <- "GS"
biomes[biomes=="Mediterranean"] <- "MED"
biomes[biomes=="Moist tropical forests"] <- "MTF"
biomes[biomes=="Savannas"] <- "SAV"
biomes[biomes=="Dry tropical forests"] <- "DTF"
biomes[biomes=="Southern temperate forests"] <- "STF"
biomes[biomes=="Tropical montane forests"] <- "TMF"
biomes[biomes=="Temperate deciduous forests"] <- "TDF"
biomes[biomes=="Temperate coniferous forests"] <- "TCF"

#create results folder

dir.create("results_moldITS2")

#########################################
#Taxonomy at class level with biome merged----
class_table <-data.frame(OTU_ID = row.names(classification_table), 
                         class = as.character(classification_table$class),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
class_table <-class_table[order(class_table$class),]
#replace empty class with "unidentified"
class_table$class[class_table$class==""] <- "unidentified"
class_table <-class_table[order(class_table$class),]
#merge all the class with the same name
class_table <- aggregate(class_table[1:sample_number+1], by=list(class_table$class), FUN = "sum")
classnames <- class_table$Group.1
class_table <- t(class_table)
colnames(class_table)=classnames
class_table=class_table[-1,]
#save class_table
write.table(class_table,file = "results_moldITS2/class_level.txt", sep = "\t", quote = FALSE, row.names = T)
#Get taxalist with respect to the order of all samples
class(class_table) <- "double"
x <- class_table/rowSums(class_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)
df <- data.frame(Sample_ID=rownames(class_table), 
                 biome=biomes,
                 class_table,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE
)
#change the number to numeric
df[,2:ncol(df)]= apply(df[,2:ncol(df)], 2, function(x) as.numeric(as.character(x)))
#merge df with biomet
df <-df[order(df$biome),]
df <- aggregate(df[2:ncol(df)], by=list(df$biome), FUN = "sum")

rownames(df)=df$Group.1
df=df[,-1]
#count all classes
s <- colSums(df[,1:ncol(df)])
df_all <- data.frame(t(s))
row.names(df_all) = c("ALL")
df <- rbind(df,df_all)

write.table(df,file = "results_moldITS2/class_level_biome.txt", sep = "\t", quote = FALSE, row.names = T)

x <- df/rowSums(df)
x<-x[,taxa_list]
taxa_list <-colnames(x)
if (is.element("unidentified",taxa_list)==TRUE){
  #remove "__unidentified__" and add it to others
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
}
new_x<-x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Classes=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
#colours[N]="#FFFF00" #yellow
colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Classes))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
p<- p + theme(text=element_text(size=20))
pdf("results_moldITS2/class_level_biome.pdf",height=6,width=21)
print(p)
dev.off()

###############################################
#Taxonomy at class level ----
class_table <-data.frame(OTU_ID = row.names(classification_table), 
                         class = as.character(classification_table$class),
                         otu_table,
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
class_table <-class_table[order(class_table$class),]
#replace empty class with "unidentified"
class_table$class[class_table$class==""] <- "unidentified"
class_table <-class_table[order(class_table$class),]
#merge all the class with the same name
class_table <- aggregate(class_table[1:sample_number+1], by=list(class_table$class), FUN = "sum")
classnames <- class_table$Group.1
class_table <- t(class_table)
colnames(class_table)=classnames
class_table=class_table[-1,]
#save class_table
write.table(class_table,file = "results_moldITS2/class_level.txt", sep = "\t", quote = FALSE, row.names = T)
class(class_table) <- "double"

x <- class_table/rowSums(class_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

#taxa_list <-colnames(x)
#if (is.element("unidentified",taxa_list)==TRUE){
#  #remove "__unidentified__" and add it to others
#  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
#  taxa_list <- c(taxa_list,"unidentified")
#}
new_x<-x[,taxa_list]
df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Classes=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
#colours[N]="#FFFF00"#yellow
colours[N]="#000000" #black
library(ggplot2)
p_c<-ggplot(df,aes(Sample,Value,fill=Classes))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p_c<-p_c+scale_fill_manual(values=colours[1:N])
p_c<-p_c+theme_bw()+ylab("Proportions")
p_c<-p_c+theme_bw()+xlab("Samples")
p_c<-p_c+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p_c<-p_c+theme(axis.text.x=element_blank())
p_c<- p_c + theme(text=element_text(size=20))
#p_c<-p_c+ theme(axis.title.y = element_text(size = 15))
#p_c<-p_c+ theme(axis.text.y = element_text(size = 15))
#p_c<-p_c+ theme(axis.title.x = element_text(size = 15))
#p_c<-p_c+ theme(legend.title = element_text(size = 15))
#p_c<-p_c+ theme(legend.text = element_text(size = 15)) 
#p_c<-p_c+ theme(strip.text = element_text(size = 15)) 
pdf("results_moldITS2/class_level.pdf",height=6,width=21)
print(p_c)
dev.off()

#########################################
#order_level_biome.pdf----
order_table <-data.frame(OTU_ID = row.names(classification_table), 
                         order = as.character(classification_table$order),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
order_table <-order_table[order(order_table$order),]
#replace empty class with "unidentified"
order_table$order[order_table$order==""] <- "unidentified"

#order again the table
order_table <-order_table[order(order_table$order),]

#merge all the order with the same name
order_table <- aggregate(order_table[1:sample_number+1], by=list(order_table$order), FUN = "sum")
ordernames <- order_table$Group.1
order_table <- t(order_table)
colnames(order_table)=ordernames
order_table=order_table[-1,]
#save order_table
write.table(order_table,file = "results_moldITS2/order_level.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the order of all samples
class(order_table) <- "double"
x<-order_table[,order(colSums(order_table),decreasing=TRUE)]
taxa_list <-colnames(x)

x <- order_table/rowSums(order_table)
x <- x[,order(colSums(x),decreasing=TRUE)]


df <- data.frame(Sample_ID=rownames(order_table), 
                 biome=biomes,
                 order_table,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE
)
#change the number to numeric
df[,2:ncol(df)]= apply(df[,2:ncol(df)], 2, function(x) as.numeric(as.character(x)))
#merge df with biomet
df <-df[order(df$biome),]
df <- aggregate(df[2:ncol(df)], by=list(df$biome), FUN = "sum")
rownames(df)=df$Group.1
df=df[,-1]

write.table(df,file = "results_moldITS2/order_level_biome.txt", sep = "\t", quote = FALSE, row.names = T)

selectedtaxa <- taxa_list
if (length(selectedtaxa) >26){
  #select only 26 taxa to show
  selectedtaxa <- taxa_list[1:26]
}
new_df <- data.frame(Sample_ID=colnames(df), 
                     taxa=colnames(df),
                     t(df),
                     row.names = 1, 
                     check.names=FALSE,
                     stringsAsFactors = FALSE
)
new_df$taxa[!(new_df$taxa %in% selectedtaxa)] <- "others"
new_df <-new_df[order(new_df$taxa),]
new_df <- aggregate(new_df[2:ncol(new_df)], by=list(new_df$taxa), FUN = "sum")
rownames(new_df)=new_df$Group.1
new_df=new_df[,-1]
new_df <- t(new_df)

#count all classes
s <- colSums(new_df[,1:ncol(new_df)])
df_all <- data.frame(t(s))
row.names(df_all) = c("ALL")
new_df <- rbind(new_df,df_all)

x <- new_df/rowSums(new_df)

if (length(taxa_list) >26){
  taxa_list <- c(selectedtaxa,"others")
} else {
  taxa_list <- selectedtaxa
}

new_x<-x[,taxa_list]
if (is.element("unidentified",taxa_list)==TRUE){
  #remove "__unidentified__" and add it to others
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
}
new_x<-x[,taxa_list]


df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Orders=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
#colours[N]="#FFFF00" #yellow
colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Orders))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
p<- p + theme(text=element_text(size=20))
pdf("results_moldITS2/order_level_biome.pdf",height=6,width=21)
print(p)
dev.off()


#########################################
#Taxonomy at order level ----
order_table <-data.frame(OTU_ID = row.names(classification_table), 
                         order = as.character(classification_table$order),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
order_table <-order_table[order(order_table$order),]
#replace empty class with "unidentified"
order_table$order[order_table$order==""] <- "unidentified"

#change the name of insignificant groups to others
#s <- rowSums(order_table[,2:sample_number+1])
#order_table$order[order_table$order != "unidentified" & s<300] <- "others"

#order again the table
order_table <-order_table[order(order_table$order),]

#merge all the order with the same name
order_table <- aggregate(order_table[1:sample_number+1], by=list(order_table$order), FUN = "sum")
ordernames <- order_table$Group.1

#save order_table
write.table(order_table,file = "results_moldITS2/order_level.txt", sep = "\t", quote = FALSE, row.names = T)

order_table <- t(order_table)
colnames(order_table)=ordernames
order_table=order_table[-1,]
class(order_table) <- "double"

x <- order_table/rowSums(order_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

#taxa_list <-colnames(x)
#selectedtaxa <- taxa_list
#if (length(selectedtaxa) >26){
#  #select only 26 taxa to show
#  selectedtaxa <- taxa_list[1:26]
#}
#change the name of the taxa not in the selectedlist to others
new_df <- data.frame(Sample_ID=colnames(order_table), 
                     taxa=colnames(order_table),
                     t(order_table),
                     row.names = 1, 
                     check.names=FALSE,
                     stringsAsFactors = FALSE
)
new_df$taxa[!(new_df$taxa %in% selectedtaxa)] <- "others"
new_df <-new_df[order(new_df$taxa),]
#change the number to numeric
new_df[,2:ncol(new_df)]= apply(new_df[,2:ncol(new_df)], 2, function(x) as.numeric(as.character(x)))
new_df <- aggregate(new_df[2:ncol(new_df)], by=list(new_df$taxa), FUN = "sum")
rownames(new_df)=new_df$Group.1
new_df=new_df[,-1]
new_df <- t(new_df)

x <- new_df/rowSums(new_df)
x <- x[,order(colSums(x),decreasing=TRUE)]

#if (length(taxa_list) >26){
#  taxa_list <- c(selectedtaxa,"others")
#} else {
#  taxa_list <- selectedtaxa
#}
#new_x<-x[,taxa_list]
#if (is.element("unidentified",taxa_list)==TRUE){
#  #remove "__unidentified__" and add it to others
#  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
#  taxa_list <- c(taxa_list,"unidentified")
#}
new_x<-x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Orders=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
#colours[N]="#FFFF00"#yellow
colours[N]="#000000" #black
library(ggplot2)
p_o<-ggplot(df,aes(Sample,Value,fill=Orders))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p_o<-p_o+scale_fill_manual(values=colours[1:N])
p_o<-p_o+theme_bw()+ylab("Proportions")
p_o<-p_o+theme_bw()+xlab("Samples")
p_o<-p_o+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p_o<-p_o+theme(axis.text.x=element_blank())
p_o<- p_o + theme(text=element_text(size=20))
#p_o<-p_o+ theme(axis.title.y = element_text(size = 15))
#p_o<-p_o+ theme(axis.text.y = element_text(size = 15))
#p_o<-p_o+ theme(legend.title = element_text(size = 15))
#p_o<-p_o+ theme(legend.text = element_text(size = 15)) 
#p_o<-p_o+ theme(strip.text = element_text(size = 15)) 
pdf("results_moldITS2/order_level.pdf",height=6,width=21)
print(p_o)
dev.off()

#########################################
#family_level_biome.pdf ----
family_table <-data.frame(OTU_ID = row.names(classification_table), 
                          family = as.character(classification_table$family),
                          otu_table,row.names = 1, 
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
family_table <-family_table[order(family_table$family),]
#replace empty class with "unidentified"
family_table$family[family_table$family==""] <- "unidentified"

#change the name of insignificant groups to others
#s <- rowSums(family_table[,2:sample_number+1])
#family_table$family[family_table$family != "unidentified" & s<300] <- "others"

#order again the table
family_table <-family_table[order(family_table$family),]

#merge all the order with the same name
family_table <- aggregate(family_table[1:sample_number+1], by=list(family_table$family), FUN = "sum")
familynames <- family_table$Group.1

family_table <- t(family_table)
colnames(family_table)=familynames
family_table=family_table[-1,]

#save family_table
write.table(family_table,file = "results_moldITS2/family_level.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the order of all samples
class(family_table) <- "double"
x<-family_table[,order(colSums(family_table),decreasing=TRUE)]
taxa_list <-colnames(x)

x <- family_table/rowSums(family_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

df <- data.frame(Sample_ID=rownames(family_table), 
                 biome=biomes,
                 family_table,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE
)
#change the number to numeric
df[,2:ncol(df)]= apply(df[,2:ncol(df)], 2, function(x) as.numeric(as.character(x)))
#merge df with biomet
df <-df[order(df$biome),]
df <- aggregate(df[2:ncol(df)], by=list(df$biome), FUN = "sum")
rownames(df)=df$Group.1
df=df[,-1]

write.table(df,file = "results_moldITS2/family_level_biome.txt", sep = "\t", quote = FALSE, row.names = T)

#select only 26 taxa to show
selectedtaxa <- taxa_list
if (length(selectedtaxa) >26){
  #select only 26 taxa to show
  selectedtaxa <- taxa_list[1:26]
}
new_df <- data.frame(Sample_ID=colnames(df), 
                     taxa=colnames(df),
                     t(df),
                     row.names = 1, 
                     check.names=FALSE,
                     stringsAsFactors = FALSE
)
new_df$taxa[!(new_df$taxa %in% selectedtaxa)] <- "others"
new_df <-new_df[order(new_df$taxa),]
new_df <- aggregate(new_df[2:ncol(new_df)], by=list(new_df$taxa), FUN = "sum")
rownames(new_df)=new_df$Group.1
new_df=new_df[,-1]
new_df <- t(new_df)

#count all classes
s <- colSums(new_df[,1:ncol(new_df)])
df_all <- data.frame(t(s))
row.names(df_all) = c("ALL")
new_df <- rbind(new_df,df_all)

x <- new_df/rowSums(new_df)
#order x by the taxalist to be consistent with all samples

if (length(taxa_list) >26){
  taxa_list <- c(selectedtaxa,"others")
} else {
  taxa_list <- selectedtaxa
}
new_x<-x[,taxa_list]

#remove "__unidentified__" 
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-new_x[,taxa_list]


df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Families=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
#colours[N]="#FFFF00" #yellow
colours[N]="#000000" #black
library(ggplot2)
p_f<-ggplot(df,aes(Sample,Value,fill=Families))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p_f<-p_f+scale_fill_manual(values=colours[1:N])
p_f<-p_f+theme_bw()+ylab("Proportions")
p_f<-p_f+theme_bw()+xlab("Biotopes")
p_f<-p_f+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p_f<-p_f+theme(axis.text.x=element_blank())
p_f<- p_f + theme(text=element_text(size=20))
pdf("results_moldITS2/family_level_biome.pdf",height=6,width=21)
print(p_f)
dev.off()

#########################################
#Taxonomy at family level ----
family_table <-data.frame(OTU_ID = row.names(classification_table), 
                          family = as.character(classification_table$family),
                          otu_table,row.names = 1, 
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
family_table <-family_table[order(family_table$family),]
#replace empty class with "unidentified"
family_table$family[family_table$family==""] <- "unidentified"
#order again the table
family_table <-family_table[order(family_table$family),]

#merge all the order with the same name
family_table <- aggregate(family_table[1:sample_number+1], by=list(family_table$family), FUN = "sum")
familynames <- family_table$Group.1

#save family_table
write.table(family_table,file = "results_moldITS2/family_level.txt", sep = "\t", quote = FALSE, row.names = T)

family_table <- t(family_table)
colnames(family_table)=familynames
family_table=family_table[-1,]

class(family_table) <- "double"
x <- family_table/rowSums(family_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
#taxa_list <-colnames(x)

##select only 26 taxa to show
#selectedtaxa <- taxa_list
#if (length(selectedtaxa) >26){
#  #select only 26 taxa to show
#  selectedtaxa <- taxa_list[1:26]
#}

#change the name of the taxa not in the selectedlist to others
new_df <- data.frame(Sample_ID=colnames(family_table), 
                     taxa=colnames(family_table),
                     t(family_table),
                     row.names = 1, 
                     check.names=FALSE,
                     stringsAsFactors = FALSE
)
new_df$taxa[!(new_df$taxa %in% selectedtaxa)] <- "others"
new_df <-new_df[order(new_df$taxa),]
#change the number to numeric
new_df[,2:ncol(new_df)]= apply(new_df[,2:ncol(new_df)], 2, function(x) as.numeric(as.character(x)))
new_df <- aggregate(new_df[2:ncol(new_df)], by=list(new_df$taxa), FUN = "sum")
rownames(new_df)=new_df$Group.1
new_df=new_df[,-1]
new_df <- t(new_df)

x <- new_df/rowSums(new_df)
x <- x[,order(colSums(x),decreasing=TRUE)]

#if (length(taxa_list) >26){
#  taxa_list <- c(selectedtaxa,"others")
#} else {
#  taxa_list <- selectedtaxa
#}
new_x<-x[,taxa_list]

if (is.element("unidentified",taxa_list)==TRUE){
  #remove "__unidentified__" and add it to others
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
}
new_x<-x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Families=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
#colours[N]="#FFFF00"#yellow
colours[N]="#000000" #black
library(ggplot2)
p_f<-ggplot(df,aes(Sample,Value,fill=Families))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p_f<-p_f+scale_fill_manual(values=colours[1:N])
p_f<-p_f+theme_bw()+ylab("Proportions")
p_f<-p_f+theme_bw()+xlab("Samples")
p_f<-p_f+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p_f<-p_f+theme(axis.text.x=element_blank())
p_f<- p_f + theme(text=element_text(size=20))
#p_f<-p_f+ theme(axis.title.y = element_text(size = 15))
#p_f<-p_f+ theme(axis.text.y = element_text(size = 15))
#p_f<-p_f+ theme(legend.title = element_text(size = 15))
#p_f<-p_f+ theme(legend.text = element_text(size = 15)) 
#p_f<-p_f+ theme(strip.text = element_text(size = 15)) 
pdf("results_moldITS2/family_level.pdf",height=6,width=21)
print(p_f)
dev.off()


###########################
#genus_level_biome.pdf----
genus_table <-data.frame(OTU_ID = row.names(classification_table), 
                         genus = as.character(classification_table$genus),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
genus_table <-genus_table[order(genus_table$genus),]
#replace empty class with "unidentified"
genus_table$genus[genus_table$genus==""] <- "unidentified"

genus_table$genus[genus_table$genus=="Mastigobasidium"] <-     "Leucosporidium" #to be in the same line with the species
genus_table$genus[genus_table$genus=="Schizoblastosporion"] <- "Nadsonia"

#genus_table$genus[genus_table$genus=="Nadsonia"] <-           "Nadsonia                 "
#genus_table$genus[genus_table$genus=="Leucosporidium"] <-     "Leucosporidium           "


#replace unknown classes with "unidentified"
genus_table$genus[startsWith(genus_table$genus,"Branch0")] <- "unidentified"
genus_table$genus[grepl("Incertae_sedis",genus_table$genus)] <- "unidentified"
genus_table$genus[startsWith(genus_table$genus,"GS")] <- "unidentified"
#change the name of insignificant groups to others
#s <- rowSums(genus_table[,2:sample_number+1])
#genus_table$genus[genus_table$genus != "unidentified" & s<300] <- "others"

#order again the table
genus_table <-genus_table[order(genus_table$genus),]

#merge all the order with the same name
genus_table <- aggregate(genus_table[1:sample_number+1], by=list(genus_table$genus), FUN = "sum")
genusnames <- genus_table$Group.1
genus_table <- t(genus_table)
colnames(genus_table)=genusnames
genus_table=genus_table[-1,]
#save genus_table
write.table(genus_table,file = "results_moldITS2/genus_level.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the colsums of all samples
class(genus_table) <- "double"
x<-genus_table[,order(colSums(genus_table),decreasing=TRUE)]
taxa_list <-colnames(x)

x <- genus_table/rowSums(genus_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

df <- data.frame(Sample_ID=rownames(genus_table), 
                 biome=biomes,
                 genus_table,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE
)
#change the number to numeric
df[,2:ncol(df)]= apply(df[,2:ncol(df)], 2, function(x) as.numeric(as.character(x)))
#merge df with biome
df <-df[order(df$biome),]
df <- aggregate(df[2:ncol(df)], by=list(df$biome), FUN = "sum")
rownames(df)=df$Group.1
df=df[,-1]
write.table(df,file = "results_moldITS2/genus_level_biome.txt", sep = "\t", quote = FALSE, row.names = T)

#select only 26 taxa to show
selectedtaxa <- taxa_list[1:26]
new_df <- data.frame(Sample_ID=colnames(df), 
                     taxa=colnames(df),
                     t(df),
                     row.names = 1, 
                     check.names=FALSE,
                     stringsAsFactors = FALSE
)
new_df$taxa[!(new_df$taxa %in% selectedtaxa)] <- "others"
new_df <-new_df[order(new_df$taxa),]
new_df <- aggregate(new_df[2:ncol(new_df)], by=list(new_df$taxa), FUN = "sum")
rownames(new_df)=new_df$Group.1
new_df=new_df[,-1]
new_df <- t(new_df)

#count all classes
s <- colSums(new_df[,1:ncol(new_df)])
df_all <- data.frame(t(s))
row.names(df_all) = c("ALL")
new_df <- rbind(new_df,df_all)

#move other to the end
taxa_list <-c(selectedtaxa,"others")

x <-new_df[,taxa_list]

#move "__unidentified__" to the end
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-x[,taxa_list]

new_x <- new_x/rowSums(new_x)
df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),
                  Genera=rep(colnames(new_x)[i],dim(new_x)[1]),
                  Value=new_x[,i],
                  Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
#species_colours <- c("#0075DC","#F0A3FF", "#993F00","#2BCE48","#4C005C","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
colours <- c("#F0A3FF", "#0075DC", "#a1d0ed","#2BCE48","#401019","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#993F00","#4C005C","#FFFF00")
N<-length(taxa_list)
#colours[N]="#FFFF00" #yellow
colours[N]="#000000" #black
library(ggplot2)
p_bg<-ggplot(df,aes(Sample,Value,fill=Genera))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")  
p_bg<-p_bg+scale_fill_manual(values=colours[1:N])
p_bg<-p_bg+theme_bw()+ylab("Proportions")
p_bg<-p_bg+theme_bw(base_size=11)+xlab("Biomes")
p_bg<-p_bg+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p_bg<-p_bg+theme(axis.text.x=element_blank())
p_bg<- p_bg + theme(text=element_text(size=20))
#p_bg<- p_bg + theme(text=element_text(size=15))
#p_bg<-p_bg+ theme(axis.title.y = element_text(size = 15))
#p_bg<-p_bg+ theme(axis.text.y = element_text(size = 15))
#p_bg<-p_bg+ theme(legend.title = element_text(size = 15))
p_bg<-p_bg+ theme(legend.text = element_text(size = 20,face='italic'))
#p_bg<-p_bg+ theme(strip.text = element_text(size = 15)) 
pdf("results_moldITS2/genus_level_biome.pdf",height=6,width=21)
print(p_bg)
dev.off()


###########################
#Taxonomy at genus level ----
genus_table <-data.frame(OTU_ID = row.names(classification_table), 
                         genus = as.character(classification_table$genus),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
genus_table <-genus_table[order(genus_table$genus),]
#replace empty class with "unidentified"
genus_table$genus[genus_table$genus==""] <- "unidentified"
genus_table$genus[genus_table$genus=="Mastigobasidium"] <- "Leucosporidium"
genus_table$genus[genus_table$genus=="Schizoblastosporion"] <- "Nadsonia"

#replace unknown classes with "unidentified"
genus_table$genus[startsWith(genus_table$genus,"Branch0")] <- "unidentified"
genus_table$genus[grepl("Incertae_sedis",genus_table$genus)] <- "unidentified"
genus_table$genus[startsWith(genus_table$genus,"GS")] <- "unidentified"

#change the name of insignificant groups to others
#s <- rowSums(genus_table[,2:sample_number+1])
#genus_table$genus[genus_table$genus != "unidentified" & s<300] <- "others"

#order again the table
genus_table <-genus_table[order(genus_table$genus),]

#merge all the order with the same name
genus_table <- aggregate(genus_table[1:sample_number+1], by=list(genus_table$genus), FUN = "sum")
genusnames <- genus_table$Group.1
genus_table <- t(genus_table)
colnames(genus_table)=genusnames
genus_table=genus_table[-1,]
#save genus_table
write.table(genus_table,file = "results_moldITS2/genus_level.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the colsums of all samples
#class(genus_table) <- "double"
#x <- genus_table/rowSums(genus_table)
#x <- x[,order(colSums(x),decreasing=TRUE)]
#taxa_list <-colnames(x)
#select only 26 taxa to show
#selectedtaxa <- taxa_list[1:26]

#change the name of the taxa not in the selectedlist to others
new_df <- data.frame(Sample_ID=colnames(genus_table), 
                     taxa=colnames(genus_table),
                     t(genus_table),
                     row.names = 1, 
                     check.names=FALSE,
                     stringsAsFactors = FALSE
)
new_df$taxa[!(new_df$taxa %in% selectedtaxa)] <- "others"
new_df <-new_df[order(new_df$taxa),]
#change the number to numeric
new_df[,2:ncol(new_df)]= apply(new_df[,2:ncol(new_df)], 2, function(x) as.numeric(as.character(x)))
new_df <- aggregate(new_df[2:ncol(new_df)], by=list(new_df$taxa), FUN = "sum")
rownames(new_df)=new_df$Group.1
new_df=new_df[,-1]
new_df <- t(new_df)

#x <- new_df
#x <- x[,order(colSums(x),decreasing=TRUE)]

#taxa_list <- c(selectedtaxa,"others")
#new_x<-x[,taxa_list]

##remove "__unidentified__" 
#taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
#taxa_list <- c(taxa_list,"unidentified")
new_x<-new_df[,taxa_list]

new_x <- new_x/rowSums(new_x)
df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Genera=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
#colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
colours <- c("#F0A3FF", "#0075DC", "#a1d0ed","#2BCE48","#401019","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#993F00","#4C005C","#FFFF00")
N<-length(taxa_list)
#colours[N]="#FFFF00"#yellow
colours[N]="#000000" #black
library(ggplot2)
p_g<-ggplot(df,aes(Sample,Value,fill=Genera))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p_g<-p_g+scale_fill_manual(values=colours[1:N])
p_g<-p_g+theme_bw()+ylab("Proportions")
p_g<-p_g+theme_bw()+xlab("Samples")
p_g<-p_g+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p_g<-p_g+theme(axis.text.x=element_blank())
p_g<- p_g + theme(text=element_text(size=20))
#p_g<-p_g+ theme(axis.title.y = element_text(size = 15))
#p_g<-p_g+ theme(axis.text.y = element_text(size = 15))
#p_g<-p_g+ theme(legend.title = element_text(size = 15))
p_g<-p_g+ theme(legend.text = element_text(size = 20,face='italic')) 
#p_g<-p_g+ theme(strip.text = element_text(size = 15)) 
pdf("results_moldITS2/genus_level.pdf",height=6,width=21)
print(p_g)
dev.off()


####################### 
#combine figures at the genus level----
library(gridExtra)
library(tidyverse)
library(rstatix)
library(ggpubr)
gp_1 = ggplotGrob(p_g)
gp_2 = ggplotGrob(p_bg)
#gp_f$widths <- gp_g$widths
figure <- grid.arrange(gp_1 , gp_2)

#figure <- ggarrange(gp_1,gp_2,
#                    ncol = 1, nrow = 2)
pdf("results_moldITS2/genus_diversity.pdf",height=12,width=21)
print(figure)
dev.off()

###########################
#Taxonomy at species level ----
species_table <-data.frame(OTU_ID = row.names(classification_table), 
                         species = as.character(classification_table$species),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
species_table <-species_table[order(species_table$species),]
#replace empty class with "unidentified"
species_table$species[species_table$species==""] <- "unidentified"
species_table$species[species_table$species=="Apiotrichum xylopini"] <- "Apiotrichum porosum"
species_table$species[species_table$species=="Solicoccozyma fuscescens"] <- "Solicoccozyma terrea"
species_table$species[species_table$species=="Cryptococcus ramirezgomezianus"] <- "Vanrija albida"
species_table$species[species_table$species=="Schizoblastosporium starkeyi-henricii"] <- "Nadsonia starkeyi-henricii"

#order again the table
species_table <-species_table[order(species_table$species),]

#merge all the order with the same name
species_table <- aggregate(species_table[1:sample_number+1], by=list(species_table$species), FUN = "sum")
speciesnames <- species_table$Group.1
species_table <- t(species_table)
colnames(species_table)=speciesnames
species_table=species_table[-1,]
#save genus_table
write.table(species_table,file = "results_moldITS2/species_level.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the colsums of all samples
class(species_table) <- "double"
x <- species_table/rowSums(species_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)
#select only 26 taxa to show
selectedtaxa <- taxa_list[1:26]

#change the name of the taxa not in the selectedlist to others
new_df <- data.frame(Sample_ID=colnames(species_table), 
                     taxa=colnames(species_table),
                     t(species_table),
                     row.names = 1, 
                     check.names=FALSE,
                     stringsAsFactors = FALSE
)
new_df$taxa[!(new_df$taxa %in% selectedtaxa)] <- "others"
new_df <-new_df[order(new_df$taxa),]
#change the number to numeric
new_df[,2:ncol(new_df)]= apply(new_df[,2:ncol(new_df)], 2, function(x) as.numeric(as.character(x)))
new_df <- aggregate(new_df[2:ncol(new_df)], by=list(new_df$taxa), FUN = "sum")
rownames(new_df)=new_df$Group.1
new_df=new_df[,-1]
new_df <- t(new_df)

x <- new_df/rowSums(new_df)
x <- x[,order(colSums(x),decreasing=TRUE)]

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]

#remove "__unidentified__" 
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-new_x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Species=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
#colours[N]="#FFFF00"#yellow
colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Species))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results_moldITS2/species_level.pdf",height=6,width=21)
print(p)
dev.off()

###########################
#species_level_biome.pdf----
species_table <-data.frame(OTU_ID = row.names(classification_table), 
                           species = as.character(classification_table$species),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
species_table <-species_table[order(species_table$species),]
#replace empty class with "unidentified"
species_table$species[species_table$species==""] <- "unidentified"
species_table$species[species_table$species=="Apiotrichum xylopini"] <- "Apiotrichum porosum"
species_table$species[species_table$species=="Solicoccozyma fuscescens"] <- "Solicoccozyma terrea"
species_table$species[species_table$species=="Cryptococcus ramirezgomezianus"] <- "Vanrija albida"
species_table$species[species_table$species=="Schizoblastosporion starkeyi-henricii"] <- "Nadsonia starkeyi-henricii"

#order again the table
species_table <-species_table[order(species_table$species),]

#merge all the order with the same name
species_table <- aggregate(species_table[1:sample_number+1], by=list(species_table$species), FUN = "sum")
speciesnames <- species_table$Group.1
species_table <- t(species_table)
colnames(species_table)=speciesnames
species_table=species_table[-1,]
#save genus_table
write.table(species_table,file = "results_moldITS2/species_level.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the colsums of all samples
class(species_table) <- "double"
x<-species_table[,order(colSums(species_table),decreasing=TRUE)]
taxa_list <-colnames(x)

x <- species_table/rowSums(species_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

df <- data.frame(Sample_ID=rownames(species_table), 
                 biome=biomes,
                 species_table,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE
)
#change the number to numeric
df[,2:ncol(df)]= apply(df[,2:ncol(df)], 2, function(x) as.numeric(as.character(x)))
#merge df with biome
df <-df[order(df$biome),]
df <- aggregate(df[2:ncol(df)], by=list(df$biome), FUN = "sum")
rownames(df)=df$Group.1
df=df[,-1]
write.table(df,file = "results_moldITS2/species_level_biome.txt", sep = "\t", quote = FALSE, row.names = T)

#select only 26 taxa to show
selectedtaxa <- taxa_list[1:26]
new_df <- data.frame(Sample_ID=colnames(df), 
                     taxa=colnames(df),
                     t(df),
                     row.names = 1, 
                     check.names=FALSE,
                     stringsAsFactors = FALSE
)
new_df$taxa[!(new_df$taxa %in% selectedtaxa)] <- "others"
new_df <-new_df[order(new_df$taxa),]
new_df <- aggregate(new_df[2:ncol(new_df)], by=list(new_df$taxa), FUN = "sum")
rownames(new_df)=new_df$Group.1
new_df=new_df[,-1]
new_df <- t(new_df)

#count all classes
s <- colSums(new_df[,1:ncol(new_df)])
df_all <- data.frame(t(s))
row.names(df_all) = c("ALL")
colnames(df_all) <- colnames(new_df)
new_df <- rbind(new_df,df_all)

taxa_list <-c(selectedtaxa,"others")

x <- new_df/rowSums(new_df)
x <-x[,taxa_list]
#remove others

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]

#remove "__unidentified__" 
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-new_x[,taxa_list]
df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),
                  Species=rep(colnames(new_x)[i],dim(new_x)[1]),
                  Value=new_x[,i],
                  Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#0075DC","#F0A3FF", "#993F00","#2BCE48","#4C005C","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
#colours[N]="#FFFF00" #yellow
colours[N]="#000000" #black
library(ggplot2)
p_s<-ggplot(df,aes(Sample,Value,fill=Species))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p_s<-p_s+scale_fill_manual(values=colours[1:N])
p_s<-p_s+theme_bw()+ylab("Proportions")
p_s<-p_s+theme_bw()+xlab("Biotopes")
p_s<-p_s+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p_s<-p_s+theme(axis.text.x=element_blank())
p_s<- p_s + theme(text=element_text(size=20))
pdf("results_moldITS2/species_level_biome.pdf",height=6,width=21)
print(p_s)
dev.off()

####################### 
#combine figures ----
library(gridExtra)
library(tidyverse)
library(rstatix)
library(ggpubr)
gp_1 = ggplotGrob(p_f)
gp_2 = ggplotGrob(p_o)
gp_3 = ggplotGrob(p_c)

#gp_f$widths <- gp_g$widths
figure <- grid.arrange(gp_1 , gp_2, gp_3)

#figure <- ggarrange(gp_1,gp_2,gp_3,
#                    ncol = 1, nrow = 3)
pdf("results_moldITS2/family_order_class_diversity.pdf",height=12,width=21)
print(figure)
dev.off()

##############################################
#Richness ----
richness <-colSums(otu_table > 0)
#biome richness
df <- data.frame(Sample_ID=colnames(otu_table), 
                 biome=biomes,
                 t(otu_table),
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE)
#change the number to numeric
df[,2:ncol(df)]= apply(df[,2:ncol(df)], 2, function(x) as.numeric(as.character(x)))
#merge all the samples of the same biome
df <- aggregate(df[2:ncol(df)], by=list(df$biome), FUN = "sum")
rownames <- df$Group.1
df <- df[2:ncol(df)]
rownames(df) <-rownames
biome_richness <- rowSums(df>0)
biome_richness_proportion <- biome_richness*100/sum(biome_richness)
biome_abundance <- rowSums(df)
biome_abundance_proportion <- biome_abundance*100/sum(biome_abundance)
df <- data.frame(Sample_ID=rownames, 
                 Richness=biome_richness,
                 Abundance=biome_abundance,
                 Richness_proportion = biome_richness_proportion,
                 Abundance_proportion = biome_abundance_proportion,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE)

cor(df$Richness,df$Abundance)
df <- df[order(df$Abundance),]
#write to file
write.table(df,file = "results/richness.txt", sep = "\t", quote = FALSE, row.names = T)

proportionlist <- c(df$Richness,df$Abundance)

#proportionlist <- c(df$Richness_proportion,df$Abundance_proportion)

df2 <- data.frame(Diversity=rep(c("Richness", "Abundance"), each=dim(df)[1]),
                  Biotope_type=rep(rownames(df),2),
                  Number=proportionlist)

library(ggplot2)
p<-ggplot(data=df2, aes(x=Biotope_type, y=Number, fill=Diversity)) +
geom_bar(stat="identity", position=position_dodge()) +  theme(text=element_text(size=20))
pdf("results_moldITS2/richness_abundance.pdf",height=6,width=21)
print(p)
dev.off()


#########################################
#Functional diversity  richness and abundance by cbs ----
#guild table
#table <-data.frame(OTU_ID = row.names(classification_table), 
#                   guild = as.character(guild_table$Guild),
#                   otu_table,
#                   row.names = 1, 
#                   check.names=FALSE,
#                   stringsAsFactors = FALSE)
#guild table

#guild_table <-read.delim("dnabarcoder/soil.unite.yeastCBSITS2_BLAST.classification.filtered.guilds.txt",row.names=1,check.names=FALSE)
guild_table <-read.delim("dnabarcoder/soil_ITS2.moldITS2_BLAST.classification.better.guilds.txt",row.names=1,check.names=FALSE)
guild_table<-guild_table[rownames(classification_table), ]
table <-data.frame(OTU_ID = row.names(guild_table), 
                   guild = as.character(guild_table$Guild),
                   otu_table,
                   row.names = 1, 
                   check.names=FALSE,
                   stringsAsFactors = FALSE)

table <-table[order(table$guild),]
uniqueguilds <-unique(table$guild)
unknown <- table[table$guild=="",]
dim(unknown)

guild_list <-c("Animal- and mycoparasites", "Animal Pathogen", "Ectomycorrhizal fungi", "Plant Pathogen",
               "Saprotrophs", "others", "unassigned")

table$guild[grepl("animal pathogen",tolower(table$guild))] <- "Animal Pathogen"
table$guild[grepl("animal_pathogen",tolower(table$guild))] <- "Animal Pathogen"
#plant pathogen  and sap
plantpathsaprotroph <- table[grepl("plant pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(plantpathsaprotroph)
plantpathsaprotroph <- table[grepl("plant_pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(plantpathsaprotroph)
table$guild[grepl("plant pathogen",tolower(table$guild))] <- "Plant Pathogen"
table$guild[grepl("plant_pathogen",tolower(table$guild))] <- "Plant Pathogen"
table$guild[grepl("mycoparasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("fungal parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("fungal_parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("lichen parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("lichen_parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("animal parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("animal_parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
#table$guild[grepl("epiphyte",tolower(table$guild))] <- "Epiphytes"
#table$guild[grepl("arbuscular mycorrhizal",tolower(table$guild))] <- "AM fungi"
#table$guild[grepl("arbuscular_mycorrhizal",tolower(table$guild))] <- "AM fungi"
#table$guild[grepl("plant pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild))] <- "Ecto-Saprotrophs"
ectosaprotroph <- table[grepl("ectomycorrhiza",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(ectosaprotroph)
table$guild[grepl("ectomycorrhizal",tolower(table$guild))] <- "Ectomycorrhizal fungi"
table$guild[grepl("saprotroph",tolower(table$guild))] <- "Saprotrophs"

table$guild[table$guild==""] <- "unassigned"
table$guild[table$guild=="-"] <- "unassigned"
table$guild[!(table$guild %in% guild_list)] <-"others"

df <- data.frame(Sample_ID=rownames(table), 
                 Guild=table$guild,
                 Abundance=rowSums(table[2:ncol(table)]),
                 Richness=1,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE)



#merge all the class with the same name
df <- aggregate(df[2:ncol(df)], by=list(df$Guild), FUN = "sum")
guildnames <- df$Group.1
rownames(df) <-guildnames
df <-df[2:ncol(df)]
richness_proportion <- df$Richness/sum(df$Richness)
abundance_proportion <- df$Abundance/sum(df$Abundance)
df_better <- data.frame(Sample_ID=guildnames, 
                 Richness=df$Richness,
                 Abundance=df$Abundance,
                 Richness_proportion =richness_proportion,
                 Abundance_proportion = abundance_proportion,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE)

df_better<-df_better[guild_list, ] #order df by guild_list

cor(df_better$Richness,df_better$Abundance)
#write to file
write.table(df_better,file = "results_moldITS2/CBS_richness_guild.txt", sep = "\t", quote = FALSE, row.names = T)

#proportionlist <- c(df$Richness_proportion,df$Abundance_proportion)
proportionlist <- c(df_better$Richness,df_better$Abundance)

#df2 <- data.frame(Diversity=rep(c("Richness", "Abundance"), each=dim(df)[1]),
#                  Functionality=rep(rownames(df),2),
#                  Number=proportionlist)

#library(ggplot2)
#p<-ggplot(data=df2, aes(x=Functionality, y=Number, fill=Diversity)) +
#  geom_bar(stat="identity",position=position_dodge()) +  theme(text=element_text(size=20))
#pdf("results/functionality_richness_abundance.pdf",height=6,width=21)
#print(p)
#dev.off()

op <- par(mar = c(12,4,4,2) + 0.2)
p_richness<-barplot(df_better$Richness, names.arg = guild_list, las=2, cex.lab=1, ylab = "", col= "#0075DC" )
par(op)
pdf("results_moldITS2/CBS_richness.pdf",height=6,width=21)
print(p_richness)
dev.off()
op <- par(mar = c(12,4,4,2) + 0.2)
p_abundance<-barplot(df_better$Abundance, names.arg = guild_list, las=2, cex.lab=1, ylab = "", col="#2BCE48")
par(op)
pdf("results_moldITS2/CBS_abundance.pdf",height=6,width=21)
print(p_abundance)
dev.off()

#########################################
#UNITE Functional diversity  richness and abundance by unite ----
#guild table
#table <-data.frame(OTU_ID = row.names(classification_table), 
#                   guild = as.character(guild_table$Guild),
#                   otu_table,
#                   row.names = 1, 
#                   check.names=FALSE,
#                   stringsAsFactors = FALSE)
#guild table

#guild_table <-read.delim("dnabarcoder/soil.unite.yeastCBSITS2_BLAST.classification.filtered.guilds.txt",row.names=1,check.names=FALSE)
guild_table <-read.delim("dnabarcoder/soil.unite.classification.guilds.txt",row.names=1,check.names=FALSE)
samplenames <-intersect(rownames(guild_table),rownames(allotu_table))
unite_otu_table <-allotu_table[samplenames, ]
guild_table <-guild_table[samplenames, ]
#guild_table<-guild_table[rownames(classification_table), ]
table <-data.frame(OTU_ID = row.names(guild_table), 
                   guild = as.character(guild_table$Guild),
                   unite_otu_table,
                   row.names = 1, 
                   check.names=FALSE,
                   stringsAsFactors = FALSE)

table <-table[order(table$guild),]
uniqueguilds <-unique(table$guild)
unknown <- table[table$guild=="",]
dim(unknown)

guild_list <-c("Animal- and mycoparasites", "Animal Pathogen", "Ectomycorrhizal fungi", "Plant Pathogen",
               "Saprotrophs", "others", "unassigned")

table$guild[grepl("animal pathogen",tolower(table$guild))] <- "Animal Pathogen"
table$guild[grepl("animal_pathogen",tolower(table$guild))] <- "Animal Pathogen"
#plant pathogen  and sap
plantpathsaprotroph <- table[grepl("plant pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(plantpathsaprotroph)
plantpathsaprotroph <- table[grepl("plant_pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(plantpathsaprotroph)
table$guild[grepl("plant pathogen",tolower(table$guild))] <- "Plant Pathogen"
table$guild[grepl("plant_pathogen",tolower(table$guild))] <- "Plant Pathogen"
table$guild[grepl("mycoparasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("fungal parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("fungal_parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("lichen parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("lichen_parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("animal parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("animal_parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
#table$guild[grepl("epiphyte",tolower(table$guild))] <- "Epiphytes"
#table$guild[grepl("arbuscular mycorrhizal",tolower(table$guild))] <- "AM fungi"
#table$guild[grepl("arbuscular_mycorrhizal",tolower(table$guild))] <- "AM fungi"
#table$guild[grepl("plant pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild))] <- "Ecto-Saprotrophs"
ectosaprotroph <- table[grepl("ectomycorrhiza",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(ectosaprotroph)
table$guild[grepl("ectomycorrhizal",tolower(table$guild))] <- "Ectomycorrhizal fungi"
table$guild[grepl("saprotroph",tolower(table$guild))] <- "Saprotrophs"

table$guild[table$guild==""] <- "unassigned"
table$guild[table$guild=="-"] <- "unassigned"
table$guild[!(table$guild %in% guild_list)] <-"others"

df <- data.frame(Sample_ID=rownames(table), 
                 Guild=table$guild,
                 Abundance=rowSums(table[2:ncol(table)]),
                 Richness=1,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE)



#merge all the class with the same name
df <- aggregate(df[2:ncol(df)], by=list(df$Guild), FUN = "sum")
guildnames <- df$Group.1
rownames(df) <-guildnames
df <-df[2:ncol(df)]
richness_proportion <- df$Richness/sum(df$Richness)
abundance_proportion <- df$Abundance/sum(df$Abundance)
df_better <- data.frame(Sample_ID=guildnames, 
                        Richness=df$Richness,
                        Abundance=df$Abundance,
                        Richness_proportion =richness_proportion,
                        Abundance_proportion = abundance_proportion,
                        row.names = 1, 
                        check.names=FALSE,
                        stringsAsFactors = FALSE)

df_better<-df_better[guild_list, ] #order df by guild_list

cor(df_better$Richness,df_better$Abundance)
#write to file
write.table(df_better,file = "results_moldITS2/UNITE_richness_guild.txt", sep = "\t", quote = FALSE, row.names = T)

#proportionlist <- c(df$Richness_proportion,df$Abundance_proportion)
proportionlist <- c(df_better$Richness,df_better$Abundance)

#df2 <- data.frame(Diversity=rep(c("Richness", "Abundance"), each=dim(df)[1]),
#                  Functionality=rep(rownames(df),2),
#                  Number=proportionlist)

#library(ggplot2)
#p<-ggplot(data=df2, aes(x=Functionality, y=Number, fill=Diversity)) +
#  geom_bar(stat="identity",position=position_dodge()) +  theme(text=element_text(size=20))
#pdf("results/functionality_richness_abundance.pdf",height=6,width=21)
#print(p)
#dev.off()

op <- par(mar = c(12,4,4,2) + 0.2)
p_richness<-barplot(df_better$Richness, names.arg = guild_list, las=2, cex.lab=1, ylab = "", col= "#0075DC" )
par(op)
pdf("results_moldITS2/UNITE_richness.pdf",height=6,width=21)
print(p_richness)
dev.off()
op <- par(mar = c(12,4,4,2) + 0.2)
p_abundance<-barplot(df_better$Abundance, names.arg = guild_list, las=2, cex.lab=1, ylab = "", col="#2BCE48")
par(op)
pdf("results_moldITS2/UNITE_abundance.pdf",height=6,width=21)
print(p_abundance)
dev.off()

#########################################
#Functional diversity  richness and abundance only by CBS classification----
classification_table_onlybyCBS <-read.delim("dnabarcoder/soil_ITS2.moldITS2_BLAST.classification.only",row.names=1,check.names=FALSE)
otu_table_onlybyCBS <- allotu_table[rownames(classification_table_onlybyCBS), ]

guild_table_onlybyCBS <-read.delim("dnabarcoder/soil_ITS2.moldITS2_BLAST.classification.only.guilds.txt",row.names=1,check.names=FALSE)
#guild_table_onlybyCBS<-guild_table[rownames(classification_table_onlybyCBS), ]
table <-data.frame(OTU_ID = row.names(guild_table_onlybyCBS), 
                   guild = as.character(guild_table_onlybyCBS$Guild),
                   otu_table_onlybyCBS,
                   row.names = 1, 
                   check.names=FALSE,
                   stringsAsFactors = FALSE)

table <-table[order(table$guild),]
uniqueguilds <-unique(table$guild)
unknown <- table[table$guild=="",]
dim(unknown)

guild_list <-c("Animal- and mycoparasites", "Animal Pathogen", "Ectomycorrhizal fungi", "Plant Pathogen",
               "Saprotrophs", "others", "unassigned")

table$guild[grepl("animal pathogen",tolower(table$guild))] <- "Animal Pathogen"
table$guild[grepl("animal_pathogen",tolower(table$guild))] <- "Animal Pathogen"
#plant pathogen  and sap
plantpathsaprotroph <- table[grepl("plant pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(plantpathsaprotroph)
plantpathsaprotroph <- table[grepl("plant_pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(plantpathsaprotroph)
table$guild[grepl("plant pathogen",tolower(table$guild))] <- "Plant Pathogen"
table$guild[grepl("plant_pathogen",tolower(table$guild))] <- "Plant Pathogen"
table$guild[grepl("mycoparasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("fungal parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("fungal_parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("lichen parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("lichen_parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("animal parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("animal_parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
#table$guild[grepl("epiphyte",tolower(table$guild))] <- "Epiphytes"
#table$guild[grepl("arbuscular mycorrhizal",tolower(table$guild))] <- "AM fungi"
#table$guild[grepl("arbuscular_mycorrhizal",tolower(table$guild))] <- "AM fungi"
#table$guild[grepl("plant pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild))] <- "Ecto-Saprotrophs"
ectosaprotroph <- table[grepl("ectomycorrhiza",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(ectosaprotroph)
table$guild[grepl("ectomycorrhizal",tolower(table$guild))] <- "Ectomycorrhizal fungi"
table$guild[grepl("saprotroph",tolower(table$guild))] <- "Saprotrophs"

table$guild[table$guild==""] <- "unassigned"
table$guild[table$guild=="-"] <- "unassigned"
table$guild[!(table$guild %in% guild_list)] <-"others"

df <- data.frame(Sample_ID=rownames(table), 
                 Guild=table$guild,
                 Abundance=rowSums(table[2:ncol(table)]),
                 Richness=1,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE)



#merge all the class with the same name
df <- aggregate(df[2:ncol(df)], by=list(df$Guild), FUN = "sum")
guildnames <- df$Group.1
rownames(df) <-guildnames
df <-df[2:ncol(df)]
richness_proportion <- df$Richness/sum(df$Richness)
abundance_proportion <- df$Abundance/sum(df$Abundance)
df_only <- data.frame(Sample_ID=guildnames, 
                 Richness=df$Richness,
                 Abundance=df$Abundance,
                 Richness_proportion =richness_proportion,
                 Abundance_proportion = abundance_proportion,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE)

df_only<-df_only[guild_list, ] #order df by guild_list

cor(df_only$Richness,df_only$Abundance)
#write to file
write.table(df_only,file = "results_moldITS2/CBS_only_richness_guild.txt", sep = "\t", quote = FALSE, row.names = T)

#proportionlist <- c(df$Richness_proportion,df$Abundance_proportion)
proportionlist <- c(df_only$Richness,df_only$Abundance)

#df2 <- data.frame(Diversity=rep(c("Richness", "Abundance"), each=dim(df)[1]),
#                  Functionality=rep(rownames(df),2),
#                  Number=proportionlist)

#library(ggplot2)
#p<-ggplot(data=df2, aes(x=Functionality, y=Number, fill=Diversity)) +
#  geom_bar(stat="identity",position=position_dodge()) +  theme(text=element_text(size=20))
#pdf("results/functionality_richness_abundance.pdf",height=6,width=21)
#print(p)
#dev.off()

op <- par(mar = c(12,4,4,2) + 0.2)
p_richness<-barplot(df_only$Richness, names.arg = guild_list, las=2, cex.lab=1, ylab = "", col= "#0075DC" )
par(op)
pdf("results/CBS_only_richness.pdf",height=6,width=21)
print(p_richness)
dev.off()
op <- par(mar = c(12,4,4,2) + 0.2)
p_abundance<-barplot(df_only$Abundance, names.arg = guild_list, las=2, cex.lab=1, ylab = "", col="#2BCE48")
par(op)
pdf("results_moldITS2/CBS_only_abundance.pdf",height=6,width=21)
print(p_abundance)
dev.off()
#########################################
#Combining figures for richness and abundance----
op <- par(mar = c(12,4,4,2) + 0.2)
values <- t(cbind(df_better$Richness, df_only$Richness))
colors <-c("#2BCE48","#993F00")
Type <- c("UNITE&CBS","CBS only")
#p_richness <- barplot(counts, names.arg = guild_list, main="",
#        xlab="", col=colors, beside=TRUE)
p_richness<-barplot(values, names.arg = guild_list, las=2, cex.lab=1, ylab = "Number of OTUs", col = colors, besid=TRUE)
par(op)
legend("topleft", Type, cex = 1, fill = colors)
pdf("results_moldITS2/functionality_richness.pdf",height=6,width=21)
print(p_richness)
dev.off()

op <- par(mar = c(12,4,4,2) + 0.2)
values <- t(cbind(df_better$Abundance, df_only$Abundance))
colors <-c("#2BCE48","#993F00")
Type <- c("UNITE&CBS","CBS only")
#p_richness <- barplot(counts, names.arg = guild_list, main="",
#        xlab="", col=colors, beside=TRUE)
p_abundance<-barplot(values, names.arg = guild_list, las=2, cex.lab=1, col = colors, besid=TRUE)
title(ylab = "Number of OTUs")
par(op)
legend("topleft", Type, cex = 1, fill = colors)
pdf("results_moldITS2/functionalityabundance.pdf",height=6,width=21)
print(p_abundance)
dev.off()


#########################################
#Comparing classification----
#take data from .compared.species
#species
#colors = c("green","red","blue")

colors <- c("#2BCE48","#993F00", "#F0A3FF","#0075DC")

c_s<-c(122,122,30,195,195,30,1364,1146)
c_g<-c(2007,2007,683,1532,1532,683,7641,1686)
c_f<-c(4123,4123,735,2306,2306,735,10708,1649)
c_o<-c(8485,8485,538,1143,1143,538,12726,742)
c_c<-c(22298,22298,1084,957,957,1084,208,112)

l_s <- c("Species_by_UNITE","Species_by_CBS")
l_g <- c("Genera_by_UNITE","Genera_by_CBS")
l_f <- c("Families_by_UNITE","Families_by_CBS")
l_o <- c("Orders_by_UNITE","Orders_by_CBS")
l_c <- c("Classes_by_UNITE","Classes_by_CBS")


#l_s <- c("UNITE_species","CBS_species")
#l_g <- c("UNITE_genus","CBS_genus")
#l_f <- c("UNITE_family","CBS_family")
#l_o <- c("UNITE_order","CBS_order")
#l_c <- c("UNITE_class","CBS_class")

labels<-cbind(l_s,l_g,l_f,l_o,l_c)


Type <- c("Same","Diff. with higher score","Diff. with lower score","Only by")
v_s <- matrix(c_s, nrow = 4, ncol = 2, byrow = TRUE)
v_g <- matrix(c_g, nrow = 4, ncol = 2, byrow = TRUE)
v_f <- matrix(c_f, nrow = 4, ncol = 2, byrow = TRUE)
v_o <- matrix(c_o, nrow = 4, ncol = 2, byrow = TRUE)
v_c <- matrix(c_c, nrow = 4, ncol = 2, byrow = TRUE)
values<-cbind(v_s,v_g,v_f,v_o,v_c)

op <- par(mar = c(12,4,4,2) + 2)
p<-barplot(values, names.arg = labels, las=2, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, col = colors)
#p<-barplot(values, names.arg = labels, las=2, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, ylab = "Number of OTUs", col = colors
#p<-barplot(values, names.arg = labels, las=2, cex.lab=1, ylab = "Number of OTUs", col = colors)
title(ylab = "Number of OTUs", cex.lab = 1.5,line = 4.5)
par(op)
legend("topleft", Type, cex = 1.5, fill = colors)



#Type <- c("Same","Diff. with higher score","Diff. with lower score","Only by")
#values <- matrix(c_s, nrow = 4, ncol = 2, byrow = TRUE)
#barplot(Values, names.arg = labels, xlab = "At the species level", ylab = "Number of identifications", col = colors)
#legend("topleft", Type, cex = 1.3, fill = colors

pdf("results_moldITS2/UNITE_CBS_comparison.pdf",height=6,width=21)
print(p)
dev.off()





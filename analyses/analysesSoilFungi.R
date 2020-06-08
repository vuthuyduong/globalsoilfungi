#setwd("/home/duong/Data/Metagenomics/soil_samplesUNITE/sciencepaper_otus/itsx2")
setwd("C:/Users/Duong Vu/Documents/CBSPapers/SoilSamples-2020/data/itsx2/metadata")
#vegan
#install.packages("vegan")
library("vegan")
library(grid)



#load otu table
otu_table <-read.delim("soil_ITS2.table",row.names=1,check.names=FALSE)
guild_table <-read.delim("soil_ITS2.guild",row.names=1,check.names=FALSE)
#meta_table <-read.delim("/home/duong/Data/Metagenomics/soil_samplesUNITE/sciencepaper_otus/soil_ITS2.metadata",row.names=1,check.names=FALSE)
meta_table <-read.delim("C:/Users/Duong Vu/Documents/CBSPapers/SoilSamples-2020/data/itsx2/metadata/soil_ITS2.metadata",row.names=1,check.names=FALSE)

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
#write to file
write.table(df,file = "results/richness.txt", sep = "\t", quote = FALSE, row.names = T)

#########################################
# Heatmap with labels ----

#Get a heatmap comapring the samples
pdf("results/heatmap_withlabels.pdf", width = 7, height = 6)

ion$heatmap(scale(t(otu_table[r_sum > 100, 1:sample_number])),
            #row_labels = colnames(otu_table)[1:sample_number], 
            row_labels = labels, 
            row_margin = 13, row_distance = "spearman")
dev.off()
#########################################
# Heatmap with colorbar ----
#library(randomcoloR)
#n <- 30
#heatmap
source("https://tvpham.github.io/ion.r")
n <-length(unique(biomes))
color_vec <- RColorBrewer::brewer.pal(n, "Set3")
color_vec[2]="#000000" #change from yellow to black

aa <- factor(biomes)
color_bar <- color_vec[aa]

pdf("results/heatmap_withcolorbar.pdf", width = 12, height = 7)

ion$heatmap(otu_table[r_sum > 100, 1:sample_number],
            col_color_bar = list("Biome" = color_bar),
            z_transform = "row",
            row_margin = 25,
            lwid = c(1, 6),
            #row_labels = colnames(otu_table)[1:sample_number], 
            #col_labels = labels, 
            #col_margin = 13, 
            col_distance = "spearman")
par(mar = c(0, 0, 0, 0), fig = c(0.75, 1, 0, 0.72), new = TRUE)
plot.new()
legend("topleft", levels(aa), col = color_vec,
       pch=15, pt.cex=1.5, bty = "n")
dev.off()

#########################################
#Plant pathogen Heatmap with colorbar ----
#library(randomcoloR)
#n <- 30
#heatmap
source("https://tvpham.github.io/ion.r")

table <-data.frame(OTU_ID = row.names(guild_table), 
                          otu_table,row.names = 1, 
                          Guild=guild_table$Guild,
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
table <- subset(table, grepl("plant pathogen",tolower(table$Guild)))

n <-length(unique(biomes))
color_vec <- RColorBrewer::brewer.pal(n, "Set3")
color_vec[2]="#000000" #change from yellow to black

aa <- factor(biomes)
color_bar <- color_vec[aa]

pdf("results/heatmap_withcolorbar_plantpathogen.pdf", width = 12, height = 7)

ion$heatmap(table[1:sample_number],
            col_color_bar = list("Biome" = color_bar),
            z_transform = "row",
            row_margin = 25,
            lwid = c(1, 6),
            #row_labels = colnames(otu_table)[1:sample_number], 
            #col_labels = labels, 
            #col_margin = 13, 
            col_distance = "spearman")
par(mar = c(0, 0, 0, 0), fig = c(0.75, 1, 0, 0.72), new = TRUE)
plot.new()
legend("topleft", levels(aa), col = color_vec,
       pch=15, pt.cex=1.5, bty = "n")
dev.off()

#########################################
#Hclust based on the labels (biome country) ----
#merge samples based on biomes and countries
tr <-data.frame(sample_id=colnames(otu_table),
               label=labels,
               t(otu_table),
               row.names = 1, 
               check.names=FALSE,
               stringsAsFactors = FALSE)
tr <- tr[order(tr$label),]
#merge all rows with the same labels
tr <-aggregate(tr[1:otu_number+1], by=list(tr$label), FUN = "sum")
#write to file
write.table(tr,file = "results/biome_country.txt", sep = "\t", quote = FALSE, row.names = T)
# Pairwise correlation between labels
m<-data.matrix(tr[1:otu_number+1])
#compute distance
dd <- 1-cor(t(m), use = "pairwise.complete.obs", method = "pearson")
#hclust
hc <- hclust(as.dist(dd),method = "ward.D2")
hc$labels <-tr$Group.1
# Convert hclust into a dendrogram and plot
hcd <- as.dendrogram(hc)
#color the labels based on biomes
n <-length(unique(biomes))
color_vec <- RColorBrewer::brewer.pal(n, "Set3")
color_vec[2]="#000000" #change from yellow to black

uniquebiomes <-unique(biomes)
uniquelabels <-hc$labels
color_labels<-rep(color_vec[1],length(uniquelabels))
for (i in 1:length(uniquelabels)){
  for (j in 1:length(uniquebiomes)){
    if (startsWith(uniquelabels[i],uniquebiomes[j])){
      color_labels[i]<-color_vec[j]
      break
    }
  }
}

#plot
pdf("results/cluster.pdf", width = 12, height = 7)
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                cex = 0.7, col = "blue")

# reduced label size
par(cex=1, mar=c(20, 8, 4, 1))
#par(cex=1, mar=c(5, 8, 4, 1))

#install.packages('dendextend')
library(dendextend)
labels_colors(hcd) <- color_labels[order.dendrogram(hcd)]

plot(hcd, ylab = "Height", nodePar = nodePar,edgePar = list(col = 2:3, lwd = 2:1))

dev.off()

#########################################
#Linear regression----
cor(meta_table$Elevation, meta_table$pH) #0.057
cor(meta_table$Elevation, meta_table$N) #0.168
cor(meta_table$Elevation, meta_table$C) #0.128

#########################################
#Mantel test----
#compute geometric distances between samples
library(geosphere)
#distHaversine()
#distMeeus()
#distRhumb()
#distVincentyEllipsoid()
#distVincentySphere()
metadata=data.frame(sample_id=rownames(meta_table),
                    lon=meta_table$Longtitude,
                    lat=meta_table$Latitude,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
metadata <- subset(metadata, metadata$lon !="" & metadata$lat !="")
samples <- rownames(metadata)
gem_dd <- matrix(nrow = length(samples), ncol = length(samples))
for (i in 1:length(samples)){
  for (j in 1:length(samples)){
    gem_dd[i,j] =distm (c(metadata$lon[i], metadata$lat[i]), 
                       c(metadata$lon[j], metadata$lat[j]), 
                       fun = distHaversine)
  }
}

#Compute distance between the samples based on env parameters
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C=meta_table$C,
                    N=meta_table$N,
                    Elevation=meta_table$Elevation,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
metadata <- metadata[samples,]
env_m<-data.matrix(metadata)
#env_dd <-  as.matrix(dist(env_m, method = "euclidean"))
env_dd <- 1-cor(t(env_m), use = "pairwise.complete.obs", method = "pearson")
#compute Mantel test between gem_dd and env_dd
library(geosphere)
mantel1 <- mantel(gem_dd, env_dd, method = "spearman", permutations = 9999, na.rm = TRUE)

#Compute distance between the samples based on otu_numbers
#select only the samples with longtitude and latitude available
metadata <-otu_table[,samples]
m<-data.matrix(metadata)
#compute distance
#dd <- as.matrix(vegdist(t(m), method = "bray"))
dd <- 1-cor(m, use = "pairwise.complete.obs", method = "pearson")

#compute Mantel test between gem_dd and env_dd
mantel2 <- mantel(dd, env_dd, method = "spearman", permutations = 9999, na.rm = TRUE)

#########################################
#Mantel test plant pathogen----
table <-data.frame(OTU_ID = row.names(guild_table), 
                   otu_table,row.names = 1, 
                   Guild=guild_table$Guild,
                   check.names=FALSE,
                   stringsAsFactors = FALSE)
table <- subset(table, grepl("plant pathogen",tolower(table$Guild)))

#compute geometric distances between samples
library(geosphere)
#distHaversine()
#distMeeus()
#distRhumb()
#distVincentyEllipsoid()
#distVincentySphere()
metadata=data.frame(sample_id=rownames(meta_table),
                    lon=meta_table$Longtitude,
                    lat=meta_table$Latitude,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
metadata <- subset(metadata, metadata$lon !="" & metadata$lat !="")
samples <- rownames(metadata)
gem_dd <- matrix(nrow = length(samples), ncol = length(samples))
for (i in 1:length(samples)){
  for (j in 1:length(samples)){
    gem_dd[i,j] =distm (c(metadata$lon[i], metadata$lat[i]), 
                        c(metadata$lon[j], metadata$lat[j]), 
                        fun = distHaversine)
  }
}

#Compute distance between the samples based on env parameters
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C=meta_table$C,
                    N=meta_table$N,
                    Elevation=meta_table$Elevation,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
metadata <- metadata[samples,]
env_m<-data.matrix(metadata)
#env_dd <-  as.matrix(dist(env_m, method = "euclidean"))
env_dd <- 1-cor(t(env_m), use = "pairwise.complete.obs", method = "pearson")
#compute Mantel test between gem_dd and env_dd
library(geosphere)
mantel1 <- mantel(gem_dd, env_dd, method = "spearman", permutations = 9999, na.rm = TRUE)

#Compute distance between the samples based on otu_numbers
#select only the samples with longtitude and latitude available
metadata <-table[,samples]
m<-data.matrix(metadata)
#compute distance
#dd <- as.matrix(vegdist(t(m), method = "bray"))
plantpathogen_dd <- 1-cor(m, use = "pairwise.complete.obs", method = "pearson")

#compute Mantel test between gem_dd and env_dd
mantel2 <- mantel(plantpathogen_dd, env_dd, method = "spearman", permutations = 9999, na.rm = TRUE)

#########################################
#NMDS based on countries and biomes ----
subotus <- otu_table

#Load onl the otus that their sum is greater > 100)
#subotus <- otu_table[b > 100, 1:sample_number]
t_subotus <- t(subotus)
#Get MDS stats
sol<-metaMDS(t_subotus,distance = "bray", k = 2, trymax = 50)
#Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2], Biome=as.factor(biomes), Country=as.factor(countries))
#Draw boudaries
#Get spread of points based on seasons
plot.new()
ord<-ordiellipse(sol, as.factor(biomes) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()

#Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#Generate ellipse points
df_ell <- data.frame()
for(g in levels(NMDS$Biome)){
  if(g!="" && (g %in% names(ord))){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Biome==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                  ,Biome=g))
  }
}

head(df_ell)

#Generate mean values from NMDS plot grouped on biomes
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Biome),mean)
NMDS.mean
#plot
locationnumber <-length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot(data=NMDS,aes(x,y,color=Biome))
#p<-ggplot(data=NMDS,aes(x,y))
p<-p+ annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)
p<-p+ geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)
p<-p+geom_point(aes(shape=Country))+scale_shape_manual(values=shape_values)+theme_bw() 
pdf("results/NMDS.pdf", width = 12, height = 9)
print(p)
dev.off()
############################
#Correlation----
abund <-colSums(otu_table[1:otu_number,])
rich <-colSums(otu_table >0)
metadata=data.frame(sample_id=rownames(meta_table),
                    biome=biomes,
                    pH=meta_table$pH,
                    n=meta_table$N,
                    c=meta_table$C,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=abs(meta_table$Latitude),
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
cor(metadata$abundance,abs(metadata$Latitude))
cor(metadata$abundance,metadata$pH)
cor(metadata$abundance,metadata$C_N_ratio)
cor(metadata$abundance,metadata$Elevation)
cor(metadata$richness,abs(metadata$Latitude))
cor(metadata$richness,metadata$pH)
cor(metadata$richness,metadata$C_N_ratio)
cor(metadata$richness,metadata$Elevation)
cor(metadata$abundance,metadata$richness)

cor(metadata$pH,metadata$C_N_ratio)
cor(metadata$pH,abs(metadata$Latitude))
cor(metadata$pH,metadata$Elevation)
#guild table
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
#plant pathogen  and sap
plantpathsaprotroph <- table[grepl("plant pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(plantpathsaprotroph)
table$guild[grepl("plant pathogen",tolower(table$guild))] <- "Plant pathogens"
table$guild[grepl("fungal parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("lichen parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("animal parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("arbuscular mycorrhizal",tolower(table$guild))] <- "AM fungi"
#table$guild[grepl("plant pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild))] <- "Ecto-Saprotrophs"
ectosaprotroph <- table[grepl("ectomycorrhiza",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(ectosaprotroph)
table$guild[grepl("ectomycorrhizal",tolower(table$guild))] <- "Ectomycorrhizal fungi"
table$guild[grepl("saprotroph",tolower(table$guild))] <- "Saprotrophs"
table$guild[table$guild==""] <- "unknown"
table$guild[!(table$guild=="Plant pathogens" 
              | table$guild== "AM fungi" 
              | table$guild=="Ectomycorrhizal fungi"
              | table$guild=="Animal- and mycoparasites"
              | table$guild=="unknown"
              | table$guild=="Saprotrophs"
)] <- "others"

ecm <- subset(table,table$guild=="Ectomycorrhizal fungi")
ecm <-ecm[2:ncol(ecm)]
sap  <- subset(table,table$guild=="Saprotrophs")
sap<-sap[2:ncol(sap)]
pp <- subset(table,table$guild=="Plant pathogens")
pp<-pp[2:ncol(pp)]
am <- subset(table,table$guild=="AM fungi")
am<-am[2:ncol(am)]
pa <- subset(table,table$guild=="Animal- and mycoparasites")
pa<-pa[2:ncol(pa)]

#AT
m<-subset(metadata,metadata$biome=="AT")
samples <- rownames(m)
subotus <- table[,samples]
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(m$pH,m$C_N_ratio)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)

#richness and abundance of Ecm in AT
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in AT
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in AT
subotus <- pp[,AT_samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in AT
subotus <- sap[,AT_samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in AT
subotus <- pa[,AT_samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)

#BF
m<-subset(metadata,metadata$biome=="BF")
samples <- rownames(m)
cor(m$pH,m$C_N_ratio)
subotus <- table[,samples]
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(m$pH,m$C_N_ratio)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)

#richness and abundance of Ecm in BF
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in BF
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in BF
subotus <- pp[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in BF
subotus <- sap[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in BF
subotus <- pa[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)

#GS
m<-subset(metadata,metadata$biome=="GS")
samples <- rownames(m)
subotus <- table[,samples]
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(m$pH,m$C_N_ratio)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)

#richness and abundance of Ecm in gs
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in GS
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in GS
subotus <- pp[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in GS
subotus <- sap[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in GS
subotus <- pa[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)

#MED
m<-subset(metadata,metadata$biome=="MED")
samples <- rownames(m)
subotus <- table[,samples]
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(m$pH,m$C_N_ratio)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)

#richness and abundance of Ecm in MED
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in MED
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in MED
subotus <- pp[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in MED
subotus <- sap[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in MED
subotus <- pa[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)

#MTF
m<-subset(metadata,metadata$biome=="MTF")
samples <- rownames(m)
subotus <- table[,samples]
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(m$pH,m$C_N_ratio)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)

#richness and abundance of Ecm in MTF
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in MTF
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in MTF
subotus <- pp[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in MTF
subotus <- sap[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in MTF
subotus <- pa[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)

#SAV
m<-subset(metadata,metadata$biome=="SAV")
samples <- rownames(m)
subotus <- table[,samples]
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(m$pH,m$C_N_ratio)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)

#richness and abundance of Ecm in SAV
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in SAV
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in SAV
subotus <- pp[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in SAV
subotus <- sap[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in SAV
subotus <- pa[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)


#DTF
m<-subset(metadata,metadata$biome=="DTF")
samples <- rownames(m)
subotus <- table[,samples]
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(m$pH,m$C_N_ratio)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)

#richness and abundance of Ecm in DTF
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in DTF
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in DTF
subotus <- pp[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in DTF
subotus <- sap[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in DTF
subotus <- pa[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)

#STF
m<-subset(metadata,metadata$biome=="STF")
samples <- rownames(m)
subotus <- table[,samples]
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(m$pH,m$C_N_ratio)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)

#richness and abundance of Ecm in STF
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in STF
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in STF
subotus <- pp[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in STF
subotus <- sap[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in STF
subotus <- pa[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)


#TMF
m<-subset(metadata,metadata$biome=="TMF")
cor(m$pH,m$C_N_ratio)
samples <- rownames(m)
#richness and abundance of Ecm in TMF
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in TMF
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in TMF
subotus <- pp[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in TMF
subotus <- sap[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in TMF
subotus <- pa[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)



#TDF
m<-subset(metadata,metadata$biome=="TDF")
samples <- rownames(m)
subotus <- table[,samples]
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(m$pH,m$C_N_ratio)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)

#richness and abundance of Ecm in TDF
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in TDF
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in TDF
subotus <- pp[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in TDF
subotus <- sap[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in TDF
subotus <- pa[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)


#TCF
m<-subset(metadata,metadata$biome=="TCF")
samples <- rownames(m)
subotus <- table[,samples]
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(m$pH,m$C_N_ratio)
cor(rich,m$n)
cor(rich,m$c)
cor(abund,m$n)
cor(abund,m$c)

#richness and abundance of Ecm in TCF
subotus <- ecm[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of am in TCF
subotus <- am[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pp in TCF
subotus <- pp[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of sap in TCF
subotus <- sap[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)
#richness and abundance of pa in TCF
subotus <- pa[,samples]
dim(subotus)
subotus <-t(subotus)
abund <-rowSums(subotus)
rich <-rowSums(subotus >0)
cor(rich,abs(m$Latitude))
cor(rich,m$pH)
cor(rich,m$C_N_ratio)
cor(rich,m$Elevation)
cor(abund,abs(m$Latitude))
cor(abund,m$pH)
cor(abund,m$C_N_ratio)
cor(abund,m$Elevation)
cor(abund,rich)



#########################################
#Species indicator----
#install.packages("indicspecies")
library(indicspecies)

df <-data.frame(OTU_ID = row.names(guild_table), 
                        taxa = as.character(guild_table$genus),
                        otu_table,row.names = 1, 
                        check.names=FALSE,
                        stringsAsFactors = FALSE)
df <-df[order(df$taxa),]
#replace empty class with "unidentified"
df$taxa[df$taxa==""] <- "unidentified"
#replace unknown classes with "unidentified"
df$taxa[startsWith(df$taxa,"Branch0")] <- "unidentified"
df$taxa[grepl("Incertae_sedis",df$taxa)] <- "unidentified"
df$taxa[startsWith(df$taxa,"GS")] <- "unidentified"


#order again the table
df <-df[order(df$taxa),]
#change the number to numeric
df[,2:ncol(df)]= apply(df[,2:ncol(df)], 2, function(x) as.numeric(as.character(x)))
#merge all the order with the same name
df <- aggregate(df[2:ncol(df)], by=list(df$taxa), FUN = "sum")
taxa <- df$Group.1
df <- t(df)
colnames(df)=taxa
df=df[-1,]


class(df) <-"double" 
df<-subset(df,rowSums(df)!=0)
#Convert to relative frequencies
df<-df/rowSums(df)
#Use multipatt to find significant environmental variables
inv = multipatt(df, biomes, func = "r.g", control = how(nperm=9999))

summary(inv)



#########################################
#Species indicator plant pathogen ----
#install.packages("indicspecies")
library(indicspecies)

df <-data.frame(OTU_ID = row.names(guild_table), 
                Guild=guild_table$Guild,
                taxa = as.character(guild_table$genus),
                otu_table,
                row.names = 1, 
                check.names=FALSE,
                stringsAsFactors = FALSE)
df <- subset(df, grepl("plant pathogen",tolower(df$Guild)))
df <-df[2:ncol(df)]
df <-df[order(df$taxa),]
#replace empty class with "unidentified"
df$taxa[df$taxa==""] <- "unidentified"
#replace unknown classes with "unidentified"
df$taxa[startsWith(df$taxa,"Branch0")] <- "unidentified"
df$taxa[grepl("Incertae_sedis",df$taxa)] <- "unidentified"
df$taxa[startsWith(df$taxa,"GS")] <- "unidentified"


#order again the table
df <-df[order(df$taxa),]
#change the number to numeric
df[,2:ncol(df)]= apply(df[,2:ncol(df)], 2, function(x) as.numeric(as.character(x)))
#merge all the order with the same name
df <- aggregate(df[2:ncol(df)], by=list(df$taxa), FUN = "sum")
taxa <- df$Group.1
df <- t(df)
colnames(df)=taxa
df=df[-1,]


class(df) <-"double" 
df<-subset(df,rowSums(df)!=0)
#Convert to relative frequencies
df<-df/rowSums(df)
#Use multipatt to find significant environmental variables
inv = multipatt(df, biomes, func = "r.g", control = how(nperm=9999))

summary(inv)


#########################################
#CCA environmental information----
#CCA
abund <-colSums(otu_table[1:otu_number,])
rich <-colSums(otu_table >0)
metadata=data.frame(sample_id=rownames(meta_table),
                   pH=meta_table$pH,
                  
                   C_N_ratio=meta_table$C/meta_table$N,
                   Elevation=meta_table$Elevation,
                   Latitude=meta_table$Latitude,
                   richness=rich,
                   abundance=abund,
                   row.names = 1, 
                   check.names=FALSE,
                   stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(otu_table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 

#abund_table.mrpp <-with(data, mrpp(abund_table,data$pH))
#abund_table.mrpp
#write.table(abund_table.adonis,file = "results/adonis_env_biome.txt")



bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
#df_sites<-data.frame(scrs$sites,biomes,countries)
#colnames(df_sites)<-c("x","y","Biome","Country")
df_sites<-data.frame(scrs$sites,biomes)
colnames(df_sites)<-c("x","y","Biome")

#Draw sites
#locationnumber <- length(unique(countries))
#shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
#p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
#  scale_shape_manual(values=shape_values)
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome))
#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),
               aes(x, y, label = rownames(df_arrows)),
               color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_biome.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_cordinates.txt", sep = "\t", quote = FALSE, row.names = T)


#########################################
#CCA plant pathogen based on countries and biomes and environmental information----
table <-data.frame(OTU_ID = row.names(guild_table), 
                   otu_table,row.names = 1, 
                   Guild=guild_table$Guild,
                   check.names=FALSE,
                   stringsAsFactors = FALSE)
table <- subset(table, grepl("plant pathogen",tolower(table$Guild)))
table <-table[1:(ncol(table)-1)]
abund <-colSums(table)
rich <-colSums(table >0)
#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)


bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_biome.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_cordinates.txt", sep = "\t", quote = FALSE, row.names = T)

#########################################
#CCA ectomycorrhizal based on countries and biomes and environmental information----
table <-data.frame(OTU_ID = row.names(guild_table), 
                   otu_table,row.names = 1, 
                   Guild=guild_table$Guild,
                   check.names=FALSE,
                   stringsAsFactors = FALSE)
table <- subset(table, grepl("ectomycorrhizal",tolower(table$Guild)))
table <-table[1:(ncol(table)-1)]
abund <-colSums(table)
rich <-colSums(table >0)
#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_biome.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_cordinates.txt", sep = "\t", quote = FALSE, row.names = T)

#########################################
#CCA arbuscular mycorrhizal based on countries and biomes and environmental information----
table <-data.frame(OTU_ID = row.names(guild_table), 
                   otu_table,row.names = 1, 
                   Guild=guild_table$Guild,
                   check.names=FALSE,
                   stringsAsFactors = FALSE)
table <- subset(table, grepl("arbuscular mycorrhizal",tolower(table$Guild)))
table <-table[1:(ncol(table)-1)]
abund <-colSums(table)
rich <-colSums(table >0)
#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_biome.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_cordinates.txt", sep = "\t", quote = FALSE, row.names = T)

#########################################
#CCA saprotroph based on countries and biomes and environmental information----
table <-data.frame(OTU_ID = row.names(guild_table), 
                   otu_table,row.names = 1, 
                   Guild=guild_table$Guild,
                   check.names=FALSE,
                   stringsAsFactors = FALSE)
table <- subset(table, grepl("saprotroph",tolower(table$Guild)))
table <-table[1:(ncol(table)-1)]
abund <-colSums(table)
rich <-colSums(table >0)
#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)


bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_biome.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_cordinates.txt", sep = "\t", quote = FALSE, row.names = T)

#########################################
#CCA parasite based on countries and biomes and environmental information----
table <-data.frame(OTU_ID = row.names(guild_table), 
                   otu_table,row.names = 1, 
                   Guild=guild_table$Guild,
                   check.names=FALSE,
                   stringsAsFactors = FALSE)
table <- subset(table, 
                grepl("fungal parasite",tolower(table$Guild)) ||  
                grepl("animal parasite",tolower(table$Guild)) |
                grepl("lichen parasite",tolower(table$Guild)))
table <-table[1:(ncol(table)-1)]
abund <-colSums(table)
rich <-colSums(table >0)
#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)


bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_biome.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_cordinates.txt", sep = "\t", quote = FALSE, row.names = T)

########################################
#CCA AT----
table <- data.frame(Sample_ID=colnames(otu_table), 
                 biome=biomes,
                 t(otu_table),
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE)
table <- subset(table, table$biome=="AT")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
AT_otunames <-rownames(table)
#select only the OTUs of AT
table <-otu_table[AT_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_AT.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_AT.txt", sep = "\t", quote = FALSE, row.names = T)

########################################
#CCA GS ----
table <- data.frame(Sample_ID=colnames(otu_table), 
                    biome=biomes,
                    t(otu_table),
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
table <- subset(table, table$biome=="GS")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
GS_otunames <-rownames(table)
table <-otu_table[GS_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)


bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_AT.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_AT.txt", sep = "\t", quote = FALSE, row.names = T)

########################################
#CCA DTF ----
table <- data.frame(Sample_ID=colnames(otu_table), 
                    biome=biomes,
                    t(otu_table),
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
table <- subset(table, table$biome=="DTF")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
DTF_otunames <-rownames(table)
table <-otu_table[DTF_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)


bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_AT.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_AT.txt", sep = "\t", quote = FALSE, row.names = T)

########################################
#CCA MED ----
table <- data.frame(Sample_ID=colnames(otu_table), 
                    biome=biomes,
                    t(otu_table),
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
table <- subset(table, table$biome=="MED")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
MED_otunames <-rownames(table)
table <-otu_table[MED_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)


########################################
#CCA BF ----
table <- data.frame(Sample_ID=colnames(otu_table), 
                    biome=biomes,
                    t(otu_table),
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
table <- subset(table, table$biome=="BF")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
BF_otunames <-rownames(table)
table <-otu_table[BF_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)


bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_AT.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_AT.txt", sep = "\t", quote = FALSE, row.names = T)


########################################
#CCA TMF ----
table <- data.frame(Sample_ID=colnames(otu_table), 
                    biome=biomes,
                    t(otu_table),
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
table <- subset(table, table$biome=="TMF")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
TMF_otunames <-rownames(table)
table <-otu_table[TMF_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)


bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_AT.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_AT.txt", sep = "\t", quote = FALSE, row.names = T)

########################################
#CCA SAV ----
table <- data.frame(Sample_ID=colnames(otu_table), 
                    biome=biomes,
                    t(otu_table),
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
table <- subset(table, table$biome=="SAV")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
SAV_otunames <-rownames(table)
table <-otu_table[SAV_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)


bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_env_otus_AT.pdf",width = 12, height = 9)
print(p)
dev.off()
write.table(df,file = "results/otus_AT.txt", sep = "\t", quote = FALSE, row.names = T)

########################################
#CCA STF ----
table <- data.frame(Sample_ID=colnames(otu_table), 
                    biome=biomes,
                    t(otu_table),
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
table <- subset(table, table$biome=="STF")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
STF_otunames <-rownames(table)
table <-otu_table[STF_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)

########################################
#CCA TCF ----
table <- data.frame(Sample_ID=colnames(otu_table), 
                    biome=biomes,
                    t(otu_table),
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
table <- subset(table, table$biome=="TCF")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
TCF_otunames <-rownames(table)
table <-otu_table[TCF_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)

########################################
#CCA TDF ----
table <- data.frame(Sample_ID=colnames(otu_table), 
                    biome=biomes,
                    t(otu_table),
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
table <- subset(table, table$biome=="TDF")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
TDF_otunames <-rownames(table)
table <-otu_table[TDF_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)
########################################
#CCA MTF ----
table <- data.frame(Sample_ID=colnames(otu_table), 
                    biome=biomes,
                    t(otu_table),
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
table <- subset(table, table$biome=="MTF")
table <-table[2:ncol(table)]
table <-t(table)
table <- table[rowSums(table)>0,]
MTF_otunames <-rownames(table)
table <-otu_table[MTF_otunames,]
#CCA
abund <-colSums(table)
rich <-colSums(table >0)

metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    C_N_ratio=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    richness=rich,
                    abundance=abund,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)

abund_table <- t(table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#correlation
cor(data$abundance,abs(data$Latitude))
cor(data$abundance,data$pH)
cor(data$abundance,data$C_N_ratio)
cor(data$abundance,data$Elevation)
cor(data$richness,abs(data$Latitude))
cor(data$richness,data$pH)
cor(data$richness,data$C_N_ratio)
cor(data$richness,data$Elevation)
cor(data$abundance,data$richness)

#########################################
#CCA Class evn----
#load plant pathogen otus
class_table <-data.frame(OTU_ID = row.names(guild_table), 
                         class = as.character(guild_table$bioclass),
                         otu_table,
                         Guild=as.character(guild_table$Guild),
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
class_table$order[class_table$class==""] <- "unidentified"
#replace unknown classes with "unidentified"
class_table$class[startsWith(class_table$class,"Branch0")] <- "unidentified"
class_table$class[grepl("Incertae_sedis",class_table$class)] <- "unidentified"
class_table$class[startsWith(class_table$class,"GS")] <- "unidentified"
#order by name
class_table <- class_table[order(class_table$class),]
#merge all the genus with the same name
class_table <- aggregate(class_table[1:sample_number+1], by=list(class_table$class), FUN = "sum")

#save genus_table
#write.table(class_table,file = "results/class_level_only.txt", sep = "\t", quote = FALSE, row.names = F)

rownames(class_table) <-class_table$Group.1
abund_table<-class_table[1:sample_number+1]
abund_table <- t(abund_table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
abund_table<-abund_table/rowSums(abund_table)

#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    "C:N"=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#write.table(abund_table.adonis,file = "results/adonis_plant_pathogen_env_class.txt")

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first 
biomes<-meta_table$Biome
countries <-meta_table$Country
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")


#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
#p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country))
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)


#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_evn_class.pdf",width = 12, height = 9)
print(p)
dev.off()


#########################################
#Plant Path CCA Family ENV ----
#load plant pathogen otus
family_table <-data.frame(OTU_ID = row.names(guild_table), 
                          family = as.character(guild_table$family),
                          otu_table,
                          Guild=as.character(guild_table$Guild),
                          row.names = 1, 
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
family_table$family[family_table$family==""] <- "unidentified"
#replace unknown classes with "unidentified"
family_table$family[startsWith(family_table$family,"Branch0")] <- "unidentified"
family_table$family[grepl("Incertae_sedis",family_table$family)] <- "unidentified"
family_table$family[startsWith(family_table$family,"GS")] <- "unidentified"

family_table <- subset(family_table,grepl("plant pathogen",tolower(family_table$Guild)))
#order by name
family_table <- family_table[order(family_table$family),]
#merge all the genus with the same name
family_table <- aggregate(family_table[1:sample_number+1], by=list(family_table$family), FUN = "sum")

#save genus_table
write.table(family_table,file = "results/family_level_plant_pathogen_only.txt", sep = "\t", quote = FALSE, row.names = F)

rownames(family_table) <-family_table$Group.1
abund_table<-family_table[1:sample_number+1]
abund_table <- t(abund_table)
write.table(data.frame(names=colnames(abund_table), sums=colSums(abund_table)),
            file = "results/family_level_evn_plantpath_sums.txt", 
            sep = "\t", quote = FALSE, row.names = F)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
abund_table<-abund_table/rowSums(abund_table)

#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    "C:N"=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#write.table(abund_table.adonis,file = "results/adonis_plant_pathogen_env_family.txt")

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first 
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")


#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
#p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country))
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)


#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_plant_pathogen_evn_family.pdf",width = 12, height = 9)
print(p)
dev.off()

#########################################
#Plant path CCA Family ----

#load plant pathogen table
abund_table <-data.frame(OTU_ID = row.names(otu_table), 
                         guild=guild_table$Guild,
                         phylum = as.character(guild_table$phylum),
                         family = as.character(guild_table$family),
                         otu_table,
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
abund_table <- subset(abund_table, grepl("plant pathogen",tolower(abund_table$guild)))

#load Basidio and Asco. plant pathogens
df <-subset(abund_table,abund_table$phylum=="Basidiomycota" | abund_table$phylum=="Ascomycota")
#df <-subset(abund_table,abund_table$phylum=="Ascomycota")
#df <-subset(abund_table,abund_table$phylum=="Basidiomycota")

#merge family name
df <-df[3:ncol(df)]
df <-df[order(df$family),]
df<- subset(df,df$family!="")
df <- subset(df,df$family!="unidentified")
df <- subset(df,!grepl("Incertae_sedis",df$family))
#change the number to numeric
df[,2:ncol(df)]= apply(df[,2:ncol(df)], 2, function(x) as.numeric(as.character(x)))
df <- aggregate(df[2:ncol(df)], by=list(df$family), FUN = "sum")
rownames(df) <- df$Group.1
df <-df[2:ncol(df)]
write.table(data.frame(names=rownames(df), sums=rowSums(df)),
            file = "results/family_level_plantpath_sums.txt", 
            sep = "\t", quote = FALSE, row.names = F)

abund_table <-abund_table[4:ncol(abund_table)]
abund_table <-t(abund_table)
#Convert to relative frequencies
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
abund_table<-abund_table/rowSums(abund_table)
df <- df[rowSums(df) > 100,]
df<-t(df)
df<-df/rowSums(df)
data=as.data.frame(df)

#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#write.table(abund_table.adonis,file = "results/adonis_plant_pathogen_basidio_family.txt")

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(data))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes)
colnames(df_sites)<-c("x","y","Biome")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome))

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_plant_pathogen_familynames.pdf",width = 12, height = 9)
print(p)
dev.off()


#########################################
#Plant path CCA Genus ----

#load plant pathogen table
abund_table <-data.frame(OTU_ID = row.names(otu_table), 
                         guild=guild_table$Guild,
                         phylum = as.character(guild_table$phylum),
                         genus = as.character(guild_table$genus),
                         otu_table,
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
abund_table <- subset(abund_table, grepl("plant pathogen",tolower(abund_table$guild)))

#load Basidio and Asco. plant pathogens
df <-subset(abund_table,abund_table$phylum=="Basidiomycota" | abund_table$phylum=="Ascomycota")
#df <-subset(abund_table,abund_table$phylum=="Ascomycota")
#df <-subset(abund_table,abund_table$phylum=="Basidiomycota")

#merge family name
df <-df[3:ncol(df)]
df <-df[order(df$genus),]
df<- subset(df,df$genus!="")
df <- subset(df,df$genus!="unidentified")
df <- subset(df,!grepl("Incertae_sedis",df$genus))
#change the number to numeric
df[,2:ncol(df)]= apply(df[,2:ncol(df)], 2, function(x) as.numeric(as.character(x)))
df <- aggregate(df[2:ncol(df)], by=list(df$genus), FUN = "sum")
rownames(df) <- df$Group.1
df <-df[2:ncol(df)]
write.table(data.frame(names=rownames(df), sums=rowSums(df)),
            file = "results/genus_level_plantpath_sums.txt", 
            sep = "\t", quote = FALSE, row.names = F)

abund_table <-abund_table[4:ncol(abund_table)]
abund_table <-t(abund_table)
#Convert to relative frequencies
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
abund_table<-abund_table/rowSums(abund_table)
df <- df[rowSums(df) > 100,]
df<-t(df)
df<-df/rowSums(df)
data=as.data.frame(df)

#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(data))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes)
colnames(df_sites)<-c("x","y","Biome")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome))

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_plant_pathogen_genusnames.pdf",width = 12, height = 9)
print(p)
dev.off()



#########################################
#Functional diversity  richness and abundance----
#guild table
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
#plant pathogen  and sap
plantpathsaprotroph <- table[grepl("plant pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(plantpathsaprotroph)
table$guild[grepl("plant pathogen",tolower(table$guild))] <- "Plant pathogens"
table$guild[grepl("fungal parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("lichen parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("animal parasite",tolower(table$guild))] <- "Animal- and mycoparasites"
table$guild[grepl("arbuscular mycorrhizal",tolower(table$guild))] <- "AM fungi"
#table$guild[grepl("plant pathogen",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild))] <- "Ecto-Saprotrophs"
ectosaprotroph <- table[grepl("ectomycorrhiza",tolower(table$guild)) & grepl("saprotroph",tolower(table$guild)),]
dim(ectosaprotroph)
table$guild[grepl("ectomycorrhizal",tolower(table$guild))] <- "Ectomycorrhizal fungi"
table$guild[grepl("saprotroph",tolower(table$guild))] <- "Saprotrophs"
table$guild[table$guild==""] <- "unknown"
table$guild[!(table$guild=="Plant pathogens" 
              | table$guild== "AM fungi" 
              | table$guild=="Ectomycorrhizal fungi"
              | table$guild=="Animal- and mycoparasites"
              | table$guild=="unknown"
              | table$guild=="Saprotrophs"
)] <- "others"

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
richness_proportion <- df$Richness *100/sum(df$Richness)
abundance_proportion <- df$Abundance *100/sum(df$Abundance)
df <- data.frame(Sample_ID=guildnames, 
                 Richness=df$Richness,
                 Abundance=df$Abundance,
                 Richness_proportion =richness_proportion,
                 Abundance_proportion = abundance_proportion,
                 row.names = 1, 
                 check.names=FALSE,
                 stringsAsFactors = FALSE)

cor(df$Richness,df$Abundance)
#write to file
write.table(df,file = "results/richness_guild.txt", sep = "\t", quote = FALSE, row.names = T)

#########################################
#Taxa plot at Phylum level ----
phylum_table <-data.frame(OTU_ID = row.names(guild_table), 
                          phylum = as.character(guild_table$phylum),
                          otu_table,row.names = 1, 
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
phylum_table <-phylum_table[order(phylum_table$phylum),]
#replace unknown phyla with ""
phylum_table$phylum[phylum_table$phylum=="Fungi_phy_Incertae_sedis"] <- ""
phylum_table$phylum[phylum_table$phylum=="GS01"] <- ""
length(phylum_table$phylum[phylum_table$phylum!=""])#38140
#replace empty phylum with "unidentified"
phylum_table$phylum[phylum_table$phylum==""] <- "unidentified"
phylum_table <-phylum_table[order(phylum_table$phylum),]
#merge all the class with the same name
phylum_table <- aggregate(phylum_table[1:sample_number+1], by=list(phylum_table$phylum), FUN = "sum")
phylumnames <- phylum_table$Group.1
#save phylum_table
write.table(phylum_table,file = "results/phylum_level.txt", sep = "\t", quote = FALSE, row.names = T)
phylum_table <- t(phylum_table)
colnames(phylum_table)=phylumnames
phylum_table=phylum_table[-1,]
class(phylum_table) <- "double"

x <- phylum_table/rowSums(phylum_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

taxa_list <-colnames(x)
#remove "__unidentified__" and add it to others
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/phylum_level.pdf",height=6,width=21)
print(p)
dev.off()
#########################################
#Taxonomy at phylum level with biomes merged----
#Taxa plot at Phylum level
phylum_table <-data.frame(OTU_ID = row.names(guild_table), 
                          phylum = as.character(guild_table$phylum),
                          otu_table,row.names = 1, 
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
phylum_table <-phylum_table[order(phylum_table$phylum),]
#replace unknown phyla with ""
phylum_table$phylum[phylum_table$phylum=="Fungi_phy_Incertae_sedis"] <- ""
phylum_table$phylum[phylum_table$phylum=="GS01"] <- ""
length(phylum_table$phylum[phylum_table$phylum!=""])#38140
#replace empty phylum with "unidentified"
phylum_table$phylum[phylum_table$phylum==""] <- "unidentified"
phylum_table <-phylum_table[order(phylum_table$phylum),]
#merge all the class with the same name
phylum_table <- aggregate(phylum_table[1:sample_number+1], by=list(phylum_table$phylum), FUN = "sum")
phylumnames <- phylum_table$Group.1
phylum_table <-t(phylum_table)
phylum_table=phylum_table[-1,]
colnames(phylum_table) <- phylumnames
samplenames <-rownames(phylum_table)
#Get taxalist with respect to the order of all samples
class(phylum_table) <- "double"
x <- phylum_table/rowSums(phylum_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)

df <- data.frame(Sample_ID=samplenames, 
                 biome=biomes,
                 phylum_table,
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
#count all phyla 
s <- colSums(df[,1:ncol(df)])
df_all <- data.frame(t(s))
row.names(df_all) = c("ALL")
df <- rbind(df,df_all)


write.table(df,file = "results/phylum_level_biome.txt", sep = "\t", quote = FALSE, row.names = T)
                
x <- df/rowSums(df)

x<-x[,taxa_list]
#remove "__unidentified__" and add it to others
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/phylum_level_biome.pdf",height=6,width=21)
print(p)
dev.off()



#########################################
#Taxonomy at class level ----
class_table <-data.frame(OTU_ID = row.names(guild_table), 
                         class = as.character(guild_table$bioclass),
                         otu_table,
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
class_table <-class_table[order(class_table$class),]
#replace empty class with "unidentified"
class_table$class[class_table$class==""] <- "unidentified"
#replace unknown classes with "unidentified"
class_table$class[class_table$class=="Fungi_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Mucoromycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Pezizomycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Rozellomycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="GS01"] <- "unidentified"
class_table$class[class_table$class=="GS18"] <- "unidentified"
class_table$class[class_table$class=="GS25"] <- "unidentified"
class_table$class[class_table$class=="GS26"] <- "unidentified"
class_table$class[class_table$class=="GS27"] <- "unidentified"
class_table$class[class_table$class=="GS35"] <- "unidentified"
class_table$class[class_table$class=="GS37"] <- "unidentified"
class_table <-class_table[order(class_table$class),]
#merge all the class with the same name
class_table <- aggregate(class_table[1:sample_number+1], by=list(class_table$class), FUN = "sum")
classnames <- class_table$Group.1
class_table <- t(class_table)
colnames(class_table)=classnames
class_table=class_table[-1,]
#save class_table
write.table(class_table,file = "results/class_level.txt", sep = "\t", quote = FALSE, row.names = T)
class(class_table) <- "double"

x <- class_table/rowSums(class_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

taxa_list <-colnames(x)
#remove "__unidentified__" and add it to others
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00"#yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/class_level.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#Taxonomy at class level with biome merged----
class_table <-data.frame(OTU_ID = row.names(guild_table), 
                         class = as.character(guild_table$bioclass),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
class_table <-class_table[order(class_table$class),]
#replace empty class with "unidentified"
class_table$class[class_table$class==""] <- "unidentified"
#replace unknown classes with "unidentified"
class_table$class[class_table$class=="Fungi_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Mucoromycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Pezizomycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Rozellomycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="GS01"] <- "unidentified"
class_table$class[class_table$class=="GS18"] <- "unidentified"
class_table$class[class_table$class=="GS25"] <- "unidentified"
class_table$class[class_table$class=="GS26"] <- "unidentified"
class_table$class[class_table$class=="GS27"] <- "unidentified"
class_table$class[class_table$class=="GS35"] <- "unidentified"
class_table$class[class_table$class=="GS37"] <- "unidentified"
class_table <-class_table[order(class_table$class),]
#merge all the class with the same name
class_table <- aggregate(class_table[1:sample_number+1], by=list(class_table$class), FUN = "sum")
classnames <- class_table$Group.1
class_table <- t(class_table)
colnames(class_table)=classnames
class_table=class_table[-1,]
#save class_table
write.table(class_table,file = "results/class_level.txt", sep = "\t", quote = FALSE, row.names = T)
#Get taxalist with respect to the order of all samples
class(class_table) <- "double"
x <- class_table/rowSums(class_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)
df <- data.frame(Sample_ID=samplenames, 
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

write.table(df,file = "results/class_level_biome.txt", sep = "\t", quote = FALSE, row.names = T)

x <- df/rowSums(df)
x<-x[,taxa_list]
#remove "__unidentified__" and add it to others
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/class_level_biome.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#Taxonomy at order level ----
order_table <-data.frame(OTU_ID = row.names(guild_table), 
                         order = as.character(guild_table$order),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
order_table <-order_table[order(order_table$order),]
#replace empty class with "unidentified"
order_table$order[order_table$order==""] <- "unidentified"
#replace unknown classes with "unidentified"
order_table$order[startsWith(order_table$order,"Branch0")] <- "unidentified"
order_table$order[grepl("Incertae_sedis",order_table$order)] <- "unidentified"
order_table$order[startsWith(order_table$order,"GS")] <- "unidentified"

#change the name of insignificant groups to others
#s <- rowSums(order_table[,2:sample_number+1])
#order_table$order[order_table$order != "unidentified" & s<300] <- "others"

#order again the table
order_table <-order_table[order(order_table$order),]

#merge all the order with the same name
order_table <- aggregate(order_table[1:sample_number+1], by=list(order_table$order), FUN = "sum")
ordernames <- order_table$Group.1

#save order_table
write.table(order_table,file = "results/order_level.txt", sep = "\t", quote = FALSE, row.names = T)

order_table <- t(order_table)
colnames(order_table)=ordernames
order_table=order_table[-1,]
class(order_table) <- "double"

x <- order_table/rowSums(order_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

taxa_list <-colnames(x)

#select only 26 taxa to show
selectedtaxa <- taxa_list[1:26]

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

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]

#remove "__unidentified__" 
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-new_x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00"#yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/order_level.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#order_level_biome.pdf----
order_table <-data.frame(OTU_ID = row.names(guild_table), 
                         order = as.character(guild_table$order),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
order_table <-order_table[order(order_table$order),]
#replace empty class with "unidentified"
order_table$order[order_table$order==""] <- "unidentified"
#replace unknown classes with "unidentified"
order_table$order[startsWith(order_table$order,"Branch0")] <- "unidentified"
order_table$order[grepl("Incertae_sedis",order_table$order)] <- "unidentified"
order_table$order[startsWith(order_table$order,"GS")] <- "unidentified"

#change the name of insignificant groups to others
#s <- rowSums(order_table[,2:sample_number+1])
#order_table$order[order_table$order != "unidentified" & s<300] <- "others"

#order again the table
order_table <-order_table[order(order_table$order),]

#merge all the order with the same name
order_table <- aggregate(order_table[1:sample_number+1], by=list(order_table$order), FUN = "sum")
ordernames <- order_table$Group.1
order_table <- t(order_table)
colnames(order_table)=ordernames
order_table=order_table[-1,]
#save order_table
write.table(order_table,file = "results/order_level.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the order of all samples
class(order_table) <- "double"
x <- order_table/rowSums(order_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)

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

write.table(df,file = "results/order_level_biome.txt", sep = "\t", quote = FALSE, row.names = T)

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

x <- new_df/rowSums(new_df)

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]

#remove "__unidentified__" 
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-new_x[,taxa_list]


df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/order_level_biome.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#Taxonomy at family level ----
family_table <-data.frame(OTU_ID = row.names(guild_table), 
                         family = as.character(guild_table$family),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
family_table <-family_table[order(family_table$family),]
#replace empty class with "unidentified"
family_table$family[family_table$family==""] <- "unidentified"
#replace unknown classes with "unidentified"
family_table$family[startsWith(family_table$family,"Branch0")] <- "unidentified"
family_table$family[grepl("Incertae_sedis",family_table$family)] <- "unidentified"
family_table$family[startsWith(family_table$family,"GS")] <- "unidentified"

#change the name of insignificant groups to others
#s <- rowSums(family_table[,2:sample_number+1])
#family_table$family[family_table$family != "unidentified" & s<300] <- "others"

#order again the table
family_table <-family_table[order(family_table$family),]

#merge all the order with the same name
family_table <- aggregate(family_table[1:sample_number+1], by=list(family_table$family), FUN = "sum")
familynames <- family_table$Group.1

#save family_table
write.table(family_table,file = "results/family_level.txt", sep = "\t", quote = FALSE, row.names = T)

family_table <- t(family_table)
colnames(family_table)=familynames
family_table=family_table[-1,]

class(family_table) <- "double"
x <- family_table/rowSums(family_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)

#select only 26 taxa to show
selectedtaxa <- taxa_list[1:26]

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

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]

#remove "__unidentified__" 
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-new_x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00"#yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/family_level.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#family_level_biome.pdf ----
family_table <-data.frame(OTU_ID = row.names(guild_table), 
                          family = as.character(guild_table$family),
                          otu_table,row.names = 1, 
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
family_table <-family_table[order(family_table$family),]
#replace empty class with "unidentified"
family_table$family[family_table$family==""] <- "unidentified"
#replace unknown classes with "unidentified"
family_table$family[startsWith(family_table$family,"Branch0")] <- "unidentified"
family_table$family[grepl("Incertae_sedis",family_table$family)] <- "unidentified"
family_table$family[startsWith(family_table$family,"GS")] <- "unidentified"

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
write.table(family_table,file = "results/family_level.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the order of all samples
class(family_table) <- "double"
x <- family_table/rowSums(family_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)

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

write.table(df,file = "results/family_level_biome.txt", sep = "\t", quote = FALSE, row.names = T)

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

x <- new_df/rowSums(new_df)
#order x by the taxalist to be consistent with all samples

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]

#remove "__unidentified__" 
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-new_x[,taxa_list]


df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/family_level_biome.pdf",height=6,width=21)
print(p)
dev.off()

###########################
#Taxonomy at genus level ----
genus_table <-data.frame(OTU_ID = row.names(guild_table), 
                         genus = as.character(guild_table$genus),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
genus_table <-genus_table[order(genus_table$genus),]
#replace empty class with "unidentified"
genus_table$genus[genus_table$genus==""] <- "unidentified"
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
write.table(genus_table,file = "results/genus_level.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the colsums of all samples
class(genus_table) <- "double"
x <- genus_table/rowSums(genus_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)
#select only 26 taxa to show
selectedtaxa <- taxa_list[1:26]

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
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00"#yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/genus_level.pdf",height=6,width=21)
print(p)
dev.off()

###########################
#genus_level_biome.pdf----
genus_table <-data.frame(OTU_ID = row.names(guild_table), 
                         genus = as.character(guild_table$genus),
                         otu_table,row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
genus_table <-genus_table[order(genus_table$genus),]
#replace empty class with "unidentified"
genus_table$genus[genus_table$genus==""] <- "unidentified"
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
write.table(genus_table,file = "results/genus_level.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the colsums of all samples
class(genus_table) <- "double"
x <- genus_table/rowSums(genus_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)


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
write.table(df,file = "results/genus_level_biome.txt", sep = "\t", quote = FALSE, row.names = T)

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
                  Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),
                  Value=new_x[,i],
                  Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/genus_level_biome.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#Plant path evironmental variables ----

#load data
table <- data.frame(otu_table,guild_table)

#load plant pathogen 
plant_pathogen_table <- subset(table, grepl("plant pathogen",tolower(table$Guild)))

#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    "C:N"=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)
abund_table <- plant_pathogen_table[1:sample_number]
abund_table <- t(abund_table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#Convert to relative frequencies
abund_table<-abund_table/rowSums(abund_table)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#write.table(abund_table.adonis,file = "results/adonis_level.txt")

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")

#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country))+
  scale_shape_manual(values=shape_values)

#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/plant_pathogen_CCA.pdf",width = 12, height = 9)
print(p)
dev.off()

#########################################
#Plant path Hclust based on the labels (biome country) ----
#merge samples based on biomes and countries
table <-data.frame(OTU_ID = row.names(guild_table), 
                   otu_table,row.names = 1, 
                   Guild=as.character(guild_table$Guild),
                   check.names=FALSE,
                   stringsAsFactors = FALSE)
table <- subset(table,grepl("plant pathogen",tolower(table$Guild)))
table <-table[1:sample_number]
tr <-data.frame(sample_id=colnames(table),
                label=labels,
                t(table),
                row.names = 1, 
                check.names=FALSE,
                stringsAsFactors = FALSE)
tr <- tr[order(tr$label),]
#merge all rows with the same labels
otu_number_plant_path <-dim(table)[1]
tr <-aggregate(tr[1:otu_number_plant_path+1], by=list(tr$label), FUN = "sum")
#write to file
write.table(tr,file = "results/biome_country_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = F)
# Pairwise correlation between labels
m<-data.matrix(tr[1:otu_number_plant_path+1])
dd <- 1-cor(t(m), use = "pairwise.complete.obs", method = "pearson")
#hclust
hc <- hclust(as.dist(dd),method = "ward.D2")
hc$labels <-tr$Group.1
# Convert hclust into a dendrogram and plot
hcd <- as.dendrogram(hc)
#color the labels based on biomes
n <-length(unique(biomes))
color_vec <- RColorBrewer::brewer.pal(n, "Set3")
color_vec[2]="#000000" #change from yellow to black

uniquebiomes <-unique(biomes)
uniquelabels <-hc$labels
color_labels<-rep(color_vec[1],length(uniquelabels))
for (i in 1:length(uniquelabels)){
  for (j in 1:length(uniquebiomes)){
    if (startsWith(uniquelabels[i],uniquebiomes[j])){
      color_labels[i]<-color_vec[j]
      break
    }
  }
}

#plot
pdf("results/cluster_plant_pathogen.pdf", width = 12, height = 7)
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                cex = 0.7, col = "blue")

# reduced label size
par(cex=1, mar=c(20, 8, 4, 1))
#par(cex=1, mar=c(5, 8, 4, 1))

#install.packages('dendextend')
library(dendextend)
labels_colors(hcd) <- color_labels[order.dendrogram(hcd)]

plot(hcd, ylab = "Height", nodePar = nodePar,edgePar = list(col = 2:3, lwd = 2:1))

dev.off()


#########################################
#Plant pathogen Phylum level ----
phylum_table <-data.frame(OTU_ID = row.names(guild_table), 
                          phylum = as.character(guild_table$phylum),
                          otu_table,row.names = 1, 
                          Guild=guild_table$Guild,
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
phylum_table <- subset(phylum_table, grepl("plant pathogen",tolower(phylum_table$Guild)))
phylum_table <-phylum_table[order(phylum_table$phylum),]
#replace unknown phyla with ""
phylum_table$phylum[phylum_table$phylum=="Fungi_phy_Incertae_sedis"] <- ""
phylum_table$phylum[phylum_table$phylum=="GS01"] <- ""
length(phylum_table$phylum[phylum_table$phylum!=""])#38140
#replace empty phylum with "unidentified"
phylum_table$phylum[phylum_table$phylum==""] <- "unidentified"
phylum_table <-phylum_table[order(phylum_table$phylum),]
#merge all the class with the same name
phylum_table <- aggregate(phylum_table[1:sample_number+1], by=list(phylum_table$phylum), FUN = "sum")
phylumnames <- phylum_table$Group.1
#save phylum_table
write.table(phylum_table,file = "results/phylum_level_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)
phylum_table <- t(phylum_table)
colnames(phylum_table)=phylumnames
phylum_table=phylum_table[-1,]
class(phylum_table) <- "double"

x <- phylum_table/rowSums(phylum_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

taxa_list <-colnames(x)
new_x<-x[,taxa_list]
if ("unidentified" %in% colnames(x)) {
  #remove "__unidentified__" and add it to others
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
  new_x<-x[,taxa_list]
}
df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/phylum_level_plant_pathogen.pdf",height=6,width=21)
print(p)
dev.off()
#########################################
#phylum_level_biome_plant_pathogen.pdf ----
phylum_table <-data.frame(OTU_ID = row.names(guild_table), 
                          phylum = as.character(guild_table$phylum),
                          otu_table,row.names = 1, 
                          Guild=guild_table$Guild,
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
phylum_table <- subset(phylum_table, grepl("plant pathogen",tolower(phylum_table$Guild)))
phylum_table <-phylum_table[order(phylum_table$phylum),]
#replace unknown phyla with ""
phylum_table$phylum[phylum_table$phylum=="Fungi_phy_Incertae_sedis"] <- ""
phylum_table$phylum[phylum_table$phylum=="GS01"] <- ""
length(phylum_table$phylum[phylum_table$phylum!=""])#38140
#replace empty phylum with "unidentified"
phylum_table$phylum[phylum_table$phylum==""] <- "unidentified"
phylum_table <-phylum_table[order(phylum_table$phylum),]
#merge all the class with the same name
phylum_table <- aggregate(phylum_table[1:sample_number+1], by=list(phylum_table$phylum), FUN = "sum")
phylumnames <- phylum_table$Group.1
phylum_table <-t(phylum_table)
phylum_table=phylum_table[-1,]
colnames(phylum_table) <- phylumnames
samplenames <-rownames(phylum_table)
#Get taxalist with respect to the order of all samples
class(phylum_table) <- "double"
x <- phylum_table/rowSums(phylum_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)

df <- data.frame(Sample_ID=samplenames, 
                 biome=biomes,
                 phylum_table,
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
#count all phyla 
s <- colSums(df[,1:ncol(df)])
df_all <- data.frame(t(s))
row.names(df_all) = c("ALL")
df <- rbind(df,df_all)


write.table(df,file = "results/phylum_level_biome_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)

x <- df/rowSums(df)

taxa_list <-colnames(x)
new_x<-x[,taxa_list]
if ("unidentified" %in% colnames(x)) {
  #remove "__unidentified__" and add it to others
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
  new_x<-x[,taxa_list]
}

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/phylum_level_biome_plant_pathogen.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#Plant Path CCA Phylum ----
#load plant pathogen otus
phylum_table <-data.frame(OTU_ID = row.names(guild_table), 
                         phylum = as.character(guild_table$phylum),
                         otu_table,
                         Guild=as.character(guild_table$Guild),
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
phylum_table$order[phylum_table$phylum==""] <- "unidentified"
#replace unknown classes with "unidentified"
phylum_table$phylum[startsWith(phylum_table$phylum,"Branch0")] <- "unidentified"
phylum_table$phylum[grepl("Incertae_sedis",phylum_table$phylum)] <- "unidentified"
phylum_table$phylum[startsWith(phylum_table$phylum,"GS")] <- "unidentified"

phylum_table <- subset(phylum_table,grepl("plant pathogen",tolower(phylum_table$Guild)))
#order by name
phylum_table <- phylum_table[order(phylum_table$phylum),]
#merge all the genus with the same name
phylum_table <- aggregate(phylum_table[1:sample_number+1], by=list(phylum_table$phylum), FUN = "sum")

#save genus_table
write.table(phylum_table,file = "results/phylum_level_plant_pathogen_only.txt", sep = "\t", quote = FALSE, row.names = F)

rownames(phylum_table) <-phylum_table$Group.1
abund_table<-phylum_table[1:sample_number+1]
abund_table <- t(abund_table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
abund_table<-abund_table/rowSums(abund_table)

#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    "C:N"=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#write.table(abund_table.adonis,file = "results/adonis_plant_pathogen_env_class.txt")

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first 
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")


#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
#p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country))
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)


#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_plant_pathogen_evn_phylum.pdf",width = 12, height = 9)
print(p)
dev.off()


################################################
#Plant pathogen Class  class level ----
class_table <-data.frame(OTU_ID = row.names(guild_table), 
                         class = as.character(guild_table$bioclass),
                         otu_table,
                         Guild=guild_table$Guild,
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
class_table <- subset(class_table, grepl("plant pathogen",tolower(class_table$Guild)))
class_table <-class_table[order(class_table$class),]
#replace empty class with "unidentified"
class_table$class[class_table$class==""] <- "unidentified"
#replace unknown classes with "unidentified"
class_table$class[class_table$class=="Fungi_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Mucoromycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Pezizomycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Rozellomycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="GS01"] <- "unidentified"
class_table$class[class_table$class=="GS18"] <- "unidentified"
class_table$class[class_table$class=="GS25"] <- "unidentified"
class_table$class[class_table$class=="GS26"] <- "unidentified"
class_table$class[class_table$class=="GS27"] <- "unidentified"
class_table$class[class_table$class=="GS35"] <- "unidentified"
class_table$class[class_table$class=="GS37"] <- "unidentified"
class_table <-class_table[order(class_table$class),]
#merge all the class with the same name
class_table <- aggregate(class_table[1:sample_number+1], by=list(class_table$class), FUN = "sum")
classnames <- class_table$Group.1
class_table <- t(class_table)
colnames(class_table)=classnames
class_table=class_table[-1,]
#save class_table
write.table(class_table,file = "results/class_level_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)
class(class_table) <- "double"

x <- class_table/rowSums(class_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

taxa_list <-colnames(x)
#remove "__unidentified__" and add it to others
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00"#yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/class_level_plant_pathogen.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#class_level_biome_plant_pathogen.pdf----
class_table <-data.frame(OTU_ID = row.names(guild_table), 
                         class = as.character(guild_table$bioclass),
                         otu_table,
                         Guild=guild_table$Guild,
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
class_table <- subset(class_table, grepl("plant pathogen",tolower(class_table$Guild)))
class_table <-class_table[order(class_table$class),]
#replace empty class with "unidentified"
class_table$class[class_table$class==""] <- "unidentified"
#replace unknown classes with "unidentified"
class_table$class[class_table$class=="Fungi_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Mucoromycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Pezizomycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="Rozellomycotina_cls_Incertae_sedis"] <- "unidentified"
class_table$class[class_table$class=="GS01"] <- "unidentified"
class_table$class[class_table$class=="GS18"] <- "unidentified"
class_table$class[class_table$class=="GS25"] <- "unidentified"
class_table$class[class_table$class=="GS26"] <- "unidentified"
class_table$class[class_table$class=="GS27"] <- "unidentified"
class_table$class[class_table$class=="GS35"] <- "unidentified"
class_table$class[class_table$class=="GS37"] <- "unidentified"
class_table <-class_table[order(class_table$class),]
#merge all the class with the same name
class_table <- aggregate(class_table[1:sample_number+1], by=list(class_table$class), FUN = "sum")
classnames <- class_table$Group.1
class_table <- t(class_table)
colnames(class_table)=classnames
class_table=class_table[-1,]
#save class_table
write.table(class_table,file = "results/class_level_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)
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

write.table(df,file = "results/class_level_biome_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)

x <- df/rowSums(df)
x<-x[,taxa_list]
#remove "__unidentified__" and add it to others
taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
taxa_list <- c(taxa_list,"unidentified")
new_x<-x[,taxa_list]

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/class_level_biome_plant_pathogen.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#Plant Path CCA Class evn----
#load plant pathogen otus
class_table <-data.frame(OTU_ID = row.names(guild_table), 
                         class = as.character(guild_table$bioclass),
                         otu_table,
                         Guild=as.character(guild_table$Guild),
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
class_table$order[class_table$class==""] <- "unidentified"
#replace unknown classes with "unidentified"
class_table$class[startsWith(class_table$class,"Branch0")] <- "unidentified"
class_table$class[grepl("Incertae_sedis",class_table$class)] <- "unidentified"
class_table$class[startsWith(class_table$class,"GS")] <- "unidentified"

class_table <- subset(class_table,grepl("plant pathogen",tolower(class_table$Guild)))
#order by name
class_table <- class_table[order(class_table$class),]
#merge all the genus with the same name
class_table <- aggregate(class_table[1:sample_number+1], by=list(class_table$class), FUN = "sum")

#save genus_table
write.table(class_table,file = "results/class_level_plant_pathogen_only.txt", sep = "\t", quote = FALSE, row.names = F)

rownames(class_table) <-class_table$Group.1
abund_table<-class_table[1:sample_number+1]
abund_table <- t(abund_table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
abund_table<-abund_table/rowSums(abund_table)

#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    "C:N"=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#write.table(abund_table.adonis,file = "results/adonis_plant_pathogen_env_class.txt")

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first 
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")


#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
#p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country))
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)


#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_plant_pathogen_evn_class.pdf",width = 12, height = 9)
print(p)
dev.off()


#########################################
#Plant pathogens order level ----
order_table <-data.frame(OTU_ID = row.names(guild_table), 
                         order = as.character(guild_table$order),
                         otu_table,
                         Guild=guild_table$Guild,
                        row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
order_table <- subset(order_table, grepl("plant pathogen",tolower(order_table$Guild)))
order_table <-order_table[order(order_table$order),]
#replace empty class with "unidentified"
order_table$order[order_table$order==""] <- "unidentified"
#replace unknown classes with "unidentified"
order_table$order[startsWith(order_table$order,"Branch0")] <- "unidentified"
order_table$order[grepl("Incertae_sedis",order_table$order)] <- "unidentified"
order_table$order[startsWith(order_table$order,"GS")] <- "unidentified"

#order again the table
order_table <-order_table[order(order_table$order),]

#merge all the order with the same name
order_table <- aggregate(order_table[1:sample_number+1], by=list(order_table$order), FUN = "sum")
ordernames <- order_table$Group.1

#save order_table
write.table(order_table,file = "results/order_level_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)


order_table <- t(order_table)
colnames(order_table)=ordernames
order_table=order_table[-1,]
class(order_table) <- "double"

x <- order_table/rowSums(order_table)
x <- x[,order(colSums(x),decreasing=TRUE)]

taxa_list <-colnames(x)

#select only 26 taxa to show including "unidentified"
selectedtaxa <- taxa_list[1:26]
if (!("unidentified" %in% selectedtaxa) & ("unidentified" %in% taxa_list)){
  selectedtaxa <- c(selectedtaxa[1:25],"unidentified")
}

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

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]

if ("unidentified" %in% taxa_list){
  #remove "__unidentified__" 
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
  new_x<-new_x[,taxa_list]
}

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00"#yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/order_level_plant_pathogen.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#order_level_biome_plant_pathogen.pdf----
order_table <-data.frame(OTU_ID = row.names(guild_table), 
                         order = as.character(guild_table$order),
                         otu_table,
                         Guild=guild_table$Guild,
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
order_table <- subset(order_table, grepl("plant pathogen",tolower(order_table$Guild)))
order_table <-order_table[order(order_table$order),]
#replace empty class with "unidentified"
order_table$order[order_table$order==""] <- "unidentified"
#replace unknown classes with "unidentified"
order_table$order[startsWith(order_table$order,"Branch0")] <- "unidentified"
order_table$order[grepl("Incertae_sedis",order_table$order)] <- "unidentified"
order_table$order[startsWith(order_table$order,"GS")] <- "unidentified"

#change the name of insignificant groups to others
#s <- rowSums(order_table[,2:sample_number+1])
#order_table$order[order_table$order != "unidentified" & s<300] <- "others"

#order again the table
order_table <-order_table[order(order_table$order),]

#merge all the order with the same name
order_table <- aggregate(order_table[1:sample_number+1], by=list(order_table$order), FUN = "sum")
ordernames <- order_table$Group.1
order_table <- t(order_table)
colnames(order_table)=ordernames
order_table=order_table[-1,]
#save order_table
write.table(order_table,file = "results/order_level_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the order of all samples
class(order_table) <- "double"
x <- order_table/rowSums(order_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)

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

write.table(df,file = "results/order_level_biome_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)

#select only 26 taxa to show
#select only 26 taxa to show including "unidentified"
selectedtaxa <- taxa_list[1:26]
if (!("unidentified" %in% selectedtaxa) & ("unidentified" %in% taxa_list)){
  selectedtaxa <- c(selectedtaxa[1:25],"unidentified")
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

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]
if ("unidentified" %in% taxa_list){
  #remove "__unidentified__" 
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
  new_x<-new_x[,taxa_list]
}

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/order_level_biome_plant_pathogen.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#Plant Path CCA Order ----
#load plant pathogen otus
order_table <-data.frame(OTU_ID = row.names(guild_table), 
                          order = as.character(guild_table$order),
                          otu_table,
                          Guild=as.character(guild_table$Guild),
                          row.names = 1, 
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
order_table$order[order_table$order==""] <- "unidentified"
#replace unknown classes with "unidentified"
order_table$order[startsWith(order_table$order,"Branch0")] <- "unidentified"
order_table$order[grepl("Incertae_sedis",order_table$order)] <- "unidentified"
order_table$order[startsWith(order_table$order,"GS")] <- "unidentified"

order_table <- subset(order_table,grepl("plant pathogen",tolower(order_table$Guild)))
#order by name
order_table <- order_table[order(order_table$order),]
#merge all the genus with the same name
order_table <- aggregate(order_table[1:sample_number+1], by=list(order_table$order), FUN = "sum")

#save genus_table
write.table(order_table,file = "results/order_level_plant_pathogen_only.txt", sep = "\t", quote = FALSE, row.names = F)

rownames(order_table) <-order_table$Group.1
abund_table<-order_table[1:sample_number+1]
abund_table <- t(abund_table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
abund_table<-abund_table/rowSums(abund_table)

#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    "C:N"=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#write.table(abund_table.adonis,file = "results/adonis_plant_pathogen_env_order.txt")

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first 
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")


#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
#p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country))
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)


#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_plant_pathogen_evn_order.pdf",width = 12, height = 9)
print(p)
dev.off()


#########################################
#Plant pathogens family level ----
family_table <-data.frame(OTU_ID = row.names(guild_table), 
                          family = as.character(guild_table$family),
                          otu_table,
                          Guild=guild_table$Guild,
                          row.names = 1, 
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
family_table <- subset(family_table, grepl("plant pathogen",tolower(family_table$Guild)))
family_table <-family_table[order(family_table$family),]
#replace empty class with "unidentified"
family_table$family[family_table$family==""] <- "unidentified"
#replace unknown classes with "unidentified"
family_table$family[startsWith(family_table$family,"Branch0")] <- "unidentified"
family_table$family[grepl("Incertae_sedis",family_table$family)] <- "unidentified"
family_table$family[startsWith(family_table$family,"GS")] <- "unidentified"

#change the name of insignificant groups to others
#s <- rowSums(family_table[,2:sample_number+1])
#family_table$family[family_table$family != "unidentified" & s<300] <- "others"

#order again the table
family_table <-family_table[order(family_table$family),]

#merge all the order with the same name
family_table <- aggregate(family_table[1:sample_number+1], by=list(family_table$family), FUN = "sum")
familynames <- family_table$Group.1

#save family_table
write.table(family_table,file = "results/family_level_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)

family_table <- t(family_table)
colnames(family_table)=familynames
family_table=family_table[-1,]

class(family_table) <- "double"
x <- family_table/rowSums(family_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)

#select only 26 taxa to show
selectedtaxa <- taxa_list[1:26]

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

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]

if ("unidentified" %in% taxa_list){
  #remove "__unidentified__" 
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
  new_x<-new_x[,taxa_list]
}  

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00"#yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/family_level.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#family_level_biome_plant_pathogen.pdf ----
family_table <-data.frame(OTU_ID = row.names(guild_table), 
                          family = as.character(guild_table$family),
                          otu_table,
                          Guild=guild_table$Guild,
                          row.names = 1, 
                          check.names=FALSE,
                          stringsAsFactors = FALSE)
family_table <- subset(family_table, grepl("plant pathogen",tolower(family_table$Guild)))
family_table <-family_table[order(family_table$family),]
#replace empty class with "unidentified"
family_table$family[family_table$family==""] <- "unidentified"
#replace unknown classes with "unidentified"
family_table$family[startsWith(family_table$family,"Branch0")] <- "unidentified"
family_table$family[grepl("Incertae_sedis",family_table$family)] <- "unidentified"
family_table$family[startsWith(family_table$family,"GS")] <- "unidentified"

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
write.table(family_table,file = "results/family_level_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the order of all samples
class(family_table) <- "double"
x <- family_table/rowSums(family_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)

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

write.table(df,file = "results/family_level_biome_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)

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

x <- new_df/rowSums(new_df)
#order x by the taxalist to be consistent with all samples

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]
if ("unidentified" %in% taxa_list){
  #remove "__unidentified__" 
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
  new_x<-new_x[,taxa_list]
}  


df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/family_level_biome_plant_pathogen.pdf",height=6,width=21)
print(p)
dev.off()


#########################################
#Plant pathogens genus level ----
genus_table <-data.frame(OTU_ID = row.names(guild_table), 
                         genus = as.character(guild_table$genus),
                         otu_table,
                         Guild=guild_table$Guild,
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
genus_table <- subset(genus_table, grepl("plant pathogen",tolower(genus_table$Guild)))
genus_table <-genus_table[order(genus_table$genus),]
#replace empty class with "unidentified"
genus_table$genus[genus_table$genus==""] <- "unidentified"
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
write.table(genus_table,file = "results/genus_level_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the colsums of all samples
class(genus_table) <- "double"
x <- genus_table/rowSums(genus_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)
#select only 26 taxa to show
selectedtaxa <- taxa_list[1:26]

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

x <- new_df/rowSums(new_df)
x <- x[,order(colSums(x),decreasing=TRUE)]

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]

if ("unidentified" %in% taxa_list){
  #remove "__unidentified__" 
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
  new_x<-new_x[,taxa_list]
}  
df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=biomes)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00"#yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Samples")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/genus_level_plant_pathogen.pdf",height=6,width=21)
print(p)
dev.off()

###########################
#genus_level_biome_plant_pathogen.pdf----
genus_table <-data.frame(OTU_ID = row.names(guild_table), 
                         genus = as.character(guild_table$genus),
                         otu_table,
                         Guild=guild_table$Guild,
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
genus_table <- subset(genus_table, grepl("plant pathogen",tolower(genus_table$Guild)))
genus_table <-genus_table[order(genus_table$genus),]
#replace empty class with "unidentified"
genus_table$genus[genus_table$genus==""] <- "unidentified"
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
write.table(genus_table,file = "results/genus_level_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)

#Get taxalist with respect to the colsums of all samples
class(genus_table) <- "double"
x <- genus_table/rowSums(genus_table)
x <- x[,order(colSums(x),decreasing=TRUE)]
taxa_list <-colnames(x)


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
write.table(df,file = "results/genus_level_biome_plant_pathogen.txt", sep = "\t", quote = FALSE, row.names = T)

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

taxa_list <-c(selectedtaxa,"others")

x <- new_df/rowSums(new_df)
x <-x[,taxa_list]
#remove others

taxa_list <- c(selectedtaxa,"others")
new_x<-x[,taxa_list]

if ("unidentified" %in% taxa_list){
  #remove "__unidentified__" 
  taxa_list<-taxa_list[!grepl("unidentified",taxa_list)]
  taxa_list <- c(taxa_list,"unidentified")
  new_x<-new_x[,taxa_list]
}  

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),
                  Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),
                  Value=new_x[,i],
                  Type=rownames(x))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#22757e","#c59799","#a1d0ed","#401019","#FFFF00")
N<-length(taxa_list)
colours[N]="#FFFF00" #yellow
#colours[N]="#000000" #black
library(ggplot2)
p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:N])
p<-p+theme_bw()+ylab("Proportions")
p<-p+theme_bw()+xlab("Biomes")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
#p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p<-p+theme(axis.text.x=element_blank())
pdf("results/genus_level_biome_plant_pathogen.pdf",height=6,width=21)
print(p)
dev.off()

#########################################
#Plant Path CCA Genus ----
#load plant pathogen otus
genus_table <-data.frame(OTU_ID = row.names(guild_table), 
                         genus = as.character(guild_table$genus),
                         otu_table,
                         Guild=as.character(guild_table$Guild),
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
genus_table$genus[genus_table$genus==""] <- "unidentified"
#replace unknown classes with "unidentified"
genus_table$genus[startsWith(genus_table$genus,"Branch0")] <- "unidentified"
genus_table$genus[grepl("Incertae_sedis",genus_table$genus)] <- "unidentified"
genus_table$genus[startsWith(genus_table$genus,"GS")] <- "unidentified"

genus_table <- subset(genus_table,grepl("plant pathogen",tolower(genus_table$Guild)))
#order by name
genus_table <- genus_table[order(genus_table$genus),]
#merge all the genus with the same name
genus_table <- aggregate(genus_table[1:sample_number+1], by=list(genus_table$genus), FUN = "sum")

#save genus_table
write.table(genus_table,file = "results/genus_level_plant_pathogen_only.txt", sep = "\t", quote = FALSE, row.names = F)

rownames(genus_table) <-genus_table$Group.1
abund_table<-genus_table[1:sample_number+1]
abund_table <- t(abund_table)
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
abund_table<-abund_table/rowSums(abund_table)

#CCA
metadata=data.frame(sample_id=rownames(meta_table),
                    pH=meta_table$pH,
                    "C:N"=meta_table$C/meta_table$N,
                    Elevation=meta_table$Elevation,
                    Latitude=meta_table$Latitude,
                    row.names = 1, 
                    check.names=FALSE,
                    stringsAsFactors = FALSE)
data=as.data.frame(metadata)
#Use adonis to find significant environmental variables
abund_table.adonis <- adonis(abund_table ~ ., data)
abund_table.adonis 
#write.table(abund_table.adonis,file = "results/adonis_plant_pathogen_env_genus.txt")

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<0.05]
#remove NA entries
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables
#do all
#sol<-cca(abund_table ~ ., data)

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=as.data.frame(metadata))",sep="")))

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
#Check the attributes
# > attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first 
df_sites<-data.frame(scrs$sites,biomes,countries)
colnames(df_sites)<-c("x","y","Biome","Country")


#Draw sites
locationnumber <- length(unique(countries))
shape_values<-seq(1,locationnumber)
library(ggplot2)
p<-ggplot()
#p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country))
p<-p+geom_point(data=df_sites,aes(x,y,colour=Biome,shape=Country)) + 
  scale_shape_manual(values=shape_values)


#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)

# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")

p<-p+theme_bw()
pdf("results/CCA_plant_pathogen_evn_genus.pdf",width = 12, height = 9)
print(p)
dev.off()

#########################################
#Plant Path Heatmap Genus ----
#load plant pathogen otus
genus_table <-data.frame(OTU_ID = row.names(guild_table), 
                         genus = as.character(guild_table$genus),
                         otu_table,
                         Guild=as.character(guild_table$Guild),
                         row.names = 1, 
                         check.names=FALSE,
                         stringsAsFactors = FALSE)
genus_table$genus[genus_table$genus==""] <- "unidentified"
#replace unknown classes with "unidentified"
genus_table$genus[startsWith(genus_table$genus,"Branch0")] <- "unidentified"
genus_table$genus[grepl("Incertae_sedis",genus_table$genus)] <- "unidentified"
genus_table$genus[startsWith(genus_table$genus,"GS")] <- "unidentified"

genus_table <- subset(genus_table,grepl("plant pathogen",tolower(genus_table$Guild)))
#order by name
genus_table <- genus_table[order(genus_table$genus),]
#merge all the genus with the same name
genus_table <- aggregate(genus_table[1:sample_number+1], by=list(genus_table$genus), FUN = "sum")

#save genus_table
#write.table(genus_table,file = "results/genus_level_plant_pathogen_only.txt", sep = "\t", quote = FALSE, row.names = F)

rownames(genus_table) <-genus_table$Group.1
abund_table<-genus_table[1:sample_number+1]
#abund_table <- t(abund_table)
#abund_table<-subset(abund_table,rowSums(abund_table)!=0)
#abund_table<-abund_table/rowSums(abund_table)

n <-length(unique(biomes))
color_vec <- RColorBrewer::brewer.pal(n, "Set3")
color_vec[2]="#000000" #change from yellow to black

aa <- factor(biomes)
color_bar <- color_vec[aa]

pdf("results/genus_heatmap_withcolorbar_plant_pathogen.pdf", width = 12, height = 7)
b <- rowSums(abund_table[,1:sample_number])
ion$heatmap(abund_table[b > 100, 1:sample_number],
            col_color_bar = list("Biome" = color_bar),
            z_transform = "row",
            row_margin = 25,
            lwid = c(1, 6),
            #row_labels = colnames(otu_table)[1:sample_number], 
            #col_labels = labels, 
            #col_margin = 13, 
            col_distance = "spearman")
par(mar = c(0, 0, 0, 0), fig = c(0.75, 1, 0, 0.72), new = TRUE)
plot.new()
legend("topleft", levels(aa), col = color_vec,
       pch=15, pt.cex=1.5, bty = "n")
dev.off()

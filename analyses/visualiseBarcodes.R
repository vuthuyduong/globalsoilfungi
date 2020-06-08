#This file is to visualize the dataset in R. 
#set the working directory
setwd("C:/Users/Duong Vu/Documents/CBSPapers/SoilSamples-2020/data/CBS/CBSITSfull_classifiers")

#We first need to cluster the dataset. After that:
#Load the clustering result by clicking on the Save result of the main windows of fMLC in output file Yeast_CBS_GB_ITS_1seqfor1strain_Species.txt. 
output <- read.delim("CBSITS_withtaxa.classlevel", header = TRUE)

myoutput <-data.frame(index=output$Sequence.id, clusterindex = output$Reference.name)
y <-myoutput[order(myoutput$index),]
colorindexes <- y$clusterindex

#Load the coordinates obtained by LargeVis by clicking on Visualize button on the main windows of fMLC
x <- read.table("CBSITS_withtaxa.outLargeVis", skip = 1, sep= " ") 

mydata <- data.frame(index=x[,1],coor1=x[,2],coor2=x[,3],coor3 =x[,4]) #3D
x <-mydata[order(mydata$index),]

#count how many times a color appear for a reference cluster
freq<-sapply(1:length(colorindexes),function(x)sum(colorindexes[1:x]==colorindexes[x])) 
mydata <- data.frame(index=x[,1],coor1=x[,2],coor2=x[,3], coor3=x[,4], colorindexes=colorindexes, freq=freq)
ordereddata <-mydata[order(mydata$freq, decreasing=TRUE), ]

#get distinct colors for different groups
#install.packages("randomcoloR")
library(randomcoloR)
n <-length(unique(colorindexes))
palette <- distinctColorPalette(n)
#palette[1]=1 #black
#palette[2]=2 #red
#palette[3]=3 #green
#palette[4]=4 #blue
#palette[5]=5 #cyan
#palette[6]=6 #pink
#palette[7]=7 #yellow
#palette[8]=8 #grey
palette[1]=3 #green, Sordariomycetes
palette[2]=4 #blue, Saccharomycetes
palette[3]=5 #cyan, Eurotiomycetes
palette[4]=6 #pink, Dothideomycetes
palette[5]=2 #red, Agaricomycetes
palette[6]=8 #grey, Tremellomycetes
palette[7]=7  #yellow Leotiomycetes
palette[8]=1 #black Microbotryomycetes

#give color for each data point
newcolorindexes = c()
mylist = c()
index = 0
for (colorindex in ordereddata$colorindexes){
  #print(colorindex)
  if (length(mylist) >0) {
    index = match(colorindex,mylist)
  }
  if (is.na(index) || index == 0) {
    index = length(mylist)+1
    mylist=append(mylist,colorindex)
  } 
  newcolorindexes=append(newcolorindexes,palette[index])
}

#visualize the data
#install.packages("rgl")
library(rgl)
rgl.open()
rgl.points(ordereddata$coor1,ordereddata$coor2,ordereddata$coor3, color = newcolorindexes) 
rgl.bg(color = "white")
rgl.snapshot(filename = "CBSITS_withtaxa.png")
rgl.postscript("CBSITS_withtaxa.png",fmt="png")
rgl.postscript("CBSITS_withtaxa.pdf",fmt="pdf")
rgl.close() 

# Tammy Elliott and Ruan van Mazijk, 2020


# Load the required libraries.
library(sp)
library(rJava)
require(raster)
require(dismo)
library(rgdal)
library(maptools)
library(rgeos)
library(tidyverse)
library(caret)
require(caret)
library(faraway)
library(olsrr)
library(SSDM)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)
library(sf)
library(dplyr)
library(mapproj)
library(rgeos)


# Set working directory
setwd("/Users/tammy/PostDoc/Smuts/Schoeneae/Phylogeny/manuscript/Maps/wgsrpd-master/level3")
# double check to see that you are working in the correct working directory
getwd()

# import your species points
# this creates an object called georef
Sch.map<-read.csv("Schoenus-areas-world.csv",header=T)
head(Sch.map)

#Give species names as rownames
rownames(Sch.map)<-Sch.map$Species

# Calculate rowsums to get a vector with values per region
Sch.map.sums<-as.data.frame(colSums(Sch.map[,c(5:129)]))

Sch.map.sums.names<-data.frame(names=rownames(Sch.map.sums), Sch.map.sums)
colnames(Sch.map.sums.names)<-c("LEVEL3_COD","Count")

#Read map
map.shp <- readOGR(dsn="/Users/tammy/PostDoc/Smuts/Schoeneae/Phylogeny/manuscript/Maps/wgsrpd-master/level3",
	layer="level3", stringsAsFactors = FALSE)
map.shp.df<-fortify(map.shp, region="LEVEL3_COD")
	glimpse(map.shp)


#Create data frame for level 3 data
map.level.3<-as.data.frame(map.shp@data$"LEVEL3_COD")
colnames(map.level.3)<-"LEVEL3_COD"
map.level.3["Count"]<-"0"

#merge my data with level 3 names
sch.map.merge<-merge(map.level.3, Sch.map.sums.names, by="LEVEL3_COD",all=TRUE)
sch.map.merge[is.na(sch.map.merge)]=0

#Create data frame for level 3 data
#map.level.3<-as.data.frame(map.shp@data$"LEVEL3_COD")
#colnames(map.level.3)<-c("id", "Count")
#map.level.3["Count"]<-"0"

colnames(sch.map.merge)<-c("id", "Count.x", "Count.y")

#merge fortified shape file with my cound data for each species
level3.map.data.merge<-merge(map.shp.df, sch.map.merge, by="id", all=TRUE)





#Create map with ggplot

my.cols <- c("#FFFFFF","#FBF3F2","#FAEFEE","#F8EBEA","#F7E7E6","#F6E3E2", "#F5E0DE","#F5E0DE","#F0D0CD","#ECC4C1",
	"#E8B9B5","#E6B1AD","#D2736B","#C9584E","#B31205")


map <-ggplot() +
  geom_polygon(data = level3.map.data.merge, aes(fill = factor(Count.y), group=group, x=long, y=lat), color = "grey30", lwd=0.1) +
  theme_void()

#Get a vector of counts
count.vector<-sort(unique(sch.map.merge$Count.y))

#Create a map object
map2<-map + scale_fill_manual(name="Count", values=my.cols,
                       breaks=count.vector)
#Format the map
map3<-map2  +
	theme(legend.key.width=unit(0.04,"inch"),legend.key.height=unit(0.10,"inch")) +
	theme(legend.title=element_text(size=7),
	legend.text=element_text(size=5)) +
	theme(legend.position=c(0.035,0.5))

#Print for publication size
dev.new(width=4.5, height=3.2)
par(mar=c(1,0.5,0.1,0.3), mai=c(0.35,0.5,0.1,0.05))
print(map3)


#to share as PDF
map4<-map2  +
	theme(legend.key.width=unit(0.1,"inch"),legend.key.height=unit(0.2,"inch")) +
	theme(legend.title=element_text(size=10),
	legend.text=element_text(size=8)) +
	theme(legend.position=c(0.035,0.6))

dev.new()
pdf("Schoenus.world.pdf", width=12, height=9)
print(map4)
dev.off()

#Making dendrograms for associations between assays
# Paul Miller
# Jan 2016 

library(reshape2)
library(ggplot2)
library(RDRToolbox)
library(scales)
library(dendextend)

# Read in data
setwd('C:/Users/paulm/CRC Paul/PROJECTS/Plerixafor/figures/correlations/source data')
vitpre<-read.csv("plerixafor_prefreeze.csv",header=T,check.names=F)
#vitrop<-read.csv("plerixafor_pools_vitro.csv",header=T,check.names=T) #there is no in vivo counterpart
vitroi<-read.csv("plerixafor_vitro_individual3.csv",header=T,check.names=F)
vivoi<-read.csv("plerixafor_vivo_individual5.csv",header=T,check.names=F)
setwd('C:/Users/paulm/CRC Paul/PROJECTS/Plerixafor/figures/correlations/')


# Anonymize & integerize donor names
levels(vitpre$Donor) <- 1:10
levels(vitroi$Donor) <- 1:10
levels(vivoi$Donor) <- 1:10

# Rename Timepoints so they can be merged
levels(vitpre$Timepoint)[levels(vitpre$Timepoint)=="Day -1"] <- "-24 h"
levels(vitpre$Timepoint)[levels(vitpre$Timepoint)=="Day 0"] <- "+4 h"
levels(vitpre$Timepoint)[levels(vitpre$Timepoint)=="Day 1"] <- "+24 h"
levels(vitroi$Timepoint)[levels(vitroi$Timepoint)=="4hrs"] <- "+4 h"
levels(vitroi$Timepoint)[levels(vitroi$Timepoint)=="-24hrs"] <- "-24 h"
levels(vivoi$Timepoint)[levels(vivoi$Timepoint)=="4hrs"] <- "+4 h"
levels(vivoi$Timepoint)[levels(vivoi$Timepoint)=="-24hrs"] <- "-24 h"

# Change low values in PB leukocyte columns to detection threshold:
vivoi[,9:12][vivoi[,9:12] < 0.25] <- 0.25
# Change low values in platelet column to detection threshold:
vivoi[,13][vivoi[,13] < 25] <- 25
# Change low values in BM columns to detection threshold:
vivoi[,14:19][vivoi[,14:19] < 0.01] <- 0.01


# Calculate vivoi geometric means for triplicate mice
vivoi[,9:19]<-log10(vivoi[,9:19])
w3 <- NULL
w6 <- NULL
w12 <- NULL
for(i in 1:10){
	BLa <- vivoi[vivoi$Donor==i & vivoi$Timepoint=="BL" & vivoi$Week==3, 1:19]
	Dm1a <- vivoi[vivoi$Donor==i & vivoi$Timepoint=="-24 h" & vivoi$Week==3, 1:19] #
	D0a <- vivoi[vivoi$Donor==i & vivoi$Timepoint=="+4 h" & vivoi$Week==3, 1:19]
	BL <- BLa[1,]
	Dm1 <- Dm1a[1,] #
	D0 <- D0a[1,]
	BL[,9:19] <- colMeans(BLa[,9:19],na.rm=T)
	Dm1[,9:19] <- colMeans(Dm1a[,9:19],na.rm=T) #
	D0[,9:19] <- colMeans(D0a[,9:19],na.rm=T) 
	w3 <- rbind(w3,BL,Dm1,D0) #Dm1 #
}
for(i in 1:10){
  BLa <- vivoi[vivoi$Donor==i & vivoi$Timepoint=="BL" & vivoi$Week==6, 1:19]
  Dm1a <- vivoi[vivoi$Donor==i & vivoi$Timepoint=="-24 h" & vivoi$Week==6, 1:19] #
  D0a <- vivoi[vivoi$Donor==i & vivoi$Timepoint=="+4 h" & vivoi$Week==6, 1:19]
  BL <- BLa[1,]
  Dm1 <- Dm1a[1,] #
  D0 <- D0a[1,]
  BL[,9:19] <- colMeans(BLa[,9:19],na.rm=T)
  Dm1[,9:19] <- colMeans(Dm1a[,9:19],na.rm=T) #
  D0[,9:19] <- colMeans(D0a[,9:19],na.rm=T) 
  w6 <- rbind(w6,BL,Dm1,D0) #Dm1 #
}
for(i in 1:10){
  BLa <- vivoi[vivoi$Donor==i & vivoi$Timepoint=="BL" & vivoi$Week==12, 1:19]
  Dm1a <- vivoi[vivoi$Donor==i & vivoi$Timepoint=="-24 h" & vivoi$Week==12, 1:19]  #
  D0a <- vivoi[vivoi$Donor==i & vivoi$Timepoint=="+4 h" & vivoi$Week==12, 1:19]
  BL <- BLa[1,]
  Dm1 <- Dm1a[1,]  #
  D0 <- D0a[1,]
  BL[,9:19] <- colMeans(BLa[,9:19],na.rm=T)
  Dm1[,9:19] <- colMeans(Dm1a[,9:19],na.rm=T) #
  D0[,9:19] <- colMeans(D0a[,9:19],na.rm=T) 
  w12 <- rbind(w12,BL,Dm1,D0) #Dm1
}
vivoi <- rbind(w3,w6,w12)
vivoi[,9:19]<-10^(vivoi[,9:19])

# Drop mouse BM analysis because it is not to be included
vivoi <- vivoi[,1:13]

# Drop rows with only NA (these were attempted -24h calculations on donors that didn't have such samples)
vivoi <- vivoi[!is.na(vivoi$ID),]

# Put each analysis timepoint next to each other
w3 <- vivoi[vivoi$Week==3,]
w6 <- vivoi[vivoi$Week==6,]
w12 <- vivoi[vivoi$Week==12,]
names(w3)[9:13] <- paste("w3", colnames(w3[9:13]), sep=" ")
names(w6)[9:13] <- paste("w6", colnames(w6[9:13]), sep=" ")
names(w12)[9:13] <- paste("w12", colnames(w12[9:13]), sep=" ")
vivoi <- cbind(w3, w6[, 9:13], w12[, 9:13]) 


# Remove BM analysis from vitpre
vitpre <- vitpre[vitpre$Sample=="Blood",]
# Prefix label to prefreeze assays so they can be distinguished
names(vitpre)[6:ncol(vitpre)] <- paste("pf", colnames(vitpre[6:ncol(vitpre)]), sep=" ")
# Remove colNA's from vitpre
vitpre <- vitpre[, colSums(is.na(vitpre)) != nrow(vitpre)] #If column NA count = number of rows, must be entirely NA.

# Remove problematic columns: "Total CFC per 10e5", anything not per mL, "Experiment", "Mouse", "Week", "CFC.units", "LTCIC.units" 
vitpre <- vitpre[,c(1:8,10:13)]
vitroi <- vitroi[,2:18]
vivoi <- vivoi[,c(2:6,9:ncol(vivoi))]

# Merge!
dat <- merge(vitpre, vitroi, c("ID","Sample","Timepoint","Donor","Group"))
dat <- merge(dat, vivoi, c("ID","Sample","Timepoint","Donor","Group"))

# Rename group factors
levels(dat$Group)[levels(dat$Group)=="A"] <- "P"
levels(dat$Group)[levels(dat$Group)=="B"] <- "G+P"

# Log10 transform
dat[,6:ncol(dat)]<-log10(dat[,6:ncol(dat)])

# Label rows with sample numbers, and take out non-numeric data
attr <- dat[,1:5]
ndat <- dat[,c(3, 5, 6:ncol(dat))]

# Rename some assays so they can be merged
names(ndat)[names(ndat) == 'pf NCC x10e6/mL'] <- 'NCC x10e6/mL'
names(ndat)[names(ndat) == 'pf CD34/µL'] <- 'CD34/µL'

# Just take assays of interest. Exclude CD3, all pf.CFC, and non-total-CFC
ndat <- ndat[ ,c(1:3, 5, 13:24, 26:29, 31:34, 36)]


# Global dendrogram parameters:
hang_value <- 0.25  
labelSize <- 0.5

# Make nicer assay names, to be mapped to later on.
items <-c(
  expression(paste("TNC x 10"^"6","/mL")), 
  expression(paste("CD34"^"+"," x 10"^"3","/mL")), 
  expression(paste("CFC x 10"^"3","/mL")), 
  "Week 3 LTC-TNC", 
  "Week 6 LTC-TNC", 
  "Week 3 LTC-CD33/15", 
  "Week 6 LTC-CD33/15", 
  "Week 3 LTC-CD34", 
  "Week 6 LTC-CD34",
  "Week 3 LTC-CFC",
  "Week 6 LTC-CFC", 
  "Week 3 In Vivo CD45", 
  "Week 3 In Vivo CD33/15", 
  "Week 3 In Vivo CD19", 
  "Week 3 In Vivo Platelets", 
  "Week 6 In Vivo CD45", 
  "Week 6 In Vivo CD33/15", 
  "Week 6 In Vivo CD19", 
  "Week 6 In Vivo Platelets", 
  "Week 12 In Vivo CD45", 
  "Week 12 In Vivo CD33/15", 
  "Week 12 In Vivo CD19", 
  "Week 12 In Vivo Platelets")

# Make color code for assays
cols <- c(rep("#8B4513",2), rep("#EE7600",9), rep("#006400",12)) #brown orange green



###*** G-only ***###
# Select one sample-type:
Gdat <- ndat[ndat$Timepoint == "-24 h" & ndat$Group == "G+P", ]
Gdat <- Gdat[ , c(3:ncol(Gdat))] 

# Rank order each assay. 1 = lowest value, 20=highest. tied values get low
rdat <- apply(Gdat,2,rank, ties.method="min") 

# Transpose
tdat <- t(rdat)

## Hierarchical clustering
ddat <- hclust(dist(tdat))  # Uses "complete" method
ddat <- as.dendrogram(ddat)

# Color the branches based on the clusters:
ddat <- color_branches(ddat, k=8) #, groupLabels=assayTypes)

# Match the labels to items (assay names)
labels(ddat) <- 
  items[sort_levels_values(row(as.matrix(tdat[,1])))[order.dendrogram(ddat)]]

# Match the labels and colors
labels_colors(ddat) <- 
  cols[sort_levels_values(row(as.matrix(tdat[,1])))[order.dendrogram(ddat)]]

# Hang the dendrogram a bit:
ddat <- hang.dendrogram(ddat,hang_height = hang_value)

# Reduce the size of the labels:
ddat <- set(ddat, "labels_cex", labelSize)

# Store for later
Gdend <- ddat




###*** P-only ***###
# Select one sample-type:
Pdat <- ndat[ndat$Timepoint == "+4 h" & ndat$Group == "P", ]
Pdat <- Pdat[ , c(3:ncol(Pdat))] 

# Rank order each assay. 1 = lowest value, 20=highest. tied values get low
rdat <- apply(Pdat,2,rank, ties.method="min") 

# Transpose
tdat <- t(rdat)

## Hierarchical clustering
ddat <- hclust(dist(tdat))  # Uses "complete" method
ddat <- as.dendrogram(ddat)

# Color the branches based on the clusters:
ddat <- color_branches(ddat, k=8) #, groupLabels=assayTypes)

# Match the labels to items (assay names)
labels(ddat) <- 
  items[sort_levels_values(row(as.matrix(tdat[,1])))[order.dendrogram(ddat)]]

# Match the labels and colors
labels_colors(ddat) <- 
  cols[sort_levels_values(row(as.matrix(tdat[,1])))[order.dendrogram(ddat)]]

# Hang the dendrogram a bit:
ddat <- hang.dendrogram(ddat,hang_height = hang_value)

# Reduce the size of the labels:
ddat <- set(ddat, "labels_cex", labelSize)

# Store for later
Pdend <- ddat




###*** G+P ***###

# Select one sample-type:
GPdat <- ndat[ndat$Timepoint == "+4 h" & ndat$Group == "G+P", ]
GPdat <- GPdat[ , c(3:ncol(GPdat))] 

# Rank order each assay. 1 = lowest value, 20=highest. tied values get low
rdat <- apply(GPdat,2,rank, ties.method="min") 

# Transpose
tdat <- t(rdat)

## Hierarchical clustering
ddat <- hclust(dist(tdat))  # Uses "complete" method
ddat <- as.dendrogram(ddat)

# Color the branches based on the clusters:
ddat <- color_branches(ddat, k=7) #, groupLabels=assayTypes)
# Had to reduce k from 8 to 7 because of message:
# "In cutree.dendrogram(dend, k = k, h = h, order_clusters_as_data = FALSE,  :
# It is impossible to produce a one-to-one cut for the k/h you specidied. 0's have been introduced."

# Match the labels to items (assay names)
labels(ddat) <- 
  items[sort_levels_values(row(as.matrix(tdat[,1])))[order.dendrogram(ddat)]]

# Match the labels and colors
labels_colors(ddat) <- 
  cols[sort_levels_values(row(as.matrix(tdat[,1])))[order.dendrogram(ddat)]]
 
# Hang the dendrogram a bit
ddat <- hang.dendrogram(ddat,hang_height = hang_value)

# Reduce the size of the labels:
ddat <- set(ddat, "labels_cex", labelSize)

# Store for later
GPdend <- ddat






###*** All except BL ***###
# Select one sample-type:
Alldat <- ndat[ndat$Timepoint != "BL", ]
Alldat <- Alldat[ , c(3:ncol(Alldat))] 

# Rank order each assay. 1 = lowest value, 20=highest. tied values get low
rdat <- apply(Alldat,2,rank, ties.method="min") 

# Transpose
tdat <- t(rdat)

## Hierarchical clustering
ddat <- hclust(dist(tdat))  # Uses "complete" method
ddat <- as.dendrogram(ddat)

# Color the branches based on the clusters:
ddat <- color_branches(ddat, k=8) #, groupLabels=assayTypes)

# Match the labels to items (assay names)
labels(ddat) <- 
  items[sort_levels_values(row(as.matrix(tdat[,1])))[order.dendrogram(ddat)]]

# Match the labels and colors
labels_colors(ddat) <- 
  cols[sort_levels_values(row(as.matrix(tdat[,1])))[order.dendrogram(ddat)]]

# Hang the dendrogram a bit:
ddat <- hang.dendrogram(ddat,hang_height = 1)

# Reduce the size of the labels:
ddat <- set(ddat, "labels_cex", labelSize)

# Store for later
Alldend <- ddat



## Make plots

pdf("dendrograms.pdf", width=11.5/2.54, height=11.5/2.54)

# Global plotting parameters:
par(mfrow = c(2, 2), mar = c(1,0.25,0.75,5.25), mgp = c(3, -0.1, 0), cex.main = 0.8)
mainLine = -0.4

plot(Gdend, title=title("G only", line = mainLine),
    horiz =  TRUE,  nodePar = list(cex = .007), 
	cex.axis = 0.5, tck = -0.02)
	 
plot(Pdend, title=title("P only", line = mainLine),
    horiz =  TRUE,  nodePar = list(cex = .007), 
	cex.axis = 0.5, tck = -0.02)

plot(GPdend, title=title("G+P only", line = mainLine),
    horiz =  TRUE,  nodePar = list(cex = .007), 
	cex.axis = 0.5, tck = -0.02)
	 
plot(Alldend, title=title("G, P & G+P", line = mainLine),
    horiz =  TRUE,  nodePar = list(cex = .007), 
	cex.axis = 0.5, tck = -0.02)

dev.off()




# Orchard Common Garden - Microbial & Chemistry data cleaning, analysis, and visualization
# Set working directory and load necessary packages####
setwd("/Users/ellehorwath/Documents/Orchard_Common_Garden/commongarden")
# Install and load necessary packages ####
library(ANCOMBC)
library(effects)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(pairwiseAdonis)
library(lme4)

# ASV analysis: includes alpha diversity analysis with glm; beta diversity analysis with NMDS plots, PERMANOVA tests, and pairwise adonis; barchart plots at asv level for just the common garden; betadispersion

# Read in cleaned data ####
## METADATA
mdITS_OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #238 of 24 vars

mdITS_OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant", "Ecoregion")] <- lapply(mdITS_OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant", "Ecoregion")], as.factor)

mdITS_OCG <- droplevels(mdITS_OCG) #check to see if there are any empty levels in the factors

## ASV DATA
asvITS_OCG <- read.csv("data_csv/asvITS.OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T)
summary(rowSums(asvITS_OCG)) #0
summary(colSums(asvITS_OCG)) #1

### Remove duplicates from ASV
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v2', 'WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1')
asvITS_OCG <- asvITS_OCG[!rownames(asvITS_OCG) %in% rows_to_remove, ] 

## Remove negative control
asvITS_OCG <- asvITS_OCG[!(row.names(asvITS_OCG) == "NEG_8-28-21"),]

## Remove MTW.3.7.R_2012
asvITS_OCG <- asvITS_OCG[!(row.names(asvITS_OCG) == "MTW.3.7.R_2012"),] 
# asvITS_OCG[asvITS_OCG < 10] <- 0 # each observation needs at least 10 seqs.
# asvITS_OCG <- asvITS_OCG[rowSums(asvITS_OCG) > 0,] #the values that are greater than zero
summary(rowSums(asvITS_OCG)) #0
summary(colSums(asvITS_OCG)) #0

asvITS_OCG <- asvITS_OCG[,colSums(asvITS_OCG) > 4] # each asv needs at least 10 seqs.
asvITS_OCG <- asvITS_OCG[rowSums(asvITS_OCG) > 499,] # each samples needs at least 500 seqs.

#histogram of asv values
# hist(rowSums(asvITS_OCG), main="Histogram of ASV counts per sample", xlab="ASV counts", breaks=20)

summary(colSums(asvITS_OCG)) #6
summary(rowSums(asvITS_OCG)) #507

asvITS_OCG <- subset(asvITS_OCG, row.names(asvITS_OCG) %in% row.names(mdITS_OCG)) ##150 of 21 var
mdITS_OCG <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(asvITS_OCG)) ##166 of 21 var

row.names(asvITS_OCG) == row.names(mdITS_OCG) # sanity check:TRUE

mdITS_OCG <- droplevels(mdITS_OCG)
str(mdITS_OCG)

## TAXONOMY
tax.ITS <- read.csv("~/Documents/Orchard_Common_Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #taxonomy read in 5983 obs of 2 variables

# Alpha diversity asv level####
## Rarefying
asvITS_OCG.r <- rrarefy(asvITS_OCG, 500) ## rarefy: Warning message
asvITS_OCG.r <- asvITS_OCG.r[,colSums(asvITS_OCG.r) > 0] # each sample needs at least x seqs.
summary(colSums(asvITS_OCG.r)) #1
summary(rowSums(asvITS_OCG.r)) #500

asvITS_OCG.shannon <- diversity(asvITS_OCG.r)
asvITS_OCG.richness <- specnumber(asvITS_OCG.r)
asvITS_OCG.ef <- exp(asvITS_OCG.shannon)
asvITS_OCG.efr <- round(asvITS_OCG.ef)
mdITS_OCG <- cbind(mdITS_OCG, effective_species = asvITS_OCG.efr)
mdITS_OCG <- cbind(mdITS_OCG, richness = asvITS_OCG.richness)

glm.OCG <- glmer(richness ~ Subsp_ploidy + Year + (1|Plant), family = poisson, data=mdITS_OCG)
summary(glm.OCG) 

plot(allEffects(glm.OCG))

plot(mdITS_OCG$Year,asvITS_OCG.richness) 
plot(mdITS_OCG$Subsp_ploidy,asvITS_OCG.richness) 
plot(mdITS_OCG$Ecoregion,asvITS_OCG.richness) 

#Beta diversity - NMDS plots#### 
set.seed(23)
asvITS_OCG.nmds <- metaMDS(asvITS_OCG.r, trymax = 1000) ### Solution reached
# save(asvITS_OCG.nmds, file = "nmds/asvITS_OCG_nmds.rda") #save the nmds so you won't need to run it again
load("nmds/asvITS_OCG_nmds.rda") #load it to use in code anytime after the initial run

ordiplot(asvITS_OCG.nmds, type = "t",display = "sites",cex = .6)
rownames(asvITS_OCG.nmds$points) == rownames(mdITS_OCG)

#SUBSPECIES
plot(asvITS_OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="Sagebrush fungal community by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS_OCG$Subspecies],
     pch=c(19,17)[mdITS_OCG$Year])
legend("topleft",
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS_OCG.nmds,groups = mdITS_OCG$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(asvITS_OCG.nmds,groups = mdITS_OCG$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(asvITS_OCG.nmds,groups = mdITS_OCG$Subspecies, show.groups = "W", col = "goldenrod")

### PERMANOVA and adonis for subspecies and year ##
asvITS_OCG.subsp_yr <- adonis2(asvITS_OCG.r ~ mdITS_OCG$Subsp_ploidy * mdITS_OCG$Year, by = "margin")
asvITS_OCG.subsp_yr

# add in spatial distance data to PERMANOVA
asvITS_OCG.subsp_yr_spatial <- adonis2(asvITS_OCG.r ~ mdITS_OCG$Subsp_ploidy * mdITS_OCG$Year + mdITS_OCG$Latitude  + mdITS_OCG$Longitude)
asvITS_OCG.subsp_yr_spatial

# #PLOIDY
# plot(asvITS_OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
#      main="Sagebrush fungal community by ploidy", 
#      col= c("red","blue")[mdITS_OCG$Ploidy],
#      pch=c(19,17)[mdITS_OCG$Year])
# legend("topleft", 
#        legend=c("2n","4n"),
#        col= c("red","blue"),
#        pch=19,
#        cex=0.8,
#        bty = "n")
# ordispider(asvITS_OCG.nmds,groups = mdITS_OCG$Ploidy, show.groups = "2n", col = "red")
# ordispider(asvITS_OCG.nmds,groups = mdITS_OCG$Ploidy, show.groups = "4n", col = "blue")
# 
# ### PERMANOVA for ploidy##
# asvITS_OCG.ploidy <- adonis2(asvITS_OCG.r ~ mdITS_OCG$Ploidy,by="margin") 
# asvITS_OCG.ploidy #ploidy is not significant
# 
#SUBSPECIES PLOIDY
plot(asvITS_OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="Fungal community by subspecies, ploidy, & year",
     col= c("pink","brown","darkgreen",'tan','lightblue')[mdITS_OCG$Subsp_ploidy],
     pch=c(17,19)[mdITS_OCG$Year])
legend("topright",
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.6,
       bty = "n")
legend("topleft",
       legend=c("2012","2021"),
       col="black",
       pch=c(17,19),
       cex=0.6,
       bty = "n")
ordispider(asvITS_OCG.nmds,groups = mdITS_OCG$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(asvITS_OCG.nmds,groups = mdITS_OCG$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(asvITS_OCG.nmds,groups = mdITS_OCG$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(asvITS_OCG.nmds,groups = mdITS_OCG$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(asvITS_OCG.nmds,groups = mdITS_OCG$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")
#
# ### PERMANOVA and adonis for subspecies ploidy ##
# asvITS_OCG.subsp_ploi_yr_eco <- adonis2(asvITS_OCG.r ~ mdITS_OCG$Subsp_ploidy + mdITS_OCG$Year + mdITS_OCG$Ecoregion, by = "margin")
# asvITS_OCG.subsp_ploi_yr_eco

# #pairwiseadonis
asvITS_OCG.subsp.pw <- pairwise.adonis(asvITS_OCG.r, mdITS_OCG$Subsp_ploidy)
asvITS_OCG.subsp.pw #p-adjusted:T_4n vs V_2n = 0.04, T_4n vs V_4n = 0.03, T_4n vs W_4n = 0.05

#YEAR
plot(asvITS_OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by year", 
     col= c("maroon","cyan")[mdITS_OCG$Year],
     pch=19)
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=16,
       cex=0.8,
       bty = "n")
# ordiarrows(asvITS_OCG.nmds,mdITS_OCG$Plant)
ordiarrows(asvITS_OCG.nmds, mdITS_OCG$Plant)

## 2012 asv NMDS ####
#NMDS
asvITS.2012.r <- subset(asvITS_OCG.r, mdITS_OCG$Year=="2012") #105 of 2135 variables
summary(rowSums(asvITS.2012.r)) # 500
summary(colSums(asvITS.2012.r)) # 0

mdITS.2012 <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(asvITS.2012.r)) #101

set.seed(127)
asvITS.2012.nmds <- metaMDS(asvITS.2012.r, trymax = 500) ###solution reached! warning message
# save(asvITS.2012.nmds, file = "nmds/asvITS.2012.nmds.rda")
# load("nmds/asvITS.2012.nmds.rda")

ordiplot(asvITS.2012.nmds, type = "t",display = "sites",cex = .6)
rownames(asvITS.2012.nmds$points) == rownames(mdITS.2012)

#SUBSPECIES
plot(asvITS.2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2012 fungal community by subspecies", 
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.2012$Subspecies],
     pch=17)
legend("topright", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=17,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subspecies, show.groups = "W", col = "goldenrod")

#SUBSP PLOIDY
plot(asvITS.2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2012 fungal community by subspecies and ploidy", 
     col= c("pink","brown","darkgreen",'tan','lightblue')[mdITS.2012$Subsp_ploidy],
     pch=c(17))
legend("topright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=17,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

### PERMANOVA and adonis for subspecies ploidy ##
asvITS.2012.subsp_ploi <- adonis2(asvITS.2012.r ~ mdITS.2012$Subsp_ploidy, by = "margin") 
asvITS.2012.subsp_ploi 

# #pairwiseadonis
asvITS.2012.subsp.pw <- pairwise.adonis(asvITS.2012.r, mdITS.2012$Subsp_ploidy)
# asvITS.2012.subsp.pw #none sig

## 2021 asv NMDS ####
asvITS.2021.r <- subset(asvITS_OCG.r, mdITS_OCG$Year=="2021") #49 of 769 var 
summary(rowSums(asvITS.2021)) # 500
summary(colSums(asvITS.2021)) # 0

mdITS.2021 <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(asvITS.2021.r)) #54 obs of 22

set.seed(78)
asvITS.2021.nmds <- metaMDS(asvITS.2021.r, trymax=500) ###solution reached!
# save(asvITS.2021.nmds, file = "nmds/asvITS.2021.nmds.rda")
# load("nmds/asvITS.2021.nmds.rda")

ordiplot(asvITS.2021.nmds, type = "t",display = "sites",cex = .6)
rownames(asvITS.2021.nmds$points) == rownames(mdITS.2021)

#SUBSPECIES
plot(asvITS.2021.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2021 fungal community by subspecies", 
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.2021$Subspecies],
     pch=19)
legend("topleft", 
       legend=c("T","V","W"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subspecies, show.groups = "W", col = "goldenrod")

###PERMANOVA and adonis for 2021 subspecies##
asvITS.2021.ad <- adonis2(asvITS.2021.r ~ mdITS.2021$Subsp_ploidy, by = "margin") # Bray-Curtis is the default metric
asvITS.2021.ad #

#PLOIDY
plot(asvITS.2021.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="Sagebrush 2021 fungal community by ploidy",
     col= c("red","blue")[mdITS.2021$Ploidy],
     pch=19)
legend("topleft",
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Ploidy, show.groups = "2n", col = "red")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Ploidy, show.groups = "4n", col = "blue")

### PERMANOVA for 2021 ploidy##
asvITS.2021.ploidy <- adonis2(asvITS.2021.r ~ mdITS.2021$Ploidy,by="margin") # Bray-Curtis is the default metric
asvITS.2021.ploidy #ploidy is significant 0.048

mdITS.2021$Subsp_ploidy <- droplevels(mdITS.2021$Subsp_ploidy)

#SUBSPECIES PLOIDY
plot(asvITS.2021.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2021 fungal community by subspecies and ploidy", 
     col= c("pink","brown",'tan',"lightblue")[mdITS.2021$Subsp_ploidy],
     pch=19)
legend("topright", 
       legend=c("T_2n","T_4n","V_4n","W_4n"),
       col= c("pink","brown",'tan',"lightblue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")
#there is no V_2n in 2021

### PERMANOVA and adonis for subspecies ploidy ##
asvITS.2021.subsp_ploi <- adonis2(asvITS.2021.r ~ mdITS.2021$Subsp_ploidy, by = "margin") 
asvITS.2021.subsp_ploi #subspecies ploidy is not sig

#pairwiseadonis
asvITS_OCG.2021_subsp.pw <- pairwise.adonis(asvITS.2021.r, mdITS.2021$Subsp_ploidy)
asvITS_OCG.2021_subsp.pw # none sig

# ## T asv nmds ####
# #NMDS
# asvITS.tri.r <- subset(asvITS_OCG.r, mdITS_OCG$Subspecies=="T") 
# summary(rowSums(asvITS.tri.r)) # 500
# summary(colSums(asvITS.tri.r)) # 0
# 
# mdITS.tri <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(asvITS.tri.r)) #112 of 21
# str(mdITS.tri)

# 
# set.seed(2)
# asvITS.tri.nmds <- metaMDS(asvITS.tri.r, trymax=500) ###solution reached! warning message
# # save(asvITS.2012.nmds, file = "nmds/asvITS.2012.nmds.rda")
# load("nmds/asvITS.2012.nmds.rda")
# 
# ordiplot(asvITS.tri.nmds, type = "t",display = "sites",cex = .6)
# rownames(asvITS.tri.nmds$points) == rownames(mdITS.tri)
# 
# mdITS.tri <- droplevels(mdITS.tri)
# 
###PERMANOVA and adonis for tri year##
asvITS.tri.yr <- adonis2(asvITS.tri.r ~ mdITS.tri$Year + mdITS.tri$Subsp_ploidy)
asvITS.tri.yr # yr sig, subsp_ploidy marginally 

#Bar chart of ASV level ####
asvITSt <- t(asvITS_OCG.r)
asvITSto <- asvITSt[order(rowSums(asvITSt),decreasing = T),]
asvITStop <- asvITSto[1:59,]
unknown <- colSums(asvITSto[60:nrow(asvITSto),])
asvITStopo <- rbind(asvITStop,unknown)

customcol <- c("cadetblue4","royalblue3","darkblue","tomato1","dodgerblue2","cyan","darkred","purple","mediumblue","palegoldenrod","lightgoldenrod","indianred","yellow","purple4","darkgreen","lightsalmon","yellow3","purple2","lightblue","firebrick","navy","red4","red","darkmagenta","mediumvioletred","violetred2","skyblue","dodgerblue4","lightseagreen","gray")

barplot(asvITStopo,col=customcol,legend.text=F,axes=F,cex.names = .5,las=2, args.legend = list(x = "topleft", bty = "n", inset=c(-0.15, 0)))
barplot(asvITStopo,col=customcol,legend.text=T,axes=F,cex.names = .3,las=2, args.legend = list(x = "topleft", cex = 0.3, bty = "n", inset=c(-0.11, 0)))

### Barchart at other levels, just OCG
tax.ITS.p <- cbind(Feature.ID=rownames(tax.ITS),tax.ITS)
tax.ITS.p <- parse_taxonomy(tax.ITS.p)
tax.ITS.OCG.p <- subset(tax.ITS.p, rownames(tax.ITS.p) %in% rownames(t(asvITS_OCG.r)))

##### CLASS LEVEL OCG
asvITS_OCGtc <- data.frame(Class=tax.ITS.OCG.p$Class,t(asvITS_OCG.r))
asvITS_OCGtc$Class[is.na(asvITS_OCGtc$Class)] <- "Unknown"
asvITS_OCGtca <- aggregate(. ~ asvITS_OCGtc$Class, asvITS_OCGtc[,2:ncol(asvITS_OCGtc)], sum) 
row.names(asvITS_OCGtca) <- asvITS_OCGtca[,1]
asvITS_OCGtca <- asvITS_OCGtca[,2:ncol(asvITS_OCGtca)]

asvITS_OCGtcao <- as.matrix(asvITS_OCGtca[order(rowSums(asvITS_OCGtca),decreasing = T),])

customcol28 <- c("cadetblue4","royalblue3","darkblue","tomato1","dodgerblue2",
                 "cyan","darkred","purple","mediumblue","palegoldenrod",
                 "lightgoldenrod","indianred","yellow","purple4","darkgreen",
                 "lightsalmon","yellow3","purple2","lightblue","firebrick",
                 "navy","red4","red","darkmagenta","mediumvioletred",
                 "violetred2","skyblue","dodgerblue4")

barplot(asvITS_OCGtcao,col=customcol28,legend.text=F,axes=F,cex.names = .3,las=2, args.legend = list(x = "topleft", bty = "n", inset=c(-0.15, 0)))
barplot(asvITS_OCGtcao,col=customcol28,legend.text=T,axes=F,cex.names = .3,las=2, args.legend = list(x = "topleft", bty = "n", inset=c(-0.11, -0.1), cex=0.4))

# Betadispersion ####
#YEAR
modyr <- betadisper(vegdist(asvITS_OCG.r), factor(mdITS_OCG$Year, labels = c("2012","2021")))
anova(modyr) #not sig between years p = 0.1552

#SUBSPECIES
modsp <- betadisper(vegdist(asvITS_OCG.r), factor(mdITS_OCG$Subspecies, labels = c("T","V","W")))
anova(modsp) #no sig p = 0.097

#SUBSPECIES PLOIDY
modsubplo <- betadisper(vegdist(asvITS_OCG.r), factor(mdITS_OCG$Subsp_ploidy, labels = c("T_2n","T_4n","V_2n","V_4n","W_4n")))
anova(modsubplo) #sig = 0.004591
boxplot(modsubplo)

ordiplot(modsubplo, type = "t",display = "sites",cex = .6)
plot(modsubplo,ellipse = TRUE,hull = FALSE ,conf = 0.80)#90% data ellipse

permutest(modsubplo, permutations = 999, pairwise = TRUE)

modsubplo.HSD<-TukeyHSD(modsubplo) #sig difference detected
plot(modsubplo.HSD)

#ANCOM: analysis of composition of microbiomes. Differential abundance analysis for common garden####
asvITS_OCG.t <- t(asvITS_OCG)
asvITS_OCG.t <- asvITS_OCG.t[,order(colnames(asvITS_OCG.t))] # order samples alphabetically

mdITS_OCG <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(asvITS_OCG)) 
row.names(asvITS_OCG) == row.names(mdITS_OCG)# sanity check :TRUE

taxITS.OCG <- subset(tax.ITS, row.names(tax.ITS) %in% rownames(asvITS_OCG.t))
row.names(taxITS.OCG) == rownames(asvITS_OCG.t)  # sanity check: TRUE

# #Create ANCOM Function ####
# 
# ancom.W = function(otu_data,var_data,
#                    adjusted,repeated,
#                    main.var,adj.formula,
#                    repeat.var,long,rand.formula,
#                    multcorr,sig){
#   
#   n_otu=dim(otu_data)[2]-1
#   
#   otu_ids=colnames(otu_data)[-1]
#   
#   if(repeated==F){
#     data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL,check.names=FALSE)
#     #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
#   }else if(repeated==T){
#     data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
#     # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
#   }
#   
#   base.formula = paste0("lr ~ ",main.var)
#   if(repeated==T){
#     repeat.formula = paste0(base.formula," | ", repeat.var)
#   }
#   if(adjusted==T & repeated==F ){
#     adjusted.formula = paste0(base.formula," + ", adj.formula)
#   }
#   if(adjusted==T & repeated==T ){
#     adjusted.formula = paste0(base.formula," + ", adj.formula," | ", repeat.var)
#   }
#   
#   if( adjusted == F & repeated == F ){
#     fformula  <- formula(base.formula)
#   } else if( adjusted == F & repeated == T & long == T ){
#     fformula  <- formula(repeat.formula)   
#   }else if( adjusted == F & repeated == T & long == F ){
#     fformula  <- formula(repeat.formula)   
#   }else if( adjusted == T & repeated == F  ){
#     fformula  <- formula(adjusted.formula)   
#   }else if( adjusted == T & repeated == T  ){
#     fformula  <- formula(adjusted.formula)   
#   }else{
#     stop("Problem with data. Dataset should contain OTU abundances, groups, 
#          and optionally an ID for repeated measures.")
#   }
#   
#   
#   
#   if( repeated==FALSE & adjusted == FALSE){
#     if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
#       tfun <- exactRankTests::wilcox.exact
#     } else{
#       tfun <- stats::kruskal.test
#     }
#   }else if( repeated==FALSE & adjusted == TRUE){
#     tfun <- stats::aov
#   }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
#     tfun <- stats::friedman.test
#   }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
#     tfun <- nlme::lme
#   }else if( repeated== TRUE & adjusted == TRUE){
#     tfun <- nlme::lme
#   }
#   
#   logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
#   for(ii in 1:(n_otu-1)){
#     for(jj in (ii+1):n_otu){
#       data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
#       lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
#       
#       lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
#       
#       if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
#         logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
#       }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
#         logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
#       }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
#         model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
#         picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
#         logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
#       }else if(repeated==TRUE&long==TRUE){ ## GEE
#         model=tfun(fixed=fformula,data = lr_dat,
#                    random = formula(rand.formula),
#                    correlation=corAR1(),
#                    na.action=na.omit)   
#         picker=which(gsub(" ","",row.names(anova(model)))==main.var)
#         logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
#       }
#       
#     }
#   } 
#   
#   ind <- lower.tri(logratio.mat)
#   logratio.mat[ind] <- t(logratio.mat)[ind]
#   
#   
#   logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
#   
#   mc.pval <- t(apply(logratio.mat,1,function(x){
#     s <- p.adjust(x, method = "BH")
#     return(s)
#   }))
#   
#   a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
#   
#   b <- matrix(0,ncol=n_otu,nrow=n_otu)
#   b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
#   diag(b)  <- NA
#   ind.1    <- lower.tri(b)
#   b[ind.1] <- t(b)[ind.1]
#   
#   ######################## Function #################
#   ### Code to extract surrogate p-value
#   surr.pval <- apply(mc.pval,1,function(x){
#     s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
#     # s0=max(x[which(as.numeric(as.character(x))<alpha)])
#     return(s0)
#   })
#   ######################## Function #################
#   ### Conservative
#   if(multcorr==1){
#     W <- apply(b,1,function(x){
#       subp <- length(which(x<sig))
#     })
#     ### Moderate
#   } else if(multcorr==2){
#     W <- apply(mc.pval,1,function(x){
#       subp <- length(which(x<sig))
#     })
#     ### No correction
#   } else if(multcorr==3){
#     W <- apply(logratio.mat,1,function(x){
#       subp <- length(which(x<sig))
#     })
#   }
#   
#   return(W)
# }
# 
# 
# 
# ANCOM.main = function(OTUdat,Vardat,
#                       adjusted,repeated,
#                       main.var,adj.formula,
#                       repeat.var,longitudinal,
#                       random.formula,
#                       multcorr,sig,
#                       prev.cut){
#   
#   p.zeroes=apply(OTUdat[,-1],2,function(x){
#     s=length(which(x==0))/length(x)
#   })
#   
#   zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
#   colnames(zeroes.dist)=c("Taxon","Proportion_zero")
#   
#   zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
#     geom_histogram(binwidth=0.1,colour="black",fill="white") + 
#     xlab("Proportion of zeroes") + ylab("Number of taxa") +
#     theme_bw()
#   
#   #print(zero.plot)
#   
#   OTUdat.thinned=OTUdat
#   OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
#   
#   otu.names=colnames(OTUdat.thinned)[-1]
#   
#   W.detected   <- ancom.W(OTUdat.thinned,Vardat,
#                           adjusted,repeated,
#                           main.var,adj.formula,
#                           repeat.var,longitudinal,random.formula,
#                           multcorr,sig)
#   
#   W_stat       <- W.detected
#   
#   
#   ### Bubble plot
#   
#   W_frame = data.frame(otu.names,W_stat,row.names=NULL)
#   W_frame = W_frame[order(-W_frame$W_stat),]
#   
#   W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
#   W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
#   W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
#   W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
#   
#   W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
#   W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
#   W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
#   W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
#   
#   final_results=list(W_frame,zero.plot)
#   names(final_results)=c("W.taxa","PLot.zeroes")
#   return(final_results)
# }
# 
# 
# #####End of functions
# 
# #Create "Sample.ID" column for all data tables
# #ANCOM requires that data be formatted so that first *column* is named "Sample.ID"
# #Sample IDs as row names does not count!
# taxITS.OCG_t <- t(taxITS.OCG) #transpose ASV and taxonomy data frames so sample IDs become row names
# asvITS_OCG_sbst <- data.frame("Sample.ID" = row.names(asvITS_OCG), asvITS_OCG, check.names = F)#create new column in each of ASV, metadata, and taxonomy data frames labelled "Sample.ID" containing sample IDs
# mdITS_OCG_sbst <- data.frame("Sample.ID" = row.names(mdITS_OCG), mdITS_OCG)
# taxITS.OCG_sbst <- data.frame("Sample.ID" = row.names(taxITS.OCG_t), taxITS.OCG_t, check.names = F)
# 
# row.names(asvITS_OCG_sbst) == row.names(mdITS_OCG_sbst) #TRUE
# 
# #SUBSPECIES 
# ANCOM_subspecies <- ANCOM.main(asvITS_OCG_sbst,mdITS_OCG_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)
# #Create objects of significant ASVs
# sigASVs_subspecies <- subset(ANCOM_subspecies$W.taxa, ANCOM_subspecies$W.taxa$W_stat > 0)[,1]
# sigASVs_subspecies <- as.data.frame(sigASVs_subspecies) #convert significant ASV return to data frame
# row.names(sigASVs_subspecies) <- sigASVs_subspecies[1:8,1] #change row names to ASV sequence
# 
# #By default ANCOM sorts ASVs in order of decreasing abundance. This section instead aims to sort by decreasing abundance, which can be cross-referenced with the .csv generated after the break
# #If you want to look at the default ANCOM output, the full code will be included at the end of this section
# 
# sigASVs_subspecies[,1] <- c(1:8) #append numerical ranking of significance in additional column
# sigASVs_subspecies_t <- t(sigASVs_subspecies) #some transposition and fiddling to get ASV sequences as usable column names to mimic format of original data frames
# sigASVs_subspecies_t <- as.data.frame(sigASVs_subspecies_t)
# colnames(sigASVs_subspecies_t) <- as.character(colnames(sigASVs_subspecies_t))                   
# print(colnames(sigASVs_subspecies_t))                                                            
# sigASVs_subspecies_t <- sigASVs_subspecies_t[,order(colnames(sigASVs_subspecies_t))]#sort alphabetically by ASV sequence, this is how we will align with original ASV data frame subsets
# rownames(sigASVs_subspecies_t) <- c("sig_rank")#rename row to something more informative
# 
# write.csv(ANCOM_subspecies$W.taxa, file = "data_csv/ANCOM/ANCOM_subspecies.csv") #write ANCOM output to .csv to observe significance cutoffs
# 
# #Subset original data by significant ASVs
# asvITS_OCG_sigsbst_subspecies <-  t(subset(t(asvITS_OCG_sbst), colnames(asvITS_OCG_sbst) %in% row.names(sigASVs_subspecies))) #subset of ASV counts to only significant taxa
# taxITS.OCG_sigsbst_subspecies <-  t(subset(t(taxITS.OCG_sbst), colnames(taxITS.OCG_sbst) %in% row.names(sigASVs_subspecies))) #subset of taxonomy to only significant taxa
# 
# #This next section is a continuation of the above work to sort ASVs by decreasing significance
# asvITS_OCG_sigsbst_subspecies <- asvITS_OCG_sigsbst_subspecies[,order(colnames(asvITS_OCG_sigsbst_subspecies))] #sort alphabetically for alignment with significance ranking
# colnames(sigASVs_subspecies_t) == colnames(asvITS_OCG_sigsbst_subspecies) #sanity check:TRUE
# 
# asvITS_OCG_sigsbst_subspecies <- rbind(asvITS_OCG_sigsbst_subspecies, sigASVs_subspecies_t)#attach significance rankings
# asvITS_OCG_sigsbst_subspecies_t <- as.data.frame(t(asvITS_OCG_sigsbst_subspecies))#transpose data frame to prepare for sorting
# asvITS_OCG_sigsbst_subspecies_t$sig_rank <- as.numeric(asvITS_OCG_sigsbst_subspecies_t$sig_rank) #transform significance rankings to numeric (else 10 comes before 2)
# asvITS_OCG_sigsbst_subspecies_t <- asvITS_OCG_sigsbst_subspecies_t[order(asvITS_OCG_sigsbst_subspecies_t$sig_rank),] #sort by significance ranking
# asvITS_OCG_sigsbst_subspecies_t <- subset(asvITS_OCG_sigsbst_subspecies_t, select=-c(sig_rank)) #remove rank column
# asvITS_OCG_sigsbst_subspecies_t <- as.data.frame(t(asvITS_OCG_sigsbst_subspecies_t))#transpose back and convert to data frame
# 
# #Build objects for plotting
# sbsplot <- data.frame(asvITS_OCG_sigsbst_subspecies_t[,1:3], "subspecies" = mdITS_OCG_sbst$Subspecies, check.names = FALSE)
# sbsplot[,1:3] <- lapply(sbsplot[,1:3], function(x) as.numeric(as.character(x)))
# sbsplot[,1:3] <- log(sbsplot[,1:3]+1)
# sbsplot_sub <- data.frame(sample=rownames(sbsplot),sbsplot, check.names = F)
# sbsplot_sublong <- melt(sbsplot_sub)
# 
# ANCOM_subspecies$W.taxa
# 
# sbsplot_sublong$subspecies <-  factor(sbsplot_sublong$subspecies, levels = c("T", "V", "W"))
# 
# #SUBSPECIES FIGURE TOP ANCOM ASVs
# ggplot(sbsplot_sublong, aes(y = value, x = subspecies, color=variable))+
#   geom_boxplot(outlier.shape = NA) + 
#   geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
#   scale_colour_manual('Species',labels=c('Dothidea sambuci', 'Endoconidioma populi','Parastagonospora novozelandica'),values=c("darkkhaki", "palevioletred", "darkcyan"))+
#   ylab("Log  rel. abundance") + xlab("Subspecies") + theme_classic() #Highest in tridentata and then two are very low in vaseyana and one low in wyomingensis. Look it up in the taxonomy in csv and then BLAST it. Dothidea sambuci=1bfe6883c827ae962725d53d567780fc (97.95% identity), 89728a516cec30e6a6c3b1e428ca5f40= Endoconidioma populi (99.60%), 3518c430b3f2d2c3844aee6b7d76990a=	Parastagonospora novozelandica (95.62%)
# 
# print(taxITS.OCG_sigsbst_subspecies)
# 
# ###SIGNIFICANCE BY SUBSPECIES - Default ANCOM sorting 
# #Run ANCOM, specify variable
# ANCOM_subspecies_default <- ANCOM.main(asvITS_OCG_sbst,mdITS_OCG_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)
# #Create objects of significant ASVs
# sigASVs_subspecies_default <- subset(ANCOM_subspecies_default$W.taxa, ANCOM_subspecies_default$W.taxa$W_stat > 0)[,1]
# sigASVs_subspecies_default <- as.data.frame(sigASVs_subspecies_default) #convert significant ASV return to data frame
# row.names(sigASVs_subspecies_default) <- sigASVs_subspecies_default[1:9,1] #change row names to ASV sequence
# 
# #Subset original data by significant ASVs
# asvITS_OCG_sigsbst_subspecies_default <-  t(subset(t(asvITS_OCG_sbst), colnames(asvITS_OCG_sbst) %in% row.names(sigASVs_subspecies_default))) #subset of ASV counts to only significant taxa
# taxITS.OCG_sigsbst_subspecies_default <-  t(subset(t(taxITS.OCG_sbst), colnames(taxITS.OCG_sbst) %in% row.names(sigASVs_subspecies_default))) #subset of taxonomy to only significant taxa
# 
# #Build objects for plotting
# sbsdfplot <- data.frame(asvITS_OCG_sigsbst_subspecies_default[,1:9], "subspecies" = mdITS_OCG_sbst$Subspecies, check.names = FALSE)
# sbsdfplot[,1:9] <- lapply(sbsdfplot[,1:9], function(x) as.numeric(as.character(x)))
# sbsdfplot[,1:9] <- log(sbsdfplot[,1:9]+1)
# 
# sbsdfplot_sub <- data.frame(sample=rownames(sbsdfplot),sbsdfplot, check.names = F)
# sbsdfplot_sublong <- melt(sbsdfplot_sub)
# 
# sbsdfplot_sublong$subspecies <-  factor(sbsdfplot_sublong$subspecies, levels = c("T", "V", "W"))
# 
# ### Figure - OCG top ANCOM ASVs by Subspecies with default sorting 
# ggplot(sbsdfplot_sublong, aes(y = value, x = subspecies, color=variable))+
#   geom_boxplot(outlier.shape = NA) + 
#   geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
#   ylab("Log  rel. abundance") + xlab("Subspecies") + theme_classic()
# 
# print(taxITS.OCG_sigsbst_subspecies_default)
# print(taxITS.OCG_sigsbst_subspecies_default)
# 
# #PLOIDY - no sig differences found
# #Run ANCOM, specify variable
# ANCOM_ploidy <- ANCOM.main(asvITS_OCG_sbst,mdITS_OCG_sbst,F,F,"Ploidy",NULL,NULL,F,NULL,2,.05,.9)
# #Create objects of significant ASVs
# sigASVs_ploidy <- subset(ANCOM_ploidy$W.taxa, ANCOM_ploidy$W.taxa$W_stat > 0)[,1] #this is empty, as there are no significant differences found
# # sigASVs_ploidy <- as.data.frame(sigASVs_ploidy)
# # row.names(sigASVs_ploidy) <- sigASVs_ploidy[1:10,1]
# #sigASVs_ploidy[,1] <- c(1:10)
# #sigASVs_ploidy_t <- t(sigASVs_ploidy)
# #sigASVs_ploidy_t <- as.data.frame(sigASVs_ploidy_t)
# #colnames(sigASVs_ploidy_t) <- as.character(colnames(sigASVs_ploidy_t))
# #print(colnames(sigASVs_ploidy_t))
# #sigASVs_ploidy_t <- sigASVs_ploidy_t[,order(colnames(sigASVs_ploidy_t))]
# #rownames(sigASVs_ploidy_t) <- c("sig_rank")
# 
# write.csv(ANCOM_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_ploidy.csv")
# 
# #Subset data by sig ASVs
# #asvA_sigsbst_ploidy <-  t(subset(t(asvA_sbst), colnames(asvA_sbst) %in% row.names(sigASVs_ploidy)))
# #taxA_sigsbst_ploidy <-  t(subset(t(taxA_sbst), colnames(taxA_sbst) %in% row.names(sigASVs_ploidy)))
# 
# #asvA_sigsbst_ploidy <- asvA_sigsbst_ploidy[,order(colnames(asvA_sigsbst_ploidy))]
# #colnames(sigASVs_ploidy_t) == colnames(asvA_sigsbst_ploidy)
# #asvA_sigsbst_ploidy <- rbind(asvA_sigsbst_ploidy, sigASVs_ploidy_t)
# #asvA_sigsbst_ploidy_t <- as.data.frame(t(asvA_sigsbst_ploidy))
# #asvA_sigsbst_ploidy_t$sig_rank <- as.numeric(asvA_sigsbst_ploidy_t$sig_rank)
# #asvA_sigsbst_ploidy_t <- asvA_sigsbst_ploidy_t[order(asvA_sigsbst_ploidy_t$sig_rank),]
# #asvA_sigsbst_ploidy_t <- subset(asvA_sigsbst_ploidy_t, select=-c(sig_rank))
# #asvA_sigsbst_ploidy <- as.data.frame(t(asvA_sigsbst_ploidy_t))
# 
# #plplot <- data.frame(asvA_sigsbst_ploidy[,1:10], "ploidy" = metaA_sbst$Ploidy, check.names = FALSE)
# #plplot[,1:10] <- lapply(plplot[,1:10], function(x) as.numeric(as.character(x)))
# #plplot[,1:10] <- log(plplot[,1:10]+1)
# 
# #plplot_sub <- data.frame(sample=rownames(plplot),plplot, check.names = F)
# #plplot_sublong <- melt(plplot_sub)
# 
# #plplot_sublong$ploidy <-  factor(plplot_sublong$ploidy, levels = c("2n", "4n"))
# 
# ### (DNE) Figure - OCG top ANCOM ASVs by ploidy 
# # ggplot(plplot_sublong, aes(y = value, x = ploidy, color=variable))+
# #   geom_boxplot(outlier.shape = NA) + 
# #   geom_point(position=position_dodge(width=0.75), aes(group=variable, shape = variable), alpha =.4) +
# #   ylab("Log  rel. abundance") + xlab("Ploidy") + theme_classic()
# # 
# # #print(taxA_sigsbst_ploidy)
# 
# ##SIGNIFICANCE BY SUBSPECIES PLOIDY
# #Run ANCOM, specify variable
# ANCOM_subsp_ploidy <- ANCOM.main(asvITS_OCG_sbst,mdITS_OCG_sbst,F,F,"Subsp_ploidy",NULL,NULL,F,NULL,2,.05,.9)
# #Create objects of significant ASVs
# sigASVs_subsp_ploidy <- subset(ANCOM_subsp_ploidy$W.taxa, ANCOM_subsp_ploidy$W.taxa$W_stat > 0)[,1]
# sigASVs_subsp_ploidy <- as.data.frame(sigASVs_subsp_ploidy)
# row.names(sigASVs_subsp_ploidy) <- sigASVs_subsp_ploidy[1:9,1]
# sigASVs_subsp_ploidy[,1] <- c(1:9)
# sigASVs_subsp_ploidy_t <- t(sigASVs_subsp_ploidy)
# sigASVs_subsp_ploidy_t <- as.data.frame(sigASVs_subsp_ploidy_t)
# colnames(sigASVs_subsp_ploidy_t) <- as.character(colnames(sigASVs_subsp_ploidy_t))
# print(colnames(sigASVs_subsp_ploidy_t))
# sigASVs_subsp_ploidy_t <- sigASVs_subsp_ploidy_t[,order(colnames(sigASVs_subsp_ploidy_t))]
# rownames(sigASVs_subsp_ploidy_t) <- c("sig_rank")
# 
# write.csv(ANCOM_subsp_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_subsp_ploidy.csv")
# 
# #Subset data by sig ASVs
# asvITS_OCG_sigsbst_subsp_ploidy <-  t(subset(t(asvITS_OCG_sbst), colnames(asvITS_OCG_sbst) %in% row.names(sigASVs_subsp_ploidy)))
# taxITS.OCG_sigsbst_subsp_ploidy <-  t(subset(t(taxITS.OCG_sbst), colnames(taxITS.OCG_sbst) %in% row.names(sigASVs_subsp_ploidy)))
# 
# asvITS_OCG_sigsbst_subsp_ploidy <- asvITS_OCG_sigsbst_subsp_ploidy[,order(colnames(asvITS_OCG_sigsbst_subsp_ploidy))]
# colnames(sigASVs_subsp_ploidy_t) == colnames(asvITS_OCG_sigsbst_subsp_ploidy) #TRUE
# asvITS_OCG_sigsbst_subsp_ploidy <- rbind(asvITS_OCG_sigsbst_subsp_ploidy, sigASVs_subsp_ploidy_t)
# asvITS_OCG_sigsbst_subsp_ploidy_t <- as.data.frame(t(asvITS_OCG_sigsbst_subsp_ploidy))
# asvITS_OCG_sigsbst_subsp_ploidy_t$sig_rank <- as.numeric(asvITS_OCG_sigsbst_subsp_ploidy_t$sig_rank)
# asvITS_OCG_sigsbst_subsp_ploidy_t <- asvITS_OCG_sigsbst_subsp_ploidy_t[order(asvITS_OCG_sigsbst_subsp_ploidy_t$sig_rank),]
# asvITS_OCG_sigsbst_subsp_ploidy_t <- subset(asvITS_OCG_sigsbst_subsp_ploidy_t, select=-c(sig_rank))
# asvITS_OCG_sigsbst_subsp_ploidy <- as.data.frame(t(asvITS_OCG_sigsbst_subsp_ploidy_t))
# 
# sbsplplot <- data.frame(asvITS_OCG_sigsbst_subsp_ploidy[,1:9], "subsp_ploidy" = mdITS_OCG_sbst$Subsp_ploidy, check.names = FALSE) 
# 
# sbsplplot[,1:9] <- lapply(sbsplplot[,1:9], function(x) as.numeric(as.character(x)))
# sbsplplot[,1:9] <- log(sbsplplot[,1:9]+1)
# 
# sbsplplot_sub <- data.frame(sample=rownames(sbsplplot),sbsplplot, check.names = F)
# sbsplplot_sublong <- melt(sbsplplot_sub)
# 
# sbsplplot_sublong$subsp_ploidy <-  factor(sbsplplot_sublong$subsp_ploidy, levels = c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))
# 
# ### FIGURE FOR SUBSPECIES PLOIDY
# ggplot(sbsplplot_sublong, aes(y = value, x = subsp_ploidy, color=variable))+
#   geom_boxplot(outlier.shape = NA) + 
#   geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
#   ylab("Log  rel. abundance") + xlab("Subsp_ploidy") + theme_classic()
# 
# print(taxITS.OCG_sigsbst_subsp_ploidy)
# 
# ##SIGNIFICANCE BY YEAR
# #Run ANCOM, specify variable
# ANCOM_year <- ANCOM.main(asvITS_OCG_sbst,mdITS_OCG_sbst,F,F,"Year",NULL,NULL,F,NULL,2,.05,.9)
# #Create objects of significant ASVs
# sigASVs_year <- subset(ANCOM_year$W.taxa, ANCOM_year$W.taxa$W_stat > 0)[,1]
# sigASVs_year <- as.data.frame(sigASVs_year)
# row.names(sigASVs_year) <- sigASVs_year[1:10,1]
# sigASVs_year[,1] <- c(1:10)
# sigASVs_year_t <- t(sigASVs_year)
# sigASVs_year_t <- as.data.frame(sigASVs_year_t)
# colnames(sigASVs_year_t) <- as.character(colnames(sigASVs_year_t))
# print(colnames(sigASVs_year_t))
# sigASVs_year_t <- sigASVs_year_t[,order(colnames(sigASVs_year_t))]
# rownames(sigASVs_year_t) <- c("sig_rank")
# 
# write.csv(ANCOM_year$W.taxa, file = "data_csv/ANCOM/ANCOM_Year.csv")
# 
# #Subset data by sig ASVs
# asvITS_OCG_sigsbst_year <-  t(subset(t(asvITS_OCG_sbst), colnames(asvITS_OCG_sbst) %in% row.names(sigASVs_year)))
# taxITS.OCG_sigsbst_year <-  t(subset(t(taxITS.OCG_sbst), colnames(taxITS.OCG_sbst) %in% row.names(sigASVs_year)))
# 
# asvITS_OCG_sigsbst_year <- asvITS_OCG_sigsbst_year[,order(colnames(asvITS_OCG_sigsbst_year))]
# colnames(sigASVs_year_t) == colnames(asvITS_OCG_sigsbst_year) #TRUE
# asvITS_OCG_sigsbst_year <- rbind(asvITS_OCG_sigsbst_year, sigASVs_year_t)
# asvITS_OCG_sigsbst_year_t <- as.data.frame(t(asvITS_OCG_sigsbst_year))
# asvITS_OCG_sigsbst_year_t$sig_rank <- as.numeric(asvITS_OCG_sigsbst_year_t$sig_rank)
# asvITS_OCG_sigsbst_year_t <- asvITS_OCG_sigsbst_year_t[order(asvITS_OCG_sigsbst_year_t$sig_rank),]
# asvITS_OCG_sigsbst_year_t <- subset(asvITS_OCG_sigsbst_year_t, select=-c(sig_rank))
# asvITS_OCG_sigsbst_year <- as.data.frame(t(asvITS_OCG_sigsbst_year_t))
# 
# yrplot <- data.frame(asvITS_OCG_sigsbst_year[,1:10], "year" = mdITS_OCG_sbst$Year, check.names = FALSE)
# yrplot[,1:10] <- lapply(yrplot[,1:10], function(x) as.numeric(as.character(x)))
# yrplot[,1:10] <- log(yrplot[,1:10]+1)
# 
# yrplot_sub <- data.frame(sample=rownames(yrplot),yrplot, check.names = F)
# yrplot_sublong <- melt(yrplot_sub)
# 
# yrplot_sublong$year <-  factor(yrplot_sublong$year, levels = c("2012", "2021"))
# 
# ### FIGURE FOR YEAR
# ggplot(yrplot_sublong, aes(y = value, x = year, color=variable))+
#   geom_boxplot(outlier.shape = NA) + 
#   geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
#   ylab("Log  rel. abundance") + xlab("Year") + theme_classic()
# 
# print(taxITS.OCG_sigsbst_year)
# 
# ## ANCOM BC 2 ####
# # Create phyloseq object
# asvITS_OCG_t <- t(asvITS_OCG)
# asvITS_OCG_t <- as.data.frame(asvITS_OCG_t)
# taxITS.OCG <- subset(taxITS.OCG, row.names(taxITS.OCG) %in% rownames(asvITS_OCG_t))
# taxITS.OCG <- taxITS.OCG[order(row.names(taxITS.OCG)),] # order samples alphabetically
# row.names(taxITS.OCG) == row.names(asvITS_OCG_t) #TRUE
# colnames(asvITS_OCG_t) == row.names(mdITS_OCG) #TRUE
# otu_asv_table <- otu_table(as.matrix(asvITS_OCG_t), taxa_are_rows = TRUE)
# tax_ITS_table <- tax_table(as.matrix(taxITS.OCG))
# md_ITS_phy <- sample_data(mdITS_OCG) #make metadata into sample data to create phyloseq
# physeqITS <- merge_phyloseq(phyloseq(otu_asv_table),md_ITS_phy,tax_ITS_table) #create phyloseq object for ANCOM BC
# physeqITS #86 taxa and 164 samples
# 
# #MAKE TSE
# tse.ITS = mia::makeTreeSummarizedExperimentFromPhyloseq(physeqITS)
# tse.ITS$Subspecies <- factor(tse.ITS$Subspecies, levels=c("T", "W", "V"))
# tse.ITS$Ploidy <- factor(tse.ITS$Ploidy, levels=c("2n", "4n"))
# tse.ITS$Subsp_ploidy <- factor(tse.ITS$Subsp_ploidy, levels=c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))
# tse.ITS$Year <- factor(tse.ITS$Year, levels =c("2012", "2021"))
# 
# #Model for Subspecies
# result_ITS_sub <- ancombc2( 
#   data = tse.ITS, assay_name = "counts", tax_level = NULL,
#   fix_formula = "Subspecies + Ploidy + Year",
#   p_adj_method = "fdr", pseudo_sens = TRUE,
#   group = "Subspecies", #change variable of interest
#   alpha = 0.05, verbose = TRUE,
#   global = TRUE, prv_cut = 0.02
# )
# 
# summary(tse.ITS$Year)
# 
# ancomITSres_df.sub <- result_ITS_sub$res
# 
# #Model for year
# result_ITS_yr <- ancombc2( 
#   data = tse.ITS, assay_name = "counts", tax_level = NULL,
#   fix_formula = "Subsp_ploidy + Year",
#   p_adj_method = "fdr", pseudo_sens = TRUE,
#   group = "Year", #change variable of interest
#   alpha = 0.05, verbose = TRUE,
#   global = TRUE, prv_cut = 0.1
# )
# 
# ancomITSres_df.yr <- result_ITS_yr$res
# 
# #iter_control = list(tol = 1e-5, max_iter = 200, 
#                     #verbose = FALSE)
# 
# #METACODER ####
# #Create taxmap object for use with Metacoder functions 
# ##Create data frame matching taxonomy information with ASV sequences ###
# row.names(asvITS_OCG) == row.names(taxITS.OCG)# sanity check:TRUE
# totMC_OCG <- asvITS_OCG # new data frame to hold combined asv + taxonomy information
# totMC_OCG$Taxon <- taxITS.OCG$Taxon # append taxonomy information to ASV data frame
# row.names(taxITS.OCG) == row.names(totMC_OCG)# sanity check for ASVs: TRUE
# taxITS.OCG$Taxon == totMC_OCG$Taxon # sanity check for taxonomy: TRUE
# 
# ## Create taxmap object ####
# obj <- parse_tax_data(totMC_OCG,
#                       class_cols = "Taxon",                   # name of column that contains input taxon data
#                       class_sep = ";",                        # character that separates taxon data
#                       class_regex = "^(.+)__(.+)$",           # regex to identify taxon entries
#                       class_key = c(tax_rank = "info",        # this is the key that labels each column pulled from the regex, since we had two sections for each identifier we need two columns
#                                     tax_name = "taxon_name"))
# 
# print(obj) #ERROR                   
# print(obj$data$tax_data)
# print(obj$data$class_data)
# obj$data$class_data <- NULL     # class_data is repetitive/unnecessary
# names(obj$data) <- "asv_counts" # rename "data" to something more relevant
# print(obj)                      
# 
# heat_tree(obj, # very basic heat tree showing the overall composition of Orchard Common Garden samples
#           node_label = taxon_names, 
#           node_size = n_obs,        
#           node_color = n_obs)
# 
# #DIFFERENTIAL HEAT TREES FOR OCG
# ##Adjust taxmap so we can investigate taxa instead of ASVs #
# #Now we need to calculate abundances based on taxon not ASV
# obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "asv_counts",
#                                                   cols = row.names(mdITS_OCG))
# print(obj$data$tax_abund)
# ##**Notice** ###
# # Metacoder did not play well with separate naming conventions. As such, all heatmaps are each written into the same places: "diff_table" and "diff_heattree_color"
# # To look at different maps, run code starting from the beginning of each section: "Heatmap matrix by _____"
# 
# ##HEATMAP MATRIX BY SUBSPECIES
# obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund",
#                                                  cols = row.names(mdITS_OCG), # What columns of sample data to use (sample ID is stored in row names of metadata)
#                                                  groups = mdITS_OCG$Subspecies) # What category each sample is assigned to from metadata
# 
# print(obj$data$diff_table)
# obj <- mutate_obs(obj, "diff_table",                                                                
#                   wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
# 
# #look at the p-values and see if there is any significance
# range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
# # [1] 0.000186618 0.99738256
# # lower range is significant
# 
# ## Focus only on significant taxa
# obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0
# 
# set.seed(1)
# diff_heattree_color <- metacoder::heat_tree_matrix(obj, data = "diff_table",
#                                                    node_size = n_obs, 
#                                                    node_label = taxon_names,
#                                                    node_color = log2_median_ratio,
#                                                    node_color_range = diverging_palette(),
#                                                    node_color_trans = "linear", 
#                                                    node_color_interval = c(-3, 3), 
#                                                    edge_color_interval = c(-3, 3), 
#                                                    node_size_axis_label = "Number of ASVs",
#                                                    node_color_axis_label = "Log2 ratio median proportions",
#                                                    layout = "davidson-harel", 
#                                                    initial_layout = "reingold-tilford")
# 
# ## This plot takes a while to load
# print(diff_heattree_color) ## Show taxonomic heat tree 
# 
# 
# 
# ##HEATMAP MATRIX BY PLOIDY
# obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund",
#                                                  cols = row.names(mdITS_OCG), # What columns of sample data to use (sample ID is stored in row names of metadata)
#                                                  groups = mdITS_OCG$Ploidy) # What category each sample is assigned to from metadata
# 
# print(obj$data$diff_table)
# obj <- mutate_obs(obj, "diff_table",                                                                
#                   wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
# 
# #lets look at the p-values and see if there is any significance
# range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
# # [1] 0.5037569 0.9958621
# # not significant
# 
# ## Focus only on significant taxa
# obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0
# 
# 
# set.seed(1)
# diff_heattree_color <- metacoder::heat_tree_matrix(obj, data = "diff_table",
#                                                    node_size = n_obs, 
#                                                    node_label = taxon_names,
#                                                    node_color = log2_median_ratio,
#                                                    node_color_range = diverging_palette(),
#                                                    node_color_trans = "linear", 
#                                                    node_color_interval = c(-3, 3), 
#                                                    edge_color_interval = c(-3, 3), 
#                                                    node_size_axis_label = "Number of ASVs",
#                                                    node_color_axis_label = "Log2 ratio median proportions",
#                                                    layout = "davidson-harel", 
#                                                    initial_layout = "reingold-tilford")
# 
# ## This plot takes a while to load
# print(diff_heattree_color) ## Show taxonomic heat tree 
# 
# #HEATMAP MATRIX BY SUBSPECIES PLOIDY
# obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund",
#                                                  cols = row.names(mdITS_OCG), # What columns of sample data to use (sample ID is stored in row names of metadata)
#                                                  groups = mdITS_OCG$Subsp_ploidy) # What category each sample is assigned to from metadata
# 
# print(obj$data$diff_table)
# obj <- mutate_obs(obj, "diff_table",                                                               
#                   wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
# 
# #lets look at the p-values and see if there is any significance
# range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
# # [1] 0.04266539 1.00000000
# # the lower range is significant
# 
# ## Focus only on significant taxa
# obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0
# 
# print(obj$data$diff_table)
# 
# set.seed(1)
# diff_heattree_color <- metacoder::heat_tree_matrix(obj, data = "diff_table",
#                                                    node_size = n_obs, 
#                                                    node_label = taxon_names,
#                                                    node_color = log2_median_ratio,
#                                                    node_color_range = diverging_palette(),
#                                                    node_color_trans = "linear", 
#                                                    node_color_interval = c(-3, 3), 
#                                                    edge_color_interval = c(-3, 3), 
#                                                    node_size_axis_label = "Number of ASVs",
#                                                    node_color_axis_label = "Log2 ratio median proportions",
#                                                    layout = "davidson-harel", 
#                                                    initial_layout = "reingold-tilford")
# 
# ## This plot takes a while to load
# print(diff_heattree_color) ## Show taxonomic heat tree
# 
# ##HEATMAP BY YEAR
# obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund",
#                                                  cols = row.names(mdITS_OCG), # What columns of sample data to use (sample ID is stored in row names of metadata)
#                                                  groups = mdITS_OCG$Year) # What category each sample is assigned to from metadata
# 
# print(obj$data$diff_table)
# obj <- mutate_obs(obj, "diff_table",                                                                
#                   wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
# 
# #lets look at the p-values and see if there is any significance
# range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
# # [1] 2.288654e-06 9.026667e-01
# # the lower range is significant
# 
# ## Focus only on significant taxa
# obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0
# 
# 
# set.seed(1)
# diff_heattree_color <- metacoder::heat_tree(obj,  # heat tree code slightly different here on account of only having two variables, use heat_tree, exclude data argument
#                                             node_size = n_obs, 
#                                             node_label = taxon_names,
#                                             node_color = log2_median_ratio,
#                                             node_color_range = diverging_palette(),
#                                             node_color_trans = "linear", 
#                                             node_color_interval = c(-3, 3), 
#                                             edge_color_interval = c(-3, 3), 
#                                             node_size_axis_label = "Number of ASVs",
#                                             node_color_axis_label = "Log2 ratio median proportions",
#                                             layout = "davidson-harel", 
#                                             initial_layout = "reingold-tilford")
# 
# ## This plot takes a while to load
# print(diff_heattree_color) ## Show taxonomic heat tree
# #after cross-referencing with ANCOM results, this heat map shows changes from 2012 -> 2021; i.e. positive ratio indicates more sequences in 2021
# 


# 2021 GC versus 2021 asv procrustes 
OCG_GC_2021_pro.nmds
asv_2021_pro.nmds

asv_GC_2021.pro <- protest(OCG_GC_2021_pro.nmds, asv_2021_pro.nmds, symmetric=T) 
asv_GC_2021.pro ## Correlation in a symmetric Procrustes rotation: 0.1994, Significance: 0.05
summary(asv_GC_2021.pro) 
plot(asv_GC_2021.pro)

## Re-plotting the procrustes 
asv_GC_2021.pro_prodat <- as.data.frame(asv_GC_2021.pro$X)
asv_GC_2021.pro_prodat <- cbind(asv_GC_2021.pro_prodat,asv_GC_2021.pro$Yrot)
colnames(asv_GC_2021.pro_prodat)[colnames(asv_GC_2021.pro_prodat)=="1"] <- "Xend"
colnames(asv_GC_2021.pro_prodat)[colnames(asv_GC_2021.pro_prodat)=="2"] <- "Yend"

ggplot() + 
  geom_segment(data=asv_GC_2021.pro_prodat, mapping=aes(x=NMDS1, y=NMDS2, xend=Xend, yend=Yend), size=0.6, color="lightgrey") + 
  geom_point(data=asv_GC_2021.pro_prodat, mapping=aes(x=NMDS1, y=NMDS2), size=2, shape=19, color = "maroon") +
  geom_point(data=asv_GC_2021.pro_prodat, mapping=aes(x=Xend, y=Yend), size=2, shape=17, color = "lightseagreen") +
  labs(x="Procrustes axis 1", y="Procrustes axis 2") +
  theme_classic() 

# dbRDA for all fungi ~ predictors ####
#coordinate data, subspecies:cytoptype data, GC data, and year

# make sure row names match in fungal data and metadata data
rownames(asvITS_OCG.r) == rownames(mdITS_OCG) # sanity check: TRUE

# fungal_dist <- vegdist(asvITS_OCG.r, method = "bray")

# Read in GC data 
OCG_GC <- read.csv("data_csv/OCG_GC_thresholded.csv", head = T, row.names = 1, check.names = F, stringsAsFactors = F)
OCG_GC[is.na(OCG_GC)] <- 0

# # Read in LCMS data 
# OCG_LCMS <- read.csv("data_csv/OCG_LCMS_thresholded.csv", head = T, row.names = 1, check.names = F, stringsAsFactors = F)
# OCG_LCMS[is.na(OCG_LCMS)] <- 0

# subset to have md data only match LCMS and GC datasets
mdITS_OCG_LCMS <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(OCG_LCMS)) # 73
mdITS_OCG_GC_LCMS <- subset(mdITS_OCG_LCMS, row.names(mdITS_OCG_LCMS) %in% row.names(OCG_GC)) # 69
OCG_GC_ITS <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(mdITS_OCG_GC_LCMS))
OCG_LCMS_ITS <- subset(OCG_LCMS, row.names(OCG_LCMS) %in% row.names(mdITS_OCG_GC_LCMS)) 

# subset to match with all asv data
asvITS_OCG_GC_LCMS <- subset(asvITS_OCG, row.names(asvITS_OCG) %in% row.names(mdITS_OCG_GC_LCMS)) #151

# order the row names
OCG_GC_ITS <- OCG_GC_ITS[order(row.names(OCG_GC_ITS)), ]
OCG_LCMS_ITS <- OCG_LCMS_ITS[order(row.names(OCG_LCMS_ITS)), ]
mdITS_OCG_GC_LCMS <- mdITS_OCG_GC_LCMS[order(row.names(mdITS_OCG_GC_LCMS)), ]
asvITS_OCG_GC_LCMS <- asvITS_OCG_GC_LCMS[order(row.names(asvITS_OCG_GC_LCMS)), ]

row.names(OCG_GC_ITS) == row.names(mdITS_OCG_GC_LCMS) # sanity check: TRUE
row.names(OCG_GC_ITS) == row.names(asvITS_OCG_GC_LCMS) # sanity check: TRUE
row.names(OCG_LCMS_ITS) == row.names(mdITS_OCG_GC_LCMS) # sanity check: TRUE
row.names(OCG_LCMS_ITS) == row.names(asvITS_OCG_GC_LCMS) # sanity check: TRUE

# run PCA on GC data
OCG_GC_ITS_no0 <- OCG_GC_ITS[, colSums(OCG_GC_ITS != 0) > 0]
scaled_OCG_GC <- scale(OCG_GC_ITS_no0)

pca_gc <- prcomp(scaled_OCG_GC, center = TRUE)
summary(pca_gc) # 26.22 0.9223

loadings_gc <- pca_gc$rotation
pc1_loadings_gc <- loadings_gc[, 1]  # loadings for PC1

# Sort by absolute value of loading
important_pc1 <- sort(abs(pc1_loadings_gc), decreasing = TRUE)

# View top 10 most influential compounds (both positive and negatively influencing)
head(important_pc1, 10) # 3, 29, 10, 6, 65, 57, 34, 22, 53, 27

barplot(sort(pc1_loadings_gc, decreasing = TRUE)[1:10],
        las = 2,
        main = "Top 10 compounds contributing to PC1",
        ylab = "Loading value",
        col = "skyblue")

gc_df <- data.frame(Compound = names(pc1_loadings_gc),
                 Loading = pc1_loadings_gc)

ggplot(gc_df, aes(x = reorder(Compound, Loading), y = Loading, fill = Loading > 0)) +
  geom_col() +
  coord_flip() +
  labs(title = "Compound contributions to PC1",
       x = "Compound",
       y = "Loading") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"))

# save GC PCs in a dataframe
pca_scores_gc <- as.data.frame(pca_gc$x[,])

# run PCA on LCMS data
OCG_LCMS_ITS_no0 <- OCG_LCMS_ITS[, colSums(OCG_LCMS_ITS != 0) > 0]
scaled_OCG_LCMS <- scale(OCG_LCMS_ITS_no0)

pca_lcms <- prcomp(scaled_OCG_LCMS, center = TRUE)
summary(pca_lcms) # 15.14 9.47

# save GC PCs in a dataframe
pca_scores_lcms <- as.data.frame(pca_lcms$x[,])

#create fungal distance matrix
fungal_dist <- vegdist(asvITS_OCG_GC_LCMS, method = "bray")

# # Procrustes between first two PCs of GC and LCMS
# GC_scores <- pca_gc$x[, 1:2]
# LCMS_scores <- pca_lcms$x[, 1:2]
# 
# pro_gc_lcms <- protest(GC_scores, LCMS_scores, permutations = 999) # sig

# Perform dbRDA
str(mdITS_OCG_GC_LCMS)
dbrda_gc_lcms_model <- vegan::dbrda(fungal_dist ~  mdITS_OCG_GC_LCMS$x + mdITS_OCG_GC_LCMS$y + pca_scores_gc$PC1 + pca_scores_gc$PC2 + pca_scores_gc$PC3 + pca_scores_lcms$PC1 + pca_scores_lcms$PC2 + pca_scores_lcms$PC3 + Condition(mdITS_OCG_GC_LCMS$Plant), data = mdITS_OCG_GC_LCMS) 
anova(dbrda_gc_lcms_model, by = "term", permutations = 999)

# plot(dbrda_model, type = "n")  # plot empty frame
# points(dbrda_model, display = "sites", col = "hotpink", pch = 19)  # sample points
# text(dbrda_model, display = "bp", col = "darkviolet")  # arrows and labels for predictors

# Procrustes for GC and fungal data 2012 ####
asvITS_OCG_GC.r <- subset(asvITS_OCG.r, row.names(asvITS_OCG.r) %in% row.names(OCG_GC)) # 91
OCG_GC_ITS <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(asvITS_OCG_GC.r)) # 91

# order the row names
OCG_GC_ITS <- OCG_GC_ITS[order(row.names(OCG_GC_ITS)), ]
asvITS_OCG_GC.r <- asvITS_OCG_GC.r[order(row.names(asvITS_OCG_GC.r)), ]

scaled_OCG_GC_ITS <- scale(OCG_GC_ITS)
pca_gc_pro <- prcomp(OCG_GC_ITS, center = TRUE)

row.names(OCG_GC_ITS) == row.names(asvITS_OCG_GC.r) # sanity check: TRUE

set.seed(4637)
procrustes_nmds_fungi <- metaMDS(asvITS_OCG_GC.r, distance = "bray", trymax = 500)
pro_gc_fun <- protest(pca_gc_pro, procrustes_nmds_fungi, permutations = 999) # Sum of squares = 0.9702, Correlation in a symmetric Procrustes rotation: 0.1725, Significance: 0.062


# dbRDA for 2012 asv data ####
# read in 2012 GC data
OCG_GC_2012 <- read.csv("data_csv/OCG_GC_2012_thresholded.csv", head = T, row.names = 1, check.names = F, stringsAsFactors = F)
OCG_GC_2012[is.na(OCG_GC_2012)] <- 0

# # read in 2012 LCMS data
# OCG_LCMS_2012 <- read.csv("data_csv/OCG_LCMS_2012_thresholded.csv", head = T, row.names = 1, check.names = F, stringsAsFactors = F)
# OCG_LCMS_2012[is.na(OCG_LCMS_2012)] <- 0

# subset to match with all asv data
asvITS_OCG_GC_2012 <- subset(asvITS.2012, row.names(asvITS.2012) %in% row.names(mdITS_OCG)) # 100
OCG_GC_2012_ITS <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(asvITS_OCG_GC_2012)) # 100
# OCG_LCMS_2012_ITS <- subset(OCG_LCMS_2012, row.names(OCG_LCMS_2012) %in% row.names(asvITS.2012)) # 24
mdITS_OCG_GC_2012 <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(OCG_GC_2012_ITS)) # 100
asvITS_OCG_GC_2012 <- subset(asvITS_OCG_GC_2012, row.names(asvITS_OCG_GC_2012) %in% row.names(mdITS_OCG_GC_2012)) # 100

# order the row names
OCG_GC_2012_ITS <- OCG_GC_2012_ITS[order(row.names(OCG_GC_2012_ITS)), ]
# OCG_LCMS_ITS <- OCG_LCMS_ITS[order(row.names(OCG_LCMS_ITS)), ]
mdITS_OCG_GC_2012 <- mdITS_OCG_GC_2012[order(row.names(mdITS_OCG_GC_2012)), ]
asvITS_OCG_GC_2012 <- asvITS_OCG_GC_2012[order(row.names(asvITS_OCG_GC_2012)), ]

row.names(OCG_GC_2012_ITS) == row.names(asvITS_OCG_GC_2012) # sanity check: TRUE
row.names(OCG_GC_2012_ITS) == row.names(mdITS_OCG_GC_2012) # sanity check: TRUE

# run PCA on GC data
OCG_GC_2012_ITS_no0 <- OCG_GC_2012_ITS[, colSums(OCG_GC_2012_ITS != 0) > 0]
scaled_OCG_GC_2012 <- scale(OCG_GC_2012_ITS_no0)

pca_gc_2012 <- prcomp(scaled_OCG_GC_2012, center = TRUE)
summary(pca_gc_2012) # 20.43 10.27

# save GC PCs in a dataframe
pca_scores_gc_2012 <- as.data.frame(pca_gc_2012$x[,])

# 0.2013 + 0.1027 + 0.08159 + 0.05515 + 0.04985 + 0.04546 + 0.03788 + 0.0351 + 0.0351 + 0.03027  # PC 1 - PC 10 # 0.6744

fungal_dist_2012 <- vegdist(asvITS_OCG_GC_2012, method = "bray")

# Perform dbRDA
dbrda_2012_gc_model <- vegan::dbrda(fungal_dist_2012 ~ mdITS_OCG_GC_2012$x + mdITS_OCG_GC_2012$y + pca_scores_gc_2012$PC1 + pca_scores_gc_2012$PC2 + pca_scores_gc_2012$PC3, data = mdITS_OCG_GC_2012) 

# save(dbrda_model_2012, file = "dbrda_model_2012.rda")

anova(dbrda_2012_gc_model, by = "term", permutations = 999)

# Procrustes for GC and fungal data 2012 ####
asvITS_OCG_GC_2012.r <- subset(asvITS.2012.r, row.names(asvITS.2012.r) %in% row.names(OCG_GC_2012_ITS)) # 91

row.names(OCG_GC_2012_ITS) == row.names(asvITS_OCG_GC_2012.r) # sanity check: TRUE

set.seed(67)
procrustes_nmds_fungi12 <- metaMDS(asvITS_OCG_GC_2012.r, distance = "bray", trymax = 500)
pro_gc_fun12 <- protest(procrustes_nmds_fungi12, pca_gc_2012, permutations = 999) # Sum of squares = 0.9702, Correlation in a symmetric Procrustes rotation: 0.1725, Significance: 0.062

# dbRDA for 2021 asv data vs GC ####
# # Read in thermal data and add it to md
# thermal_data <- read.csv("data_csv/orchard_2019_thermal_data.csv", header = T, check.names = F)

# #subset to just June season
# thermal_data <- thermal_data[thermal_data$Season == "June", ]
# 
# colnames(thermal_data)[colnames(thermal_data) == "ID"] <- "Garden Plant ID"
# # thermal_data <- thermal_data[thermal_data, ]
# 
# #remove rows with NAs from thermal data
# thermal_data <- thermal_data[!is.na(thermal_data$UAS_Thermal), ]
# 
# # subset thermal need to be matched with metadata.
# mdITS_OCG_2021 <- mdITS_OCG[mdITS_OCG$Year == "2021", ]
# mdITS_OCG_2021_thermal <- merge(mdITS_OCG_2021, thermal_data, by = "Garden Plant ID")
# # mdITS_OCG_2021_thermal <- mdITS_OCG_2021_thermal[!duplicated(mdITS_OCG_2021_thermal$`Garden Plant ID`), ] 
# 
# # remove rows with NAs 
# 
# #make description row names
# mdITS_OCG_2021_thermal$Description -> rownames(mdITS_OCG_2021_thermal)

# read in 2021 GC data
OCG_GC_2021 <- read.csv("data_csv/OCG_GC_2021_thresholded.csv", head = T, row.names = 1, check.names = F, stringsAsFactors = F)

# make NAs zero
OCG_GC_2021[is.na(OCG_GC_2021)] <- 0

# Procrustes for GC and fungal data 2021 ####
asvITS_OCG_GC_2021.r <- subset(asvITS.2021.r, row.names(asvITS.2021.r) %in% row.names(OCG_GC_2021_ITS)) # 43
OCG_GC_2021_ITS.r <- subset(OCG_GC_2021_ITS, row.names(OCG_GC_2021_ITS) %in% row.names(asvITS_OCG_GC_2021.r)) # 43
row.names(OCG_GC_2021_ITS.r) == row.names(asvITS_OCG_GC_2021.r) # sanity check: TRUE

# run PCA on GC data
OCG_GC_2021_ITS_no0.r <- OCG_GC_2021_ITS.r[, colSums(OCG_GC_2021_ITS.r != 0) > 0]
scaled_OCG_GC_2021.r <- scale(OCG_GC_2021_ITS_no0.r)

set.seed(95)
procrustes_nmds_fungi21 <- metaMDS(asvITS_OCG_GC_2021.r, distance = "bray", trymax = 500)
pca_gc_2021.r <- prcomp(scaled_OCG_GC_2021.r, center = TRUE)
pro_gc_fun21 <- protest(procrustes_nmds_fungi21, pca_gc_2021.r, permutations = 999)
# 0.9363, 0.2523, 0.077

# ## 2021 LC-MS presence threshold defined ####
# # T_4n
# md_OCG_2021_LCMS_t4n <- subset(md_OCG_2021, md_OCG_2021$Subsp_ploidy=="T_4n") # 16
# OCG_LCMS_2021_t4n <- subset(OCG_LCMS_2021, row.names(OCG_LCMS_2021) %in% row.names(md_OCG_2021_LCMS_t4n)) # 11 of 308
# 
# OCG_LCMS_2021_t4n[OCG_LCMS_2021_t4n == 0] <- NA
# colSums(is.na(OCG_LCMS_2021_t4n)) #checking for null values since they are used for filtering
# 
# #Calculate proportion of NA values that are in each column
# na_proportion_lcms_2021_t4n <- colMeans(is.na(OCG_LCMS_2021_t4n))
# print(na_proportion_lcms_2021_t4n)
# 
# threshold <- 0.8 
# #identify which columns to keep
# columns_to_keep_lcms_2021_t4n <- na_proportion_lcms_2021_t4n <= threshold 
# 
# # Subset dataframe to keep only columns with NA proportion <= threshold
# OCG_LCMS_2021_t4n_subset <- OCG_LCMS_2021_t4n[, columns_to_keep_lcms_2021_t4n] # 16 of 279
# 
# # Replace NA values with zeroes
# OCG_LCMS_2021_t4n_subset[is.na(OCG_LCMS_2021_t4n_subset)] <- 0
# colSums(is.na(OCG_LCMS_2021_t4n_subset))
# 
# # T_2n
# md_OCG_LCMS_2021_t2n <- subset(md_OCG_2021, md_OCG_2021$Subsp_ploidy=="T_2n") # 28
# OCG_LCMS_2021_t2n <- subset(OCG_LCMS_2021, row.names(OCG_LCMS_2021) %in% row.names(md_OCG_LCMS_2021_t2n)) # 24 of 307
# 
# OCG_LCMS_2021_t2n[OCG_LCMS_2021_t2n == 0] <- NA
# colSums(is.na(OCG_LCMS_2021_t2n)) 
# 
# #Calculate proportion of NA values that are in each column
# na_proportion_lcms_2021_t2n <- colMeans(is.na(OCG_LCMS_2021_t2n))
# print(na_proportion_lcms_2021_t2n)
# 
# #identify which columns to keep
# columns_to_keep_lcms_2021_t2n <- na_proportion_lcms_2021_t2n <= threshold 
# 
# # Subset dataframe to keep only columns with NA proportion <= threshold
# OCG_LCMS_2021_t2n_subset <- OCG_LCMS_2021_t2n[, columns_to_keep_lcms_2021_t2n] # 24 of 250
# 
# # Replace NA values with zeroes
# OCG_LCMS_2021_t2n_subset[is.na(OCG_LCMS_2021_t2n_subset)] <- 0
# colSums(is.na(OCG_LCMS_2021_t2n_subset))
# 
# # V_4n
# md_OCG_LCMS_2021_v4n <- subset(md_OCG_2021, md_OCG_2021$Subsp_ploidy=="V_4n") # 5
# OCG_LCMS_2021_v4n <- subset(OCG_LCMS_2021, row.names(OCG_LCMS_2021) %in% row.names(md_OCG_LCMS_2021_v4n)) # 4 of 307
# 
# OCG_LCMS_2021_v4n[OCG_LCMS_2021_v4n == 0] <- NA
# colSums(is.na(OCG_LCMS_2021_v4n))
# 
# #Calculate proportion of NA values that are in each column
# na_proportion_lcms_2021_v4n <- colMeans(is.na(OCG_LCMS_2021_v4n))
# print(na_proportion_lcms_2021_v4n)
# 
# #identify which columns to keep
# columns_to_keep_lcms_2021_v4n <- na_proportion_lcms_2021_v4n <= threshold
# 
# # Subset dataframe to keep only columns with NA proportion <= threshold
# OCG_LCMS_2021_v4n_subset <- OCG_LCMS_2021_v4n[, columns_to_keep_lcms_2021_v4n] # 4 of 290
# 
# # Replace NA values with zeroes
# OCG_LCMS_2021_v4n_subset[is.na(OCG_LCMS_2021_v4n_subset)] <- 0
# colSums(is.na(OCG_LCMS_2021_v4n_subset))
# 
# # # V_2n
# md_OCG_LCMS_2021_v2n <- subset(md_OCG_2021, md_OCG_2021$Subsp_ploidy=="V_2n") # 1
# OCG_LCMS_2021_v2n <- subset(OCG_LCMS_2021, row.names(OCG_LCMS_2021) %in% row.names(md_OCG_LCMS_2021_v2n)) # 1 of 307
# 
# OCG_LCMS_2021_v2n[OCG_LCMS_2021_v2n == 0] <- NA
# colSums(is.na(OCG_LCMS_2021_v2n))
# 
# # Calculate proportion of NA values that are in each column
# na_proportion_lcms_2021_v2n <- colMeans(is.na(OCG_LCMS_2021_v2n))
# print(na_proportion_lcms_2021_v2n)
# 
# #identify which columns to keep
# columns_to_keep_lcms_2021_v2n <- na_proportion_lcms_2021_v2n <= threshold
# 
# # Subset dataframe to keep only columns with NA proportion <= threshold
# OCG_LCMS_2021_v2n_subset <- OCG_LCMS_2021_v2n[, columns_to_keep_lcms_2021_v2n] # 1 of 234
# 
# # Replace NA values with zeroes
# OCG_LCMS_2021_v2n_subset[is.na(OCG_LCMS_2021_v2n_subset)] <- 0
# colSums(is.na(OCG_LCMS_2021_v2n_subset))
# 
# # W_4n
# md_OCG_LCMS_2021_w4n <- subset(md_OCG_2021, md_OCG_2021$Subsp_ploidy=="W_4n") # 26
# OCG_LCMS_2021_w4n <- subset(OCG_LCMS_2021, row.names(OCG_LCMS_2021) %in% row.names(md_OCG_LCMS_2021_w4n)) # 26 of 307
# 
# OCG_LCMS_2021_w4n[OCG_LCMS_2021_w4n == 0] <- NA
# colSums(is.na(OCG_LCMS_2021_w4n))
# 
# #Calculate proportion of NA values that are in each column
# na_proportion_lcms_2021_w4n <- colMeans(is.na(OCG_LCMS_2021_w4n))
# print(na_proportion_lcms_2021_w4n)
# 
# #identify which columns to keep
# columns_to_keep_lcms_2021_w4n <- na_proportion_lcms_2021_w4n <= threshold
# 
# # Subset dataframe to keep only columns with NA proportion <= threshold
# OCG_LCMS_2021_w4n_subset <- OCG_LCMS_2021_w4n[, columns_to_keep_lcms_2021_w4n] # 26 of 279
# 
# # Replace NA values with zeroes
# OCG_LCMS_2021_w4n_subset[is.na(OCG_LCMS_2021_w4n_subset)] <- 0
# colSums(is.na(OCG_LCMS_2021_w4n_subset))
# 
# # Combine all 5 subspecies groups into one LCMS dataset
# OCG_LCMS_2021_subset <- bind_rows(
#   OCG_LCMS_2021_t4n_subset,
#   OCG_LCMS_2021_t2n_subset,
#   OCG_LCMS_2021_v2n_subset,
#   OCG_LCMS_2021_v4n_subset,
#   OCG_LCMS_2021_w4n_subset
# ) # 71 of 306

# write LCMS data as csv
# write.csv(OCG_LCMS_2021_subset, "data_csv/OCG_LCMS_2021_thresholded.csv")

# subset to match all data
OCG_GC_2021_ITS <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(asvITS_OCG)) # 48
asvITS_OCG_GC_2021 <- subset(asvITS_OCG, row.names(asvITS_OCG) %in% row.names(OCG_GC_2021_ITS)) # 48
mdITS_OCG_GC_2021 <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(asvITS_OCG_GC_2021)) # 48

# order the row names
OCG_GC_2021_ITS <- OCG_GC_2021_ITS[order(row.names(OCG_GC_2021_ITS)), ]
mdITS_OCG_GC_2021 <- mdITS_OCG_GC_2021[order(row.names(mdITS_OCG_GC_2021)), ]
asvITS_OCG_GC_2021 <- asvITS_OCG_GC_2021[order(row.names(asvITS_OCG_GC_2021)), ]

row.names(OCG_GC_2021_ITS) == row.names(asvITS_OCG_GC_2021) # sanity check: TRUE
row.names(mdITS_OCG_GC_2021) == row.names(asvITS_OCG_GC_2021) # sanity check: TRUE

# run PCA on GC 2021 data
OCG_GC_2021_ITS_no0 <- OCG_GC_2021_ITS[, colSums(OCG_GC_2021_ITS != 0) > 0]
scaled_OCG_GC_2021 <- scale(OCG_GC_2021_ITS_no0)

pca_gc_2021 <- prcomp(scaled_OCG_GC_2021, center = TRUE)
summary(pca_gc_2021) # 21.21 13.95

# save GC PCs in a dataframe
pca_scores_gc_2021 <- as.data.frame(pca_gc_2021$x[,])

fungal_dist_2021 <- vegdist(asvITS_OCG_GC_2021, method = "bray")

# Perform dbRDA
dbrda_2021_gc_model <- vegan::dbrda(fungal_dist_2021 ~ mdITS_OCG_GC_2021$x + mdITS_OCG_GC_2021$y + pca_scores_gc_2021$PC1 + pca_scores_gc_2021$PC2 + pca_scores_gc_2021$PC3, data = mdITS_OCG_GC_2021) 

# mdITS_OCG_2021_thermal_GC$AREA + mdITS_OCG_2021_thermal_GC$Variance_Plant_Temp + mdITS_OCG_2021_thermal_GC$Tleaf +

anova(dbrda_2021_gc_model, by = "term")

# dbrda_model_lcms_fun <- vegan::dbrda(fungal_dist_lcms ~ PC1 + PC2 + PC3 + PC4 + PC5 +  PC6 + PC7 + PC8 + PC9 + PC10 + x + y, data = combined_df_lcms) 
# anova(dbrda_model_lcms_fun, by = "term")
# adonis2_model_lcms <- adonis2(fungal_dist_lcms ~ PC1 + PC2 + PC3 + PC4 + PC5 +  PC6 + PC7 + PC8 + PC9 + PC10 + x + y, data = combined_df_lcms, by = "margin") 
# 
# 
# scores_pcoa <- scores(dbrda_model_2012, display = "sites")
# scores_pcoa <- as.data.frame(scores_pcoa)
# scores_pcoa$SampleID <- rownames(scores_pcoa)
# 
# arrow_scores <- scores(dbrda_model_2012, display = "bp")
# row.names(arrow_scores)
# rownames(arrow_scores) <- c("PC1", "PC7", "PC10")
# 
# # combined_df_2012 <- combined_df_2012[, -78]
# 
# plot(scores_pcoa[,1], scores_pcoa[,2], xlab = "PCoA1", ylab = "PCoA2",
#      main="Sagebrush 2012 fungal community by subspecies ploidy", 
#      col= c("pink","brown",'darkgreen','tan','lightblue')[combined_df_2012$Subsp_ploidy],
#      pch=17)
# scaling_factor <- 1.5
# arrows(0, 0, 
#        arrow_scores[,1] * scaling_factor, 
#        arrow_scores[,2] * scaling_factor, 
#        length = 0.1, col = "black")
# text(arrow_scores[,1]* scaling_factor, arrow_scores[,2]* scaling_factor, 
#      labels = rownames(arrow_scores), 
#      col = "black", pos = 3 , cex = 0.8)
# legend("topleft", 
#        legend=c("T_2n","T_4n","V_2n", "V_4n","W_4n"),
#        col= c("pink","brown",'darkgreen','tan','lightblue'),
#        pch=17,
#        cex=0.8,
#        bty = "n")
# 
# # Plot dbRDA using base R
# plot(scores_pcoa, xlab="PCoA1", ylab="PCoA2", 
#      main="Sagebrush 2012 fungal community by subspecies", 
#      col= c("pink","brown",'darkgreen','tan','lightblue')[combined_df_2012$Subsp_ploidy],
#      pch=17)
# arrow.plot(0, 0, col = "black", cex = 1.5, lwd = 2)
# legend("topleft", 
#        legend=c("T_2n","T_4n","V_2n", "V_4n","W_4n"),
#        col= c("pink","brown",'darkgreen','tan','lightblue'),
#        pch=17,
#        cex=0.8,
#        bty = "n")
# # ordispider(scores_pcoa,groups = combined_df_2012$Subsp_ploidy, show.groups = "T_2n", col = "olivedrab")
# # ordispider(scores_pcoa,groups = combined_df_2012$Subsp_ploidy, show.groups = "T_4n", col = "cadetblue")
# # ordispider(scores_pcoa,groups = combined_df_2012$Subsp_ploidy, show.groups = "V_2n", col = "goldenrod")
# # ordispider(scores_pcoa,groups = combined_df_2012$Subsp_ploidy, show.groups = "V_4n", col = "pink")
# # ordispider(scores_pcoa,groups = combined_df_2012$Subsp_ploidy, show.groups = "W_4n", col = "magenta")
# 
# # plot(dbrda_model_2012, type = "n")  # plot empty frame
# # points(dbrda_model_2012, display = "sites", col = "yellowgreen", pch = 17)  # sample points
# # text(dbrda_model_2012, display = "bp", col = "springgreen4")  # arrows and labels for predictors
# # title("dbRDA of 2012 fungal data")
# 

# # Mantel test for fungi and spatial dist ####
# # Extract coordinates matrix for populations
# site_coords <- md.OCG[, c("Longitude", "Latitude")]
# site_coords <- site_coords[order(row.names(site_coords)),]
# 
# # Subset and reorder both
# fungi_geo_samples <- intersect(rownames(asvITS_OCG.r), rownames(site_coords)) 
# fungi_geo <- asvITS_OCG.r[fungi_geo_samples, ]
# site_coords_fungi <- site_coords[fungi_geo_samples, ]
# 
# # create gc distance matrix
# fungi_dist <- vegdist(fungi_geo, method = "bray")
# 
# # create geographic distance matrix for GC samples
# geo_dist_fungi <- geodist(site_coords_fungi, measure = "geodesic")
# geo_dist_fungi <- as.dist(geo_dist_fungi)
# 
# # mantel test 
# mantel_geo_fungi <- mantel(geo_dist_fungi, fungi_dist, method = "spearman", permutations = 9999)
# print(mantel_geo_fungi) # not significant. Mantel statistic r: -0.08847; Significance: 0.9911
# 
# # Mantel test with common garden coordinate dist and fungi ####
# # Extract garden coordinates matrix 
# garden_coords <- md.OCG[, c("x", "y")]
# garden_coords <- garden_coords[order(row.names(garden_coords)),]
# 
# # Subset and reorder both
# fungi_garden_samples <- intersect(rownames(asvITS_OCG.r), rownames(garden_coords)) # 155 samples
# fungi_garden <- asvITS_OCG.r[fungi_garden_samples, ]
# garden_coords_fungi <- garden_coords[fungi_garden_samples, ]
# 
# # create gc distance matrix
# fungi_garden_dist <- vegdist(fungi_garden, method = "bray")
# 
# # create distance matrix in the garden
# garden_dist_fungi <- geodist(garden_coords_fungi, measure = "haversine")
# garden_dist_fungi <- as.dist(garden_dist_fungi)
# 
# # mantel test 
# mantel_garden_fungi <- mantel(garden_dist_fungi, fungi_garden_dist, method = "spearman", permutations = 9999)
# print(mantel_garden_fungi) # not significant. Mantel statistic r: -0.08847; Significance: 0.9911

# Mantel test for fungi and gc ####
# match row names for mantel test with fungal distance matrix and gc distance matrix
# GC
# read in cleaned and thresholded GC data 
OCG_GC <- read.csv("data_csv/OCG_GC_full_clean.csv", head = T, row.names = 1, check.names = F, stringsAsFactors = F)
OCG_GC[is.na(OCG_GC)] <- 0

gc_fun_samples <- intersect(rownames(asvITS_OCG.r), rownames(OCG_GC)) # 151 samples

# Subset and reorder both
gc_fun <- OCG_GC[gc_fun_samples, ]
fungal <- asvITS_OCG.r[gc_fun_samples, ]

row.names(gc_fun) == row.names(fungal) # sanity check: TRUE)

# create fungal distance matrix
fungi_dist <- vegdist(fungal, method = "bray")

# create gc data matrix
gc_fun_dist <- vegdist(gc_fun, method = "euclidean")

# mantel test 
mantel_gc_fungi <- mantel(fungi_dist, gc_fun_dist, method = "spearman", permutations = 999)
print(mantel_gc_fungi) 

# 2012 GC 
OCG_GC_2012 <- read.csv("data_csv/OCG_GC_2012_cleaned.csv", head = T, row.names = 1, check.names = F, stringsAsFactors = F)
OCG_GC_2012[is.na(OCG_GC_2012)] <- 0

gc12_fun_samples <- intersect(rownames(asvITS_OCG.r), rownames(OCG_GC_2012)) # 100 samples

# Subset and reorder both
gc_fun12 <- OCG_GC[gc12_fun_samples, ]
fungal12 <- asvITS_OCG.r[gc12_fun_samples, ]

row.names(gc_fun12) == row.names(fungal12) # sanity check: TRUE)

# create fungal distance matrix
fungi_dist12 <- vegdist(fungal12, method = "bray")

# create lcms data matrix
gc_fun_dist12 <- vegdist(gc_fun12, method = "euclidean")

# mantel test 
mantel_gc_fungi12 <- mantel(fungi_dist12, gc_fun_dist12, method = "spearman", permutations = 999)
print(mantel_gc_fungi12) # 0.0949, 0.023

# 2021 GC 
OCG_GC_2021 <- read.csv("data_csv/OCG_GC_2021.csv", head = T, row.names = 1, check.names = F, stringsAsFactors = F)
OCG_GC_2021[is.na(OCG_GC_2021)] <- 0

gc21_fun_samples <- intersect(rownames(asvITS_OCG.r), rownames(OCG_GC_2021)) # 48 samples

# Subset and reorder both
gc_fun21 <- OCG_GC_2021[gc21_fun_samples, ]
fungal21 <- asvITS_OCG.r[gc21_fun_samples, ]

row.names(gc_fun21) == row.names(fungal21) # sanity check: TRUE)

# create fungal distance matrix
fungi_dist21 <- vegdist(fungal21, method = "bray")

# create lcms data matrix
gc_fun_dist21 <- vegdist(gc_fun21, method = "euclidean")

# mantel test 
mantel_gc_fungi21 <- mantel(fungi_dist21, gc_fun_dist21, method = "spearman", permutations = 9999)
print(mantel_gc_fungi21) 

# Mantel test for fungi and lcms ####
# match row names for mantel test with fungal distance matrix and gc distance matrix
# LCMS
OCG_LCMS <- read.csv("data_csv/OCG_LCMS_thresholded.csv", head = T, row.names = 1, check.names = F, stringsAsFactors = F)
OCG_LCMS[is.na(OCG_LCMS)] <- 0

lcms_fun_samples <- intersect(rownames(asvITS_OCG.r), rownames(OCG_LCMS)) # 73 samples

# Subset and reorder both
lcms_fun <- OCG_LCMS[lcms_fun_samples, ]
fungal_lcms <- asvITS_OCG.r[lcms_fun_samples, ]
row.names(lcms_fun) == row.names(fungal_lcms) # sanity check: TRUE)

# create fungal distance matrix
fungi_lcms_dist <- vegdist(fungal_lcms, method = "bray")

# create lcms data matrix
lcms_fun_dist <- vegdist(lcms_fun, method = "euclidean")

# mantel test 
mantel_lcms_fungi <- mantel(fungi_lcms_dist, lcms_fun_dist, method = "spearman", permutations = 999)
print(mantel_lcms_fungi) # not significant. Mantel statistic r: 0.03324; Significance: 0.2863

# 2012 LCMS vs fungi
lcms_fun12_samples <- intersect(rownames(asvITS.2012.r), rownames(OCG_LCMS)) # 24 samples

# Subset and reorder both
lcms_fun12 <- OCG_LCMS[lcms_fun12_samples, ]
fungal_lcms12 <- asvITS.2012.r[lcms_fun12_samples, ]
row.names(lcms_fun12) == row.names(fungal_lcms12) # sanity check: TRUE)

# create fungal distance matrix
fungi_lcms12_dist <- vegdist(fungal_lcms12, method = "bray")

# create lcms data matrix
lcms_fun12_dist <- vegdist(lcms_fun12, method = "euclidean")

# mantel test 
mantel_lcms_fungi12 <- mantel(fungi_lcms12_dist, lcms_fun12_dist, method = "spearman", permutations = 9999)
print(mantel_lcms_fungi12) #-0.1513, 0.9127

# 2021 LCMS vs fungi
lcms_fun21_samples <- intersect(rownames(asvITS.2021.r), rownames(OCG_LCMS)) # 24 samples

# Subset and reorder both
lcms_fun21 <- OCG_LCMS[lcms_fun21_samples, ]
fungal_lcms21 <- asvITS.2021.r[lcms_fun21_samples, ]
row.names(lcms_fun21) == row.names(fungal_lcms21) # sanity check: TRUE)

# create fungal distance matrix
fungi_lcms21_dist <- vegdist(fungal_lcms21, method = "bray")

# create lcms data matrix
lcms_fun21_dist <- vegdist(lcms_fun21, method = "euclidean")

# mantel test 
mantel_lcms_fungi21 <- mantel(fungi_lcms21_dist, lcms_fun21_dist, method = "spearman", permutations = 999)
print(mantel_lcms_fungi21) #-0.04258, 0.6984

# Mantel test with thermal data and 2021 fungal data
thermal_data <- read.csv("data_csv/orchard_2019_thermal_data.csv", header = T, check.names = F)

#subset to just June season
thermal_data <- thermal_data[thermal_data$Season == "June", ]

colnames(thermal_data)[colnames(thermal_data) == "ID"] <- "Garden Plant ID"

#remove rows with NAs from thermal data
thermal_data <- thermal_data[!is.na(thermal_data$UAS_Thermal), ]

# subset thermal need to be matched with metadata.
mdITS_OCG_2021 <- mdITS_OCG[mdITS_OCG$Year == "2021", ]
mdITS_OCG_2021_thermal <- merge(mdITS_OCG_2021, thermal_data, by = "Garden Plant ID")
# mdITS_OCG_2021_thermal <- mdITS_OCG_2021_thermal[!duplicated(mdITS_OCG_2021_thermal$`Garden Plant ID`), ]

#make description row names
mdITS_OCG_2021_thermal$Description -> rownames(mdITS_OCG_2021_thermal)

therm_fun21_samples <- intersect(rownames(asvITS.2021.r), rownames(mdITS_OCG_2021_thermal)) #29 samples

# subset thermal data to columns 31-39
mdITS_OCG_2021_thermal <- mdITS_OCG_2021_thermal[, c(31:39)]


# Subset and reorder both
therm_fun21 <- mdITS_OCG_2021_thermal[therm_fun21_samples, ]
fungal_therm21 <- asvITS.2021.r[therm_fun21_samples, ]
row.names(therm_fun21) == row.names(fungal_therm21) # sanity check: TRUE)

# create fungal distance matrix
fungal_therm21_dist <- vegdist(fungal_therm21, method = "bray")

# create lcms data matrix
therm_fun21_dist <- vegdist(therm_fun21, method = "euclidean")

# mantel test 
mantel_therm_fungi21 <- mantel(fungal_therm21_dist, therm_fun21_dist, method = "spearman", permutations = 9999)
print(mantel_therm_fungi21) #-0.1413 0.8649


#stable isotope analysis####
stable_iso_data_2012<-read.csv("data_csv/Stable_isotope_2012.csv")

stable_iso_data_2021<-read.csv("data_csv/Stable_isotope_2021.csv")
str(stable_iso_data_2021)
stable_iso_data_2021$Subspecies<-as.factor(stable_iso_data_2021$Subspecies)
levels(stable_iso_data_2021$Subspecies)

boxplot(stable_iso_data_2012$Delta15N~stable_iso_data_2012$Subspecies)

boxplot(stable_iso_data_2012$Delta13C~stable_iso_data_2012$Subspecies)

boxplot(stable_iso_data_2021$d15N~stable_iso_data_2021$Subspecies) #extra level in here for "T/W"
boxplot(stable_iso_data_2021$d13C~stable_iso_data_2021$Subspecies)

## Remove V-134677 (not sure what subspecies since it says "T/W"). look into later
stable_iso_data_2021 <- stable_iso_data_2021[!(row.names(stable_iso_data_2021) == "31"),]

stable_12_N15_dist <- vegdist(stable_iso_data_2012$Delta15N) #Bray-curtis distance between samples: quantifies differences in the overall taxonomic composition between two samples

groups_12 <- factor(stable_iso_data_2012$Subspecies, labels = c("T","V","W")) #using the md for year 

## Calculate multivariate dispersions
mod_12 <- betadisper(stable_12_N15_dist, groups_12)
mod_12

#Perform test
anova(mod_12)

#Permutation test for F
permutest(mod_12, permutations = 999,pairwise = TRUE) 

# stable_12_C13_dist <- vegdist(stable_iso_data_2012$Delta13C) #Bray-curtis distance between samples: quantifies differences in the overall taxonomic composition between two samples
# 
# ## Calculate multivariate dispersions
# mod_12_C13 <- betadisper(stable_12_C13_dist, groups_12)
# mod_12_C13
# 
# #Perform test
# anova(mod_12_C13)
# 
# TukeyHSD(p1)

pw.comparison.15N_12 <- aov(Delta15N ~ Subspecies, data = stable_iso_data_2012)
summary(pw.comparison.15N_12)
TukeyHSD(pw.comparison.15N)

pw.comparison.13C_12 <- aov(Delta13C ~ Subspecies, data = stable_iso_data_2012)
summary(pw.comparison.13C_12)
TukeyHSD(pw.comparison.13C)

# 2021
pw.comparison.15N_21 <- aov(d15N ~ Subspecies, data = stable_iso_data_2021)
summary(pw.comparison.15N_21)
TukeyHSD(pw.comparison.15N_21)

pw.comparison.13C_21 <- aov(d13C ~ Subspecies, data = stable_iso_data_2021)
summary(pw.comparison.13C_21)
TukeyHSD(pw.comparison.13C_21)

# #pairwiseadonis
# stable_iso.subsp.pw <- pairwise.adonis(stable_iso_data_2012$Delta15N, stable_iso_data_2012$Subspecies)
# # # stable_iso.subsp.pw #T vs V= < 2e-16, T vs W= 5.25e-05, and W vs V= 0.0391.
# 
# pairwise.adonis<-pairwise.adonis2(stable_iso_data_2012$Delta15N ~ stable_iso_data_2012$Subspecies)
# pairwise.adonis

ggplot(stable_iso_data_2012,aes(Subspecies,Delta15N))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+ggtitle("2012")+ theme(legend.position = "none")

ggplot(stable_iso_data_2021,aes(Subspecies,d15N))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+
  ylab("Delta15N")+ggtitle("2021")+theme(legend.position = "none")

library(gridExtra)
grid.arrange(p1, p2, nrow = 1)

ggplot(stable_iso_data_2012,aes(Delta13C,Delta15N))+
  geom_point(aes(color=Subspecies))+theme_classic()+
  scale_color_manual(values=c("olivedrab","cadetblue","goldenrod"))+
  ggtitle("2012")+theme(legend.position = "none")+geom_text(aes(label=ID))

ggplot(stable_iso_data_2021,aes(d13C,d15N))+
  geom_point(aes(color=Subspecies))+theme_classic()+
  scale_color_manual(values=c("olivedrab","cadetblue","goldenrod"))+ ggtitle("2021")+theme(legend.position = "none")+
  geom_text(aes(label=Lab))

library(ggpubr)
ggscatterhist(
  stable_iso_data_2012, x = "Delta13C", y = "Delta15N", group="Subspecies",
  color = "Subspecies", fill= "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.plot = "boxplot",
  margin.params = list(fill = "Subspecies", color = c("olivedrab", "cadetblue", "goldenrod"), size = 0.2),
  ggtheme = theme_bw()
)+
  geom_text(aes(label=))

ggscatterhist(
  stable_iso_data_2021, x = "d13C", y = "d15N", group="Subspecies",
  color = "Subspecies", fill= "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.plot = "boxplot",
  margin.params = list(fill = "Subspecies", color = c("olivedrab", "cadetblue", "goldenrod"), size = 0.2),
  ggtheme = theme_bw()
)+
  geom_text(aes(label=Sample.ID))

# Grouped Scatter plot with marginal density plots
ggscatterhist(
  stable_iso_data_2012, x = "Delta13C", y = "Delta15N", group = "Subspecies",
  color = "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.params = list(fill = "Subspecies", color = "black", size = 0.2)
)

# combine the 2012 stable isotope data and 2021 stable isotope data using "Sample.ID" column

# subset the 2012 isotope data to just have Delta15N and Delta13C and Lab columns
stable_iso_data_2012_subset <- stable_iso_data_2012[, c("Sample.ID", "Delta15N", "Delta13C", "Subspecies", "Year")]

stable_iso_data_2021_subset <- stable_iso_data_2021[, c("Sample.ID", "d15N", "d13C", "Subspecies", "Year")]

# # rename the columns in the 2021 isotope data to match the 2012 isotope data
# colnames(stable_iso_data_2021_subset) <- c("Sample.ID", "Delta15N", "Delta13C", "Subspecies")

length(intersect(stable_iso_data_2012_subset$Sample.ID, stable_iso_data_2021_subset$Sample.ID))

merged_isotope_df <- merge(stable_iso_data_2012_subset, stable_iso_data_2021_subset, by = "Sample.ID", all = TRUE) # 141

# write csv 
write.csv(merged_isotope_df, "data_csv/merged_stable_isotope_data.csv", row.names = FALSE)


# match sample id from isotope 2012 to asv 2012
iso_2012_asv_samples <- intersect(rownames(asvITS.2012), stable_iso_data_2012_subset$Sample.ID) # 61 samples

#####














#Trash ####
# gc_dist <- vegdist(OCG_GC, method = "euclidean")

# pca_scores_GC.f <- data.pca_GC_ID$x[,]

# pca_GC_f <- as.data.frame(pca_scores_GC.f)
# pca_GC_f$SampleID <- rownames(pca_GC_f)

# pca_21GC_f <- data.pca_2021_ID$x[,1]
# pca_21_GC_f <- as.data.frame(pca_21GC_f)
# pca_21_GC_f$SampleID <- rownames(pca_21_GC_f)
# 
# pca_GC_f <- data.pca_2021_ID$x[,1]
# pca_21_GC_f <- as.data.frame(pca_21GC_f)
# pca_21_GC_f$SampleID <- rownames(pca_21_GC_f)

# pca_scores_LCMS.f <- pca_LCMS$x[,1]
# pca_LCMS_f <- as.data.frame(pca_scores_LCMS.f)
# pca_LCMS_f$SampleID <- rownames(pca_LCMS_f)


# # Read in thermal data and add it to md
# thermal_data <- read.csv("data_csv/orchard_2019_thermal_data.csv", header = T, check.names = F)
# thermal_data <- thermal_data[thermal_data, ]
# colnames(thermal_data)[colnames(thermal_data) == "ID"] <- "Garden Plant ID"
# 
# # subset to just the season == 'June'
# thermal_data <- thermal_data[thermal_data$Season == "June", ]

# md_OCG_subset <- mdITS_OCG[, c("Garden Plant ID", "x", "y")]
# mdITS_OCG <- merge(mdITS_OCG, md_OCG_subset, by = "Garden Plant ID", all.x = TRUE)
# mdITS_OCG <- mdITS_OCG[!duplicated(mdITS_OCG$Description), ]
# mdITS_OCG_thermal <- merge(mdITS_OCG, thermal_data, by = "Garden Plant ID", all.x = TRUE)
# 
# #make description row names 
# mdITS_OCG_thermal$SampleID -> rownames(mdITS_OCG_thermal)
# mdITS_OCG.df <- mdITS_OCG_thermal[, c("Garden Plant ID", "Subsp_ploidy", "Year", "UAS_Thermal", "x","y")]
# 
# # Merge all on 'SampleID'
# rownames(combined_df) <- combined_df$SampleID
# combined_df$SampleID <- NULL
# combined_df <- na.omit(combined_df)
# 
# # Get sample names after NA removal
# keep_samples <- rownames(mdITS_OCG_thermal)
# 
# # Subset the distance matrix
# fungal_dist_clean <- as.dist(as.matrix(fungal_dist)[keep_samples, keep_samples])

# # Create spatial distance matrix
# coords <- mdITS_OCG[, c("x", "y")]
# spatial_dist <- dist(coords)

# # Build PCNM (Principal Coordinates of Neighbor Matrices)
# pcnm_model <- pcnm(spatial_dist)

# so coords, metadata, and fungal data all match
# mdITS_OCG <- mdITS_OCG[mdITS_OCG$Garden Plant ID %in% rownames(asvITS_OCG.r), ] # 166 samples

#Read in project data for ANCOM
# #ASV
# asvITS <- read.csv("data_csv/asv-table-dada2-ITS-sagebrush.csv",head=T,row.names=1, check.names = F) #5983 obs of 463 var
# asvITS <- asvITS[,order(colnames(asvITS))] # order samples alphabetically
# 
# #METADATA
# mdITS <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) # 505 obs of 21 var
# mdITS <- mdITS[order(row.names(mdITS)),] 
# mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS)) 
# colnames(asvITS) == row.names(mdITS) # sanity check:TRUE
# 
# #TAXONOMY
# tax.ITS <- read.csv("~/Documents/Orchard_Common_Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #5983 obs of 2 variables
# row.names(asvITS) == row.names(tax.ITS) #TRUE

#Cleaning
##Remove entries with insufficient sequences
# asvITS_OCG[asvITS < 10] <- 0
# asvITS <- asvITS[rowSums(asvITS) > 0,] # each observation needs at least 10 seqs
# 
# summary(rowSums(asvITS)) #10
# summary(colSums(asvITS)) #0
# sort(colSums(asvITS)) 
# 
# 
# asvITS <- asvITS[,colSums(asvITS) > 499] # each sample needs at least 1000 seqs 
# 
# mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS)) # remove trimmed samples from metadata #361 obs o f 16 var
# 
# mdITS_OCG <- subset(mdITS, mdITS$Project == "OCG") #146 obs of 16 var   
# mdITS_OCG <- subset(mdITS_OCG, mdITS_OCG$Location != "NEG")  
# asvITS_OCG.s <- subset(asvITS_OCG, row.names(asvITS_OCG) %in% row.names(mdITS_OCG)) 
# asvITS_OCG <- as.data.frame(asvITS_OCG.s) 
# 
# asvITS_OCG[asvITS_OCG < 10] <- 0  # repeat cleaning after trimming
# asvITS_OCG <- asvITS_OCG[rowSums(asvITS_OCG) > 0,] # each observation needs at least 10 seqs
# 
# summary(rowSums(asvITS_OCG)) #16
# summary(colSums(asvITS_OCG)) #1002
# 
# asvITS_OCG <- asvITS_OCG[,colSums(asvITS_OCG) > 999] # each sample needs at least 1000 seqs

#READ in project data for metacoder 
#Read in project data 
#ASV
# asvITS<- read.csv("data_csv/asv-table-dada2-ITS-sagebrush.csv",head=T,row.names=1, check.names = F) #5983 obs of 463 variable
# asvITS<- asvITS[,order(colnames(asvITS))]
# 
# ##METADATA
# mdITS <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #505 obs of 16 variables.
# mdITS <- mdITS[order(row.names(mdITS)),]
# mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS)) #361 of 16 variables
# colnames(asvITS) == row.names(mdITS) # sanity check true.
# 
# ##TAXONOMY
# tax.ITS <- read.csv("~/Documents/Orchard_Common_Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #5983 obs of 2 variables
# row.names(asvITS) == row.names(tax.ITS) #TRUE
# 
# #Clean up data for analysis
# asvITS[asvITS < 10] <- 0
# asvITS <- asvITS[rowSums(asvITS) > 0,]
# summary(rowSums(asvITS)) #10
# summary(colSums(asvITS)) #0
# asvITS <- asvITS[,colSums(asvITS) > 999]
# mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS)) # remove trimmed samples from metadata
# 
# ##Subset data frame to only include samples from the Orchard Common Garden project #
# mdITS_OCG <- subset(mdITS, mdITS$Project == "OCG")  # trim metadata to only OCG entries
# mdITS_OCG <- subset(mdITS_OCG, mdITS_OCG$Location != "NEG")# remove negative controls
# print(mdITS_OCG)
# asvITS.t <- t(asvITS) # transpose asv data frame before subsetting from metadata (subsetting didn't play well when comparing rows to columns or vice versa)
# asvITS_OCG.t <- subset(asvITS.t, row.names(asvITS.t) %in% row.names(mdITS_OCG))  # trim transposed asv data frame based on trimmed metadata
# asvITS_OCG <- t(asvITS_OCG.t)  # transpose asv table back
# asvITS_OCG <- as.data.frame(asvITS_OCG)# return asv table to data frame format after transposing
# 
# ##Trimming of improper samples from OCG subset #
# asvITS_OCG[asvITS_OCG < 10] <- 0                  
# asvITS_OCG <- asvITS_OCG[rowSums(asvITS_OCG) > 0,] # each observation needs at least 10 seqs
# 
# summary(rowSums(asvITS_OCG)) #10
# summary(colSums(asvITS_OCG)) #1018
# 
# asvITS_OCG <- asvITS_OCG[,colSums(asvITS_OCG) > 999] # each sample needs at least 1000 seqs
# 
# mdITS_OCG <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% colnames(asvITS_OCG)) # remove trimmed samples from metadata

# # average effective species number of neighboring plants for each plant in dataset
# 
# # set a distance threshold to define neighbors (same units are garden coordinates)
# neighbor_radius <- 1
# 
# #create a distance matrix between all plants
# coords <- md.OCG.its[, c("x", "y")]
# distance_matrix <- as.matrix(dist(coords))
# 
# # create an empty vector to store average effective species number of neighbors
# avg_effective_species_neighbors <- numeric(nrow(md.OCG.its))
# 
# # Loop through each plant
# for (i in 1:nrow(distance_matrix)) {
#   # Identfiy neighbors within the throwhold (excluding the plant itself)
#   neighbors <- which(distance_matrix[i, ] > 0 & distance_matrix[i, ] <= neighbor_radius)
#   
#   # if neighbors are found, compute their mean effective species number
#   if (length(neighbors) > 0) {
#     avg_effective_species_neighbors[i] <- mean(md.OCG.its$effective_species[neighbors], na.rm = TRUE)
#   } else {
#     avg_effective_species_neighbors[i] <- 0 # or NA depending on your preference
#   }
# }
# 
# md.OCG.its$avg_eff_sp_neighbors <- avg_effective_species_neighbors
# 
# ggplot(md.OCG.its, aes(x = x, y = y)) +
#   geom_point(aes(color = avg_eff_sp_neighbors, size = effective_species)) +
#   scale_color_viridis_c(option = "plasma", na.value = "grey80") +
#   scale_size(range = c(2, 6)) +
#   theme_minimal() +
#   labs(
#     title = "Spatial Map of Average Effective Species Number (Neighbors)",
#     x = "Garden X Position",
#     y = "Garden Y Position",
#     color = "Avg Effective\nSpecies (Neighbors)",
#     size = "Effective\nSpecies (Plant)"
#   ) + coord_fixed()  # Ensures equal x/y scaling
# 
# # create an empty vector to store average compound number of neighbors
# avg_compounds_neighbors <- numeric(nrow(md.OCG.its))
# 
# # Loop through each plant
# for (i in 1:nrow(distance_matrix)) {
#   # Identfiy neighbors within the throwhold (excluding the plant itself)
#   neighbors <- which(distance_matrix[i, ] > 0 & distance_matrix[i, ] <= neighbor_radius)
#   
#   # if neighbors are found, compute their mean effective species number
#   if (length(neighbors) > 0) {
#     avg_compounds_neighbors[i] <- mean(md.OCG.its$compounds[neighbors], na.rm = TRUE)
#   } else {
#     avg_compounds_neighbors[i] <- 0 # or NA depending on your preference
#   }
# }
# 
# md.OCG.its$avg_comp_neighbors <- avg_compounds_neighbors
# 
# # linear model 
# 
# eco_mod <- lm(effective_species ~ Ecoregion + avg_eff_sp_neighbors + avg_comp_neighbors, data = md.OCG.its)
# summary(eco_mod)
# 
# mod <- lm(effective_species ~ Subspecies + Ploidy + Year + avg_eff_sp_neighbors + avg_comp_neighbors, data = md.OCG.its)
# summary(mod)
# 
# modglm <- glm(effective_species ~ Subspecies + Ploidy  + Year + avg_eff_sp_neighbors + avg_comp_neighbors, family = poisson, data = md.OCG.its)
# summary(modglm)
# 
# hist(md.OCG.its$effective_species)
# 
# modme <- lme4::glmer(effective_species ~ Subspecies + Ploidy + Year + Ecoregion + avg_eff_sp_neighbors + avg_comp_neighbors + (1|Plant), family = poisson, data = md.OCG.its)
# 
# summary(modme)



# k- means clustering 
#read in lcms with cluster data saved
md_LCMS_cluster <- read.csv("data_csv/md.OCG.LCMS_cluster.csv", row.names = 1)
md_LCMS_cluster <- md_LCMS_cluster[order(row.names(md_LCMS_cluster)),]
md_LCMS_cluster <- subset(md_LCMS_cluster, row.names(md_LCMS_cluster) %in% row.names(asvITS_OCG))
asvITS_OCG_lcms <- subset(asvITS_OCG, row.names(asvITS_OCG) %in% row.names(md_LCMS_cluster))
row.names(asvITS_OCG_lcms) == row.names(md_LCMS_cluster) # sanity check:TRUE

# rarefy the asv data
asvITS_OCG_lcms <- asvITS_OCG_lcms[,colSums(asvITS_OCG_lcms) > 999]
summary(colSums(asvITS_OCG_lcms)) #1002
summary(rowSums(asvITS_OCG_lcms)) #17
asvITS_OCG_lcms.r <- rrarefy(asvITS_OCG_lcms, 65)

# model
k_means_asv_fit <- adonis2(asvITS_OCG_lcms.r ~ md_LCMS_cluster$cluster_assignments, by = "margin")
k_means_asv_fit #subspecies ploidy is significant 0.001

# read in GC clustering data
md_gc_2012_cluster <- read.csv("data_csv/md.OCG.GC.2012_cluster.csv")
md_gc_2021_cluster <- read.csv("data_csv/md.OCG.GC.2021_cluster.csv")
md_gc_cluster <- merge(md_gc_2012_cluster, md_gc_2021_cluster, all = TRUE)
rownames(md_gc_cluster) <- md_gc_cluster[, 1]
md_gc_cluster <- md_gc_cluster[, -1]

md_gc_cluster <- md_gc_cluster[order(row.names(md_gc_cluster)),]
md_gc_cluster <- subset(md_gc_cluster, row.names(md_gc_cluster) %in% row.names(asvITS_OCG)) ##166 of 21 var
asvITS_OCG_gc <- subset(asvITS_OCG, row.names(asvITS_OCG) %in% row.names(md_gc_cluster)) ##166 of 21 var
row.names(asvITS_OCG_gc) == row.names(md_gc_cluster) # sanity check:TRUE

# rarefy the asv data
asvITS_OCG_gc <- asvITS_OCG_gc[,colSums(asvITS_OCG_gc) > 999]
summary(colSums(asvITS_OCG_gc)) #1002
summary(rowSums(asvITS_OCG_gc)) #16
asvITS_OCG_gc.r <- rrarefy(asvITS_OCG_gc, 16)

# model
k_means_asv.gc_fit <- adonis2(asvITS_OCG_gc.r ~ md_gc_cluster$cluster_assignments, by = "margin")
k_means_asv.gc_fit #subspecies ploidy is significant 0.001

asv_subsp_GC <- adonis2(asvITS_OCG_gc.r ~ md_gc_cluster$Subspecies, by = "margin")

asvnmds_gc_cluster <- metaMDS(asvITS_OCG_gc.r, trymax = 500)
ordiplot(asvnmds_gc_cluster, type = "t", display = "sites",cex = .6)

plot(asvnmds_gc_cluster$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2012 fungal community by subspecies", 
     col= c("olivedrab","cadetblue","goldenrod", "magenta", "green")[md_gc_cluster$cluster_assignments],
     pch=19)
legend("topright", 
       legend=c("T","V","W"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subspecies, show.groups = "W", col = "goldenrod")

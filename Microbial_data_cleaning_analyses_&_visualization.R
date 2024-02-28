# Orchard Common Garden - Microbial data cleaning, analysis, and visualization####
# Set working directory and load necessary packages####
setwd("/Users/ellehorwath/Documents/Orchard_Common _Garden/commongarden")
if (!require("readr")) {install.packages("readr"); require("readr")}
if (!require("dplyr")) {install.packages("dplyr"); require("dplyr")}
if (!require("tidyr")) {install.packages("tidyr"); require("tidyr")}
if (!require("ggplot2")) {install.packages("ggplot2"); require("ggplot2")}
if (!require("MASS")) {install.packages("MASS"); require("MASS")}
if (!require("vegan")) {install.packages("vegan"); require("vegan")}
if (!require("effects")) {install.packages("effects"); require("effects")}
if (!require("pairwiseAdonis")) {devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"); require("pairwiseAdonis")}

#Read in data ####
## Read ASV data in ####
#Amplicon Sequence Variant table read in.
#Data has not been filtered and is not yet clean to include just observations with at least 10 seqs and each sample needs at least 1000 seq. The data has not been transposed and is ordered alphabetically. 

asvITS<- read.csv("data_csv/asv-table-dada2-ITS-sagebrush.csv",head=T,row.names=1, check.names = F) #read in Amplicon Sequence Variant table. 5983 obs of 463 variable
asvITS<- asvITS[,order(colnames(asvITS))] # order samples alphabetically
summary(rowSums(asvITS)) #1.0
summary(colSums(asvITS)) #0

##Full metadata read in####
#Data includes metadata from numerous projects. I will be subsetting for just the orchard common garden plants. The data will include plant ID number, location of origin, year of sampling (2012 or 2021), subspecies, ploidy, and subspecies ploidy, status (dead or alive) in 2020.

mdITS <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #set to correct file path 505 obs of 16 variables.
mdITS <- mdITS[order(row.names(mdITS)),] # order samples alphabetically
mdITS <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS)) #463 of 16 variables
colnames(asvITS) == row.names(mdITS) # sanity check: a check to make sure something does not contain elementary mistakes or impossibilities and is not based on invalid assumptions
#This sanity check reads true.

#change plant description ID to end in 2021 not 2020
# Replace "2020" with "2021" at the end of the strings
mdITS$Description <- sub("_2020$", "_2021", mdITS$Description)

##Taxonomy table read in####
#taxonomy table is used to match to amplicon sequence variant table to fungal ID.
tax.ITS <- read.csv("~/Documents/Orchard_Common _Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #5983 obs of 2 variables
row.names(asvITS) == row.names(tax.ITS) # sanity check: TRUE


#Cleaning data####
##ASV table cleaning ####
#asvITS[asvITS < 10] <- 0 # each observation needs at least 10 seqs.
asvITS <- asvITS[rowSums(asvITS) > 0,] #showing me the values that are greater than zero
summary(rowSums(asvITS)) #1
summary(colSums(asvITS)) #0

asvITS <- asvITS[,colSums(asvITS) > 499] # each sample needs at least 500 seqs. how low can I go for the ordiarrows to have them all have a pair. #5983 of 385 var

#asvITS <- asvITS[,colSums(asvITS) > 999] # each sample needs at least 1000 seqs. how low can I go for the ordiarrows to have them all have a pair. #5983 of 361 var

summary(colSums(asvITS)) #507
summary(rowSums(asvITS)) #0

asvITS.t <- t(asvITS) # transpose rows and columns
asvITS.t <- asvITS.t[order(row.names(asvITS.t)),] # order samples alphabetically
asvITS.t <- asvITS.t[,order(colnames(asvITS.t))] # order asvs alphabetically

summary(rowSums(asvITS.t)) #making sure that it transposed
summary(colSums(asvITS.t))

mdITS <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.t)) #subsetting metadata to match what is left in the asv table

asvITS.t2 <- asvITS.t[!(row.names(asvITS.t) %in% c("NEG_9-30-21","AH1919","AHM20207","AHM20125","UTW.1.4_2021")),] #outliers that need to be removed

asvITS.t2 <- asvITS.t2[,colSums(asvITS.t2) > 0] #keeping samples greater than 0

summary(rowSums(asvITS.t2)) #507
summary(colSums(asvITS.t2)) #2.0

mdITS2 <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.t2)) # again subsetting metadata to match what is in the asv table #380 obs of 16 var

### Subsetting to just the plant in the common garden (OCG)####
asvITS.OCG <- subset(asvITS.t2, mdITS2$Project=="OCG") #subsetting to just the Orchard common garden plants
asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 0]

summary(rowSums(asvITS.OCG)) #507
summary(colSums(asvITS.OCG)) #2.0

### Remove duplicates from asv####
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1')
asvITS.OCG <- asvITS.OCG[!rownames(asvITS.OCG) %in% rows_to_remove, ]

# asvITS.OCG <- asvITS.OCG[grep("v2", row.names(asvITS.OCG), invert = T),]

## Remove negative control
asvITS.OCG <- asvITS.OCG[!(row.names(asvITS.OCG) == "NEG_8-28-21"),]

## Remove MTW.3.7.R_2012
asvITS.OCG <- asvITS.OCG[!(row.names(asvITS.OCG) == "MTW.3.7.R_2012"),] #we arent sure what the R represents. 

## subset metadata to match asv table samples
asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 0]
summary(rowSums(asvITS.OCG)) #507
summary(colSums(asvITS.OCG)) #2.0

mdITS.OCG <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.OCG)) #subset md to match asv table samples #154 of 16 var

#Write csv for cleaned metadata and asv table####
write.csv(asvITS.OCG, file = "data_csv/asvITS.OCG.csv") #write csv for asv table of OCG only
write.csv(mdITS.OCG, file = "data_csv/metadata_OCG.csv") #write csv for metadata of OCG only

#Clear Global Environment #############################

#Read in cleaned data ####
mdITS_OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #metadata read in #154 of 16 variables
head(mdITS_OCG)

asvITS.OCG <- read.csv("data_csv/asvITS.OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #asv table read in 154 obs of 2135 variables
head(asvITS.OCG)
summary(rowSums(asvITS.OCG)) #507
summary(colSums(asvITS.OCG)) #2

row.names(asvITS.OCG) == row.names(mdITS_OCG) # sanity check:TRUE

tax.ITS <- read.csv("~/Documents/Orchard_Common _Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #taxonomy read in

##Rarefying####

asvITS.OCG.r <- rrarefy(asvITS.OCG,507) ## rarefy: Warning message

mdITS_OCG <- droplevels(mdITS_OCG) #drop levels for NMDS

set.seed(41)
#asvITS.OCG.nmds <- metaMDS(asvITS.OCG.r, trymax=500) ### Solution reached! gives me a warning
# save(asvITS.OCG.nmds, file = "nmds/asvITS_OCG_nmds.rda") #save the nmds so you won't need to run it again
load("nmds/asvITS_OCG_nmds.rda") #load it to use in code anytime after the initial run

#NMDS ####
ordiplot(asvITS.OCG.nmds, type = "t",display = "sites",cex = .6)

#make variables factor to plot
mdITS_OCG$Ploidy <- as.factor(mdITS_OCG$Ploidy)
mdITS_OCG$Subspecies <- as.factor(mdITS_OCG$Subspecies)
mdITS_OCG$Subsp_ploidy <- as.factor(mdITS_OCG$Subsp_ploidy)
mdITS_OCG$Year <- as.factor(mdITS_OCG$Year)
mdITS_OCG$Plant <- as.factor(mdITS_OCG$Plant)

## NMDS for fungal community by subspecies ####
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="Sagebrush fungal community by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS_OCG$Subspecies],
     pch=c(19,17,19)[mdITS_OCG$Year])
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.OCG.nmds,groups = mdITS_OCG$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(asvITS.OCG.nmds,groups = mdITS_OCG$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(asvITS.OCG.nmds,groups = mdITS_OCG$Subspecies, show.groups = "W", col = "goldenrod")

### PERMANOVA and adonis for subspecies ####
asvITS.OCG.subsp <- adonis2(asvITS.OCG.r ~ mdITS_OCG$Subspecies,by="margin") # Bray-Curtis is the default metric
asvITS.OCG.subsp #subspecies significant

#pairwiseadonis
asvITS.OCG.subsp.pw <- pairwise.adonis(asvITS.OCG.r, mdITS_OCG$Subspecies)
asvITS.OCG.subsp.pw #T vs V= 0.003, T vs W= 0.084, and W vs V= 0.540.

## NMDS for fungal community by ploidy ####
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by ploidy", 
     col= c("red","blue")[mdITS_OCG$Ploidy],
     pch=c(19,17)[mdITS_OCG$Year])
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.OCG.nmds,groups = mdITS_OCG$Ploidy, show.groups = "2n", col = "red")
ordispider(asvITS.OCG.nmds,groups = mdITS_OCG$Ploidy, show.groups = "4n", col = "blue")

### PERMANOVA for ploidy####
asvITS.OCG.ploidy <- adonis2(asvITS.OCG.r ~ mdITS_OCG$Ploidy,by="margin") # Bray-Curtis is the default metric
asvITS.OCG.ploidy #ploidy is not significant

## NMDS for fungal community by subsp ploidy ####
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS_OCG$Subsp_ploidy],
     pch=c(19,17)[mdITS_OCG$Year])
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.OCG.nmds,groups = mdITS_OCG$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(asvITS.OCG.nmds,groups = mdITS_OCG$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(asvITS.OCG.nmds,groups = mdITS_OCG$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(asvITS.OCG.nmds,groups = mdITS_OCG$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(asvITS.OCG.nmds,groups = mdITS_OCG$Subsp_ploidy, show.groups = "W_4n", col = "purple")

### PERMANOVA and adonis for subspecies ploidy ####
asvITS.OCG.subsp_ploi <- adonis2(asvITS.OCG.r ~ mdITS_OCG$Subsp_ploidy) 
asvITS.OCG.subsp_ploi #subspecies ploidy is significant 0.002

#pairwiseadonis
asvITS.OCG.subsp.pw <- pairwise.adonis(asvITS.OCG.r, mdITS_OCG$Subsp_ploidy)
asvITS.OCG.subsp.pw #p-adjusted:T_4n vs V_2n = 0.04, T_4n vs V_4n = 0.03, T_4n vs W_4n = 0.05

## NMDS for fungal community by year ####
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by year", 
     col= rainbow(2)[mdITS_OCG$Year],
     pch=19)
legend("topleft", 
       legend=c("2012","2021"),
       col= c("red","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")
ordiarrows(asvITS.OCG.nmds, mdITS_OCG$Plant)

levels(mdITS_OCG$Plant)

### PERMANOVAs for year ####
asvITS.OCG.subsp_loc_yr <- adonis2(asvITS.OCG.r ~ mdITS_OCG$Subspecies + mdITS_OCG$Year + mdITS_OCG$Location, by = "margin") 
asvITS.OCG.subsp_loc_yr #Year is significant

asvITS.OCG.subsp_yr <- adonis2(asvITS.OCG.r ~ mdITS_OCG$Subspecies*mdITS_OCG$Year) 
asvITS.OCG.subsp_yr #year and subspecies are sig

asvITS.OCG.yr <- adonis2(asvITS.OCG.r ~ mdITS_OCG$Year) 
asvITS.OCG.yr #year is significant


## 2012 NMDS fungal community by subspecies####
asvITS.2012 <- subset(asvITS.OCG, mdITS_OCG$Year=="2012") #105 of 1440 variables
asvITS.2012 <- asvITS.2012[,colSums(asvITS.2012) > 0]
summary(rowSums(asvITS.2012)) #507
summary(colSums(asvITS.2012)) #2.0

mdITS.2012 <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(asvITS.2012)) #105 of 16
mdITS.2012 <- droplevels(mdITS.2012)

#rarefy
asvITS.2012.r <- rrarefy(asvITS.2012,507) ## rarefy. warning message

#nmds
set.seed(12)
#asvITS.2012.nmds <- metaMDS(asvITS.2012.r, trymax=500) ###solution reached! warning message
#save(asvITS.2012.nmds, file = "nmds/asvITS.2012.nmds.rda")
load("nmds/asvITS.2012.nmds.rda")

ordiplot(asvITS.2012.nmds, type = "t",display = "sites",cex = .6)

plot(asvITS.2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2012 fungal community by subspecies", 
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.2012$Subspecies],
     pch=19)
legend("topleft", 
       legend=c("T","V","W"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subspecies, show.groups = "W", col = "goldenrod")

###PERMANOVA and adonis for 2012 subspecies####
asvITS.2012.ad <- adonis2(asvITS.2012.r ~ mdITS.2012$Subspecies) # Bray-Curtis is the default metric
asvITS.2012.ad #subspecies significant 0.008

#pairwiseadonis
asvITS.2012.subsp.pw <- pairwise.adonis(asvITS.2012.r, mdITS.2012$Subspecies)
asvITS.2012.subsp.pw #T vs V= 0.027, T vs W= 0.027, and W vs V= 1.000.

## 2012 NMDS for fungal community by ploidy ####
plot(asvITS.2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by ploidy", 
     col= c("red","blue")[mdITS.2012$Ploidy],
     pch=c(19,17)[mdITS.2012$Year])
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Ploidy, show.groups = "2n", col = "red")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Ploidy, show.groups = "4n", col = "blue")

### PERMANOVA for 2012 ploidy####
asvITS.2012.ploidy <- adonis2(asvITS.2012.r ~ mdITS.2012$Ploidy,by="margin") # Bray-Curtis is the default metric
asvITS.2012.ploidy #ploidy is not significant

## 2012 NMDS for fungal community by subsp ploidy ####
plot(asvITS.2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS.2012$Subsp_ploidy],
     pch=c(19,17)[mdITS.2012$Year])
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "W_4n", col = "purple")

### PERMANOVA and adonis for subspecies ploidy ####
asvITS.2012.subsp_ploi <- adonis2(asvITS.2012.r ~ mdITS.2012$Subsp_ploidy) 
asvITS.2012.subsp_ploi #subspecies ploidy is not sig

#pairwiseadonis
asvITS.2012.subsp.pw <- pairwise.adonis(asvITS.2012.r, mdITS.2012$Subsp_ploidy)
asvITS.2012.subsp.pw #none sig

## 2021 NMDS fungal community by subspecies####
asvITS.2021 <- subset(asvITS.OCG, mdITS_OCG$Year=="2021") #49 of 769 var 
asvITS.2021 <- asvITS.2021[,colSums(asvITS.2021) > 0]

summary(rowSums(asvITS.2021)) #648
summary(colSums(asvITS.2021)) #2

mdITS.2021 <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(asvITS.2021)) #49 obs of 16
mdITS.2021 <- droplevels(mdITS.2021)

#rarefying
asvITS.2021.r <- rrarefy(asvITS.2021,648) ## rarefy. warning message

#nmds
set.seed(78)
#asvITS.2021.nmds <- metaMDS(asvITS.2021.r, trymax=500) ###solution reached!
#save(asvITS.2021.nmds, file = "nmds/asvITS.2021.nmds.rda")
load("nmds/asvITS.2021.nmds.rda")

ordiplot(asvITS.2021.nmds, type = "t",display = "sites",cex = .6)

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

###PERMANOVA and adonis for 2021 subspecies####
asvITS.2021.ad <- adonis2(asvITS.2021.r ~ mdITS.2021$Subspecies) # Bray-Curtis is the default metric
asvITS.2021.ad #subspecies not significant

## 2021 NMDS for fungal community by ploidy ####
plot(asvITS.2021.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by ploidy", 
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

### PERMANOVA for 2021 ploidy####
asvITS.2021.ploidy <- adonis2(asvITS.2021.r ~ mdITS.2021$Ploidy,by="margin") # Bray-Curtis is the default metric
asvITS.2021.ploidy #ploidy is significant 0.024

## 2021 NMDS for fungal community by subsp ploidy ####
plot(asvITS.2021.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS.2021$Subsp_ploidy],
     pch=19)
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "W_4n", col = "purple")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "V_2n", col = "green")

#there is no V_2n in 2021

### PERMANOVA and adonis for subspecies ploidy ####
asvITS.2021.subsp_ploi <- adonis2(asvITS.2021.r ~ mdITS.2021$Subsp_ploidy) 
asvITS.2021.subsp_ploi #subspecies ploidy is not quite sig 0.059

#pairwiseadonis
asvITS.2021.subsp.pw <- pairwise.adonis(asvITS.2021.r, mdITS.2021$Subsp_ploidy)
asvITS.2021.subsp.pw #


# Mantel measuring correlation between the distance matrices ####
asvITS.2012.dist <- vegdist(asvITS.2012) 
OCG_AUC_2012.dist <- vegdist(OCG_AUC_2012_subset)

asv_GC.mant <- mantel(asvITS.2012.dist, OCG_AUC_2012.dist, permutations = 9999)



#Procrustes analyses####

## 2012 GC and fungal data asv####

load("nmds/asvITS.2012.nmds.rda") # asv 2012 nmds 
load("nmds/OCG_AUC_2012_ID.nmds.rda") # 2012 GC plant nmds 
load("nmds/OCG_AUC_2012_subset.nmds.rda") # GC 2012 compound nmds

# Subset the larger dataset (nmds1) to match the size of the smaller dataset (nmds2)
asvITS.2012_subset <- asvITS.2012[1:nrow(OCG_AUC_2012_subset), ]

#rarefy
asvITS.2012_subset.r <- rrarefy(asvITS.2012_subset,1046) ## rarefy. warning message
set.seed(7)
asvITS.2012_subset.nmds <- metaMDS(asvITS.2012_subset.r, trymax=500) ###solution reached! warning message
save(asvITS.2012_subset.nmds, file = "nmds/asvITS.2012_subset.nmds.rda")
load("nmds/asvITS.2012_subset.nmds.rda")


asv_GC.2012.pro <- protest(asvITS.2012_subset.nmds, Plant_ID_OCG_AUC_2012.nmds, symmetric=T) 
asv_GC.2012.pro ## Correlation in a symmetric Procrustes rotation: 0.07975, Significance: 0.958
summary(asv_GC.2012.pro) 
plot(asv_GC.2012.pro)

#GC_asv.2012.pro <- protest(Plant_ID_OCG_AUC_2012.nmds$points, asvITS.2012.nmds$points, symmetric=T)
#GC_asv.2012.pro #it is the same
#summary(GC_asv.2012.pro)

asv_GC.2012.pro$points #null

asv_GC_2012.prodat <- as.data.frame(asv_GC.2012.pro$X) #96obs of 2 variables
asv_GC_2012.prodat <- cbind(asv_GC.2012.pro, asv_GC.2012.pro$Yrot) #warning:In cbind(...):number of rows of result is not a multiple of vector length (arg 1)

colnames(asv_GC_2012.prodat)[colnames(asv_GC_2012.prodat)=="1"] <- "Xend" #NULL
colnames(asv_GC_2012.prodat)[colnames(asv_GC_2012.prodat)=="2"] <- "Yend"

ggplot() + 
  geom_segment(data=asv_GC_2012.prodat, mapping=aes(x=NMDS1, y=NMDS2, xend=Xend, yend=Yend)) + 
  geom_point(data=asv_GC_2012.prodat, mapping=aes(x=NMDS1, y=NMDS2)) +
  geom_point(data=asv_GC_2012.prodat, mapping=aes(x=Xend, y=Yend)) +
  labs(x="Procrustes axis 1", y="Procrustes axis 2") +
  theme_classic() 

##2021 GC and fungal procrustes ####
load("nmds/Plant_ID_OCG_AUC_2021.nmds.rda")
load("nmds/asvITS.2021.nmds.rda") # asv 2021 nmds

asv_GC.2021.pro <- protest(Plant_ID_OCG_AUC_2021.nmds, asvITS.2021.nmds, symmetric=T) 
asv_GC.2021.pro ## Correlation in a symmetric Procrustes rotation: 0.1325, Significance: 0.338
summary(asv_GC.2012.pro) 
plot(asv_GC.2012.pro)

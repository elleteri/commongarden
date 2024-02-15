# Orchard Common Garden - Microbial data cleaning, analysis, and visualizaition####
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
summary(rowSums(asvITS))
summary(colSums(asvITS)) 

##Full metadata read in####
#Data includes metadata from numerous projects. I will be subsetting for just the orchard common garden plants. The data will include plant ID number, location of origin, year of sampling (2012 or 2021), subspecies, ploidy, and subspecies ploidy.

mdITS <- read.csv("~/Documents/Orchard_Common _Garden/Shared_OCG_Code/data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #set to correct file path 463obs of 22 variables.
mdITS <- mdITS[order(row.names(mdITS)),] # order samples alphabetically
mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS))
colnames(asvITS) == row.names(mdITS) # sanity check: a check to make sure something does not contain elementary mistakes or impossibilities and is not based on invalid assumptions
#This sanity check reads true.

##Taxonomy table read in####
#taxonomy table is used to match to amplicon sequence variant table to fungal ID.
tax.ITS <- read.csv("~/Documents/Orchard_Common _Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #5983 obs of 2 variables
row.names(asvITS) == row.names(tax.ITS) # sanity check: TRUE
#Cleaning data####
##ASV table cleaning ####
#asvITS[asvITS < 10] <- 0 # each observation needs at least 10 seqs.
asvITS <- asvITS[rowSums(asvITS) > 0,] #showing me the values that are greater than zero
summary(rowSums(asvITS))
summary(colSums(asvITS)) 

asvITS <- asvITS[,colSums(asvITS) > 999] # each sample needs at least 1000 seqs. how low can I go for the ordiarrows to have them all have a pair. 

asvITS.t <- t(asvITS) # transpose rows and columns
asvITS.t <- asvITS.t[order(row.names(asvITS.t)),] # order samples alphabetically
asvITS.t <- asvITS.t[,order(colnames(asvITS.t))] # order asvs alphabetically

summary(rowSums(asvITS.t)) #making sure that it transposed
summary(colSums(asvITS.t))

mdITS <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.t)) #subsetting metadata to match what is left in the asv table

asvITS.t2 <- asvITS.t[!(row.names(asvITS.t) %in% c("NEG_9-30-21","AH1919","AHM20207","AHM20125")),] #outliers that need to be removed

asvITS.t2 <- asvITS.t2[,colSums(asvITS.t2) > 0] #keeping samples greater than 0

summary(rowSums(asvITS.t2)) #1021
summary(colSums(asvITS.t2)) #2.0

mdITS2 <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.t2)) # again subsetting metadata to match what is in the asv table

### Subsetting to just the plant in the common garden (OCG)####
asvITS.OCG <- subset(asvITS.t2, mdITS2$Project=="OCG") #subsetting to just the Orchard common garden plants
asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 0]

summary(rowSums(asvITS.OCG)) #1021
summary(colSums(asvITS.OCG)) #2.0

### Remove duplicates####
asvITS.OCG <- asvITS.OCG[grep("v2", row.names(asvITS.OCG), invert = T),]

## Remove negative control
asvITS.OCG <- asvITS.OCG[!(row.names(asvITS.OCG) == "NEG_8-28-21"),]

## Remove MTW.3.7.R_2012
asvITS.OCG <- asvITS.OCG[!(row.names(asvITS.OCG) == "MTW.3.7.R_2012"),]

## subset metadata to match asv table samples
asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 0]
summary(rowSums(asvITS.OCG)) #1021
summary(colSums(asvITS.OCG)) #2.0

mdITS.OCG <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.OCG)) #subset md to match asv table samples

#Write csv for cleaned metadata and asv table####
#write.csv(asvITS.OCG, file = "data_csv/asvITS.OCG.csv") #write csv for asv table of OCG only
#write.csv(mdITS.OCG, file = "data_csv/metadata_OCG.csv") #write csv for metadata of OCG only




#Clear Global Environment #############################


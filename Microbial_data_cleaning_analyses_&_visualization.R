# Orchard Common Garden - Microbial data cleaning, analysis, and visualization

# Set working directory and load necessary packages####
setwd("/Users/ellehorwath/Documents/Orchard_Common _Garden/commongarden")
if (!require("readr")) {install.packages("readr"); require("readr")}
if (!require("ggplot2")) {install.packages("ggplot2"); require("ggplot2")}
if (!require("dplyr")) {install.packages("dplyr"); require("dplyr")}
if (!require("tidyr")) {install.packages("tidyr"); require("tidyr")}
if (!require("MASS")) {install.packages("MASS"); require("MASS")}
if (!require("vegan")) {install.packages("vegan"); require("vegan")}
if (!require("effects")) {install.packages("effects"); require("effects")}
if (!require("metacoder")) {install.packages("metacoder"); require("metacoder")}
if (!require("exactRankTests")) {install.packages("exactRankTests"); require("exactRankTests")}
if (!require("nlme")) {install.packages("nlme"); require("nlme")}
if (!require("picante")) {install.packages("picante"); require("picante")}
if (!require("nVennR")) {install.packages("nVennR"); require("nVennR")}
if (!require("lme4")) {install.packages("lme4"); require("lme4")}
if (!require("reshape2")) {install.packages("reshape2"); require("reshape2")}
if (!require("devtools")) {install.packages("devtools"); require("devtools")}
if (!require("pairwiseAdonis")) {devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"); require("pairwiseAdonis")}
if (!require("iNEXT")) {install.packages("iNEXT"); require("iNEXT")}
if (!require("raster")) {install.packages("raster"); require("raster")}
if (!require("BiocManager")) {install.packages("BiocManager"); require("BiocManager")}
if (!require("phyloseq")) {BiocManager::install("phyloseq"); require("phyloseq")}
if (!require("qiime2R")) {BiocManager::install("qiime2R"); require("qiime2R")}

#Read in data : DONT START HERE ####
## ASV DATA
#Data has not been filtered and is not yet clean to include just observations with at least 10 seqs and each sample needs at least 1000 seq. The data has not been transposed and is ordered alphabetically. 

asvITS<- read.csv("data_csv/asv-table-dada2-ITS-sagebrush.csv",head=T,row.names=1, check.names = F) #5983 obs of 463 variable
asvITS<- asvITS[,order(colnames(asvITS))] # order samples alphabetically
summary(rowSums(asvITS)) #1.0
summary(colSums(asvITS)) #0

##METADATA
#Data includes metadata from numerous projects. I will be subsetting for just the orchard common garden plants. The data will include plant ID number, location of origin, year of sampling (2012 or 2021), subspecies, ploidy, and subspecies ploidy, status (dead or alive) in 2020.

mdITS <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #505 obs of 16 variables.
mdITS <- mdITS[order(row.names(mdITS)),]
mdITS <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS)) #463 of 16 variables
colnames(asvITS) == row.names(mdITS) # sanity check true.

mdITS$Description <- sub("_2020$", "_2021", mdITS$Description) # Replace "2020" with "2021" at the end of the strings

##TAXONOMY
#taxonomy table is used to match to amplicon sequence variant table to fungal ID.
tax.ITS <- read.csv("~/Documents/Orchard_Common_Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #5983 obs of 2 variables
row.names(asvITS) == row.names(tax.ITS) #TRUE


#Cleaning data####
##ASV 
#asvITS[asvITS < 10] <- 0 # each observation needs at least 10 seqs.
asvITS <- asvITS[rowSums(asvITS) > 0,] #the values that are greater than zero
summary(rowSums(asvITS)) #1
summary(colSums(asvITS)) #0

asvITS <- asvITS[,colSums(asvITS) > 499] # each sample needs at least 500 seqs. #5983 of 385 var

summary(colSums(asvITS)) #507
summary(rowSums(asvITS)) #0

asvITS.t <- t(asvITS) # transpose rows and columns
asvITS.t <- asvITS.t[order(row.names(asvITS.t)),] # order samples alphabetically
asvITS.t <- asvITS.t[,order(colnames(asvITS.t))] # order asvs alphabetically

summary(rowSums(asvITS.t)) 
summary(colSums(asvITS.t))

mdITS <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.t)) 

asvITS.t2 <- asvITS.t[!(row.names(asvITS.t) %in% c("NEG_9-30-21","AH1919","AHM20207","AHM20125","UTW.1.4_2021")),] #outliers removed

asvITS.t2 <- asvITS.t2[,colSums(asvITS.t2) > 0] #keeping samples greater than 0

summary(rowSums(asvITS.t2)) #507
summary(colSums(asvITS.t2)) #2.0

mdITS2 <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.t2)) #380 obs of 16 var

### Subsetting to just the plant in the common garden (OCG)
asvITS.OCG <- subset(asvITS.t2, mdITS2$Project=="OCG") 
asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 0]

summary(rowSums(asvITS.OCG)) #507
summary(colSums(asvITS.OCG)) #2.0

### Remove duplicates from ASV
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1')
asvITS.OCG <- asvITS.OCG[!rownames(asvITS.OCG) %in% rows_to_remove, ]

## Remove negative control
asvITS.OCG <- asvITS.OCG[!(row.names(asvITS.OCG) == "NEG_8-28-21"),]

## Remove MTW.3.7.R_2012
asvITS.OCG <- asvITS.OCG[!(row.names(asvITS.OCG) == "MTW.3.7.R_2012"),] 

asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 0]
summary(rowSums(asvITS.OCG)) #507
summary(colSums(asvITS.OCG)) #2.0

mdITS.OCG <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.OCG)) ##154 of 16 var

#Write csv for cleaned metadata and asv table
write.csv(asvITS.OCG, file = "data_csv/asvITS.OCG.csv") 
write.csv(mdITS.OCG, file = "data_csv/metadata_OCG.csv") 

#Clear Global Environment 
rm(list = ls())

#Read in cleaned data : START HERE ####
#METADATA
mdITS.OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #154 of 16 variables
head(mdITS.OCG)
str(mdITS.OCG)

#ASV
asvITS.OCG <- read.csv("data_csv/asvITS.OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #154 obs of 2135 variables
head(asvITS.OCG)
summary(rowSums(asvITS.OCG)) #507
summary(colSums(asvITS.OCG)) #2

row.names(asvITS.OCG) == row.names(mdITS.OCG) # sanity check:TRUE

#TAXONOMY
tax.ITS <- read.csv("~/Documents/Orchard_Common_Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #taxonomy read in

#Alpha diversity OCG####
##Rarefying
asvITS.OCG.r <- rrarefy(asvITS.OCG,507) ## rarefy: Warning message
asvITS.OCG.shannon <- diversity(asvITS.OCG.r)
asvITS.OCG.ef <- exp(asvITS.OCG.shannon)
asvITS.OCG.efr <- round(asvITS.OCG.ef)
mdITS.OCG <- cbind(mdITS.OCG, effective_species = asvITS.OCG.efr)

glm.OCG <- glm(effective_species ~ Subspecies + Year, family = poisson, data=mdITS.OCG)
summary(glm.OCG)

plot(allEffects(glm.OCG))

plot(mdITS.OCG$Year,asvITS.OCG.ef)
plot(mdITS.OCG$Subspecies,asvITS.OCG.ef)

ggplot(mdITS.OCG, aes(Year, effective_species))+
  geom_boxplot(aes(group = Year, fill = Year))+
  theme_classic()

ggplot(data = mdITS.OCG, mapping = aes(x = Subspecies, y = effective_species, fill = Subspecies)) +
  geom_boxplot() +
  theme_classic()

#Beta diversity#### 
#NMDS plots
set.seed(41)
#asvITS.OCG.nmds <- metaMDS(asvITS.OCG.r, trymax=500) ### Solution reached! gives me a warning
# save(asvITS.OCG.nmds, file = "nmds/asvITS.OCG_nmds.rda") #save the nmds so you won't need to run it again
load("nmds/asvITS.OCG_nmds.rda") #load it to use in code anytime after the initial run

ordiplot(asvITS.OCG.nmds, type = "t",display = "sites",cex = .6)

mdITS.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant")] <- lapply(mdITS.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant")], as.factor)
str(mdITS.OCG)

#SUBSPECIES
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="Sagebrush fungal community by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.OCG$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subspecies, show.groups = "W", col = "goldenrod")

### PERMANOVA and adonis for subspecies ##
asvITS.OCG.subsp <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Subspecies,by="margin") 
asvITS.OCG.subsp #subspecies significant

#pairwiseadonis
asvITS.OCG.subsp.pw <- pairwise.adonis(asvITS.OCG.r, mdITS.OCG$Subspecies)
asvITS.OCG.subsp.pw #T vs V= 0.003, T vs W= 0.075, and W vs V= 0.540.

#PLOIDY
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by ploidy", 
     col= c("red","blue")[mdITS.OCG$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Ploidy, show.groups = "2n", col = "red")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Ploidy, show.groups = "4n", col = "blue")

### PERMANOVA for ploidy##
asvITS.OCG.ploidy <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Ploidy,by="margin") 
asvITS.OCG.ploidy #ploidy is not significant

#SUBSPECIES PLOIDY
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS.OCG$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subsp_ploidy, show.groups = "W_4n", col = "purple")

### PERMANOVA and adonis for subspecies ploidy ##
asvITS.OCG.subsp_ploi <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Subsp_ploidy) 
asvITS.OCG.subsp_ploi #subspecies ploidy is significant 0.002

#pairwiseadonis
asvITS.OCG.subsp.pw <- pairwise.adonis(asvITS.OCG.r, mdITS.OCG$Subsp_ploidy)
asvITS.OCG.subsp.pw #p-adjusted:T_4n vs V_2n = 0.04, T_4n vs V_4n = 0.03, T_4n vs W_4n = 0.05

#YEAR
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by year", 
     col= rainbow(2)[mdITS.OCG$Year],
     pch=19)
# text(asvITS.OCG.nmds$points[, 1], 
#      asvITS.OCG.nmds$points[, 2], 
#      labels = mdITS.OCG$Plant, 
#      pos = 3, 
#      cex = 0.5)
legend("topleft", 
       legend=c("2012","2021"),
       col= c("red","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")
ordiarrows(asvITS.OCG.nmds, mdITS.OCG$Plant)

### PERMANOVAs for year ##
asvITS.OCG.subsp_loc_yr <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Subspecies + mdITS.OCG$Year + mdITS.OCG$Location, by = "margin") 
asvITS.OCG.subsp_loc_yr #Year is significant

asvITS.OCG.subsp_yr <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Subspecies*mdITS.OCG$Year) 
asvITS.OCG.subsp_yr #year and subspecies are sig

asvITS.OCG.yr <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Year) 
asvITS.OCG.yr #year is significant

## 2012 
#NMDS
asvITS.2012 <- subset(asvITS.OCG, mdITS.OCG$Year=="2012") #105 of 1440 variables
asvITS.2012 <- asvITS.2012[,colSums(asvITS.2012) > 0]
summary(rowSums(asvITS.2012)) #507
summary(colSums(asvITS.2012)) #2.0

mdITS.2012 <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(asvITS.2012)) #105 of 17
str(mdITS.2012)

#rarefy
asvITS.2012.r <- rrarefy(asvITS.2012,507) ## rarefy. warning message

set.seed(12)
#asvITS.2012.nmds <- metaMDS(asvITS.2012.r, trymax=500) ###solution reached! warning message
#save(asvITS.2012.nmds, file = "nmds/asvITS.2012.nmds.rda")
load("nmds/asvITS.2012.nmds.rda")

ordiplot(asvITS.2012.nmds, type = "t",display = "sites",cex = .6)

#SUBSPECIES
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

###PERMANOVA and adonis for 2012 subspecies##
asvITS.2012.ad <- adonis2(asvITS.2012.r ~ mdITS.2012$Subspecies)
asvITS.2012.ad #subspecies significant 0.008

#pairwiseadonis
asvITS.2012.subsp.pw <- pairwise.adonis(asvITS.2012.r, mdITS.2012$Subspecies)
asvITS.2012.subsp.pw #T vs V= 0.027, T vs W= 0.027, and W vs V= 1.000.

#PLOIDY
plot(asvITS.2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2012 fungal community by ploidy", 
     col= c("red","blue")[mdITS.2012$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Ploidy, show.groups = "2n", col = "red")
ordispider(asvITS.2012.nmds,groups = mdITS.2012$Ploidy, show.groups = "4n", col = "blue")

### PERMANOVA for 2012 ploidy##
asvITS.2012.ploidy <- adonis2(asvITS.2012.r ~ mdITS.2012$Ploidy,by="margin") 
asvITS.2012.ploidy #ploidy is not significant

#SUBSP PLOIDY
plot(asvITS.2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2012 fungal community by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS.2012$Subsp_ploidy],
     pch=c(19))
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

### PERMANOVA and adonis for subspecies ploidy ##
asvITS.2012.subsp_ploi <- adonis2(asvITS.2012.r ~ mdITS.2012$Subsp_ploidy) 
asvITS.2012.subsp_ploi #subspecies ploidy is sig

#pairwiseadonis
asvITS.2012.subsp.pw <- pairwise.adonis(asvITS.2012.r, mdITS.2012$Subsp_ploidy)
asvITS.2012.subsp.pw #none sig

## 2021
asvITS.2021 <- subset(asvITS.OCG, mdITS.OCG$Year=="2021") #49 of 769 var 
asvITS.2021 <- asvITS.2021[,colSums(asvITS.2021) > 0]

summary(rowSums(asvITS.2021)) #648
summary(colSums(asvITS.2021)) #2

mdITS.2021 <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(asvITS.2021)) #49 obs of 17

#rarefying
asvITS.2021.r <- rrarefy(asvITS.2021,648) ## rarefy. warning message

#NMDS
set.seed(78)
#asvITS.2021.nmds <- metaMDS(asvITS.2021.r, trymax=500) ###solution reached!
#save(asvITS.2021.nmds, file = "nmds/asvITS.2021.nmds.rda")
load("nmds/asvITS.2021.nmds.rda")

ordiplot(asvITS.2021.nmds, type = "t",display = "sites",cex = .6)

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
asvITS.2021.ad <- adonis2(asvITS.2021.r ~ mdITS.2021$Subspecies) 
asvITS.2021.ad #subspecies not significant

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
asvITS.2021.ploidy <- adonis2(asvITS.2021.r ~ mdITS.2021$Ploidy,by="margin") 
asvITS.2021.ploidy #ploidy is significant 0.024

#SUBSP PLOIDY
plot(asvITS.2021.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2021 fungal community by subspecies and ploidy", 
     col= c("red","orange","cyan","purple")[mdITS.2021$Subsp_ploidy],
     pch=19)
legend("topleft", 
       legend=c("T_2n","T_4n","V_4n","W_4n"),
       col= c("red","orange","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(asvITS.2021.nmds,groups = mdITS.2021$Subsp_ploidy, show.groups = "W_4n", col = "purple")

#there is no V_2n in 2021

### PERMANOVA and adonis for subspecies ploidy ##
asvITS.2021.subsp_ploi <- adonis2(asvITS.2021.r ~ mdITS.2021$Subsp_ploidy) 
asvITS.2021.subsp_ploi #subspecies ploidy is not quite sig 0.059

#pairwiseadonis
asvITS.2021.subsp.pw <- pairwise.adonis(asvITS.2021.r, mdITS.2021$Subsp_ploidy)
asvITS.2021.subsp.pw #none sig

#Bar chart of ASV level ####
asvITSt <- t(asvITS.OCG.r)
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
tax.ITS.OCG.p <- subset(tax.ITS.p, rownames(tax.ITS.p) %in% rownames(t(asvITS.OCG.r)))

##### Class level OCG
asvITS.OCGtc <- data.frame(Class=tax.ITS.OCG.p$Class,t(asvITS.OCG.r))
asvITS.OCGtc$Class[is.na(asvITS.OCGtc$Class)] <- "Unknown"
asvITS.OCGtca <- aggregate(. ~ asvITS.OCGtc$Class, asvITS.OCGtc[,2:ncol(asvITS.OCGtc)], sum) 
row.names(asvITS.OCGtca) <- asvITS.OCGtca[,1]
asvITS.OCGtca <- asvITS.OCGtca[,2:ncol(asvITS.OCGtca)]

asvITS.OCGtcao <- as.matrix(asvITS.OCGtca[order(rowSums(asvITS.OCGtca),decreasing = T),])

customcol28 <- c("cadetblue4","royalblue3","darkblue","tomato1","dodgerblue2",
                 "cyan","darkred","purple","mediumblue","palegoldenrod",
                 "lightgoldenrod","indianred","yellow","purple4","darkgreen",
                 "lightsalmon","yellow3","purple2","lightblue","firebrick",
                 "navy","red4","red","darkmagenta","mediumvioletred",
                 "violetred2","skyblue","dodgerblue4")

barplot(asvITS.OCGtcao,col=customcol28,legend.text=F,axes=F,cex.names = .3,las=2, args.legend = list(x = "topleft", bty = "n", inset=c(-0.15, 0)))
barplot(asvITS.OCGtcao,col=customcol28,legend.text=T,axes=F,cex.names = .3,las=2, args.legend = list(x = "topleft", bty = "n", inset=c(-0.11, -0.1), cex=0.4))

#Betadispersion ####
#YEAR
modyr <- betadisper(vegdist(asvITS.OCG.r), factor(mdITS.OCG$Year, labels = c("2012","2021")))
anova(modyr) #not sig between years p = 0.064

#SUBSPECIES
modsp <- betadisper(vegdist(asvITS.OCG.r), factor(mdITS.OCG$Subspecies, labels = c("T","V","W")))
anova(modsp) #no sig between years. p = 0.097

#SUBSPECIES PLOIDY
modsubplo <- betadisper(vegdist(asvITS.OCG.r), factor(mdITS.OCG$Subsp_ploidy, labels = c("T_2n","T_4n","V_2n","V_4n","W_4n")))
anova(modsubplo) #sig = 0.004591
boxplot(modsubplo)

ordiplot(modsubplo, type = "t",display = "sites",cex = .6)
plot(modsubplo,ellipse = TRUE,hull = FALSE ,conf = 0.80)#90% data ellipse

permutest(modsubplo, permutations = 99,pairwise = TRUE)

modsubplo.HSD<-TukeyHSD(modsubplo) #sig difference detected
plot(modsubplo.HSD)

#Clear Global environment
rm(list = )

#ANCOM####
#Read in project data 
#ASV
asvITS <- read.csv("data_csv/asv-table-dada2-ITS-sagebrush.csv",head=T,row.names=1, check.names = F) 
asvITS <- asvITS[,order(colnames(asvITS))] # order samples alphabetically

#METADATA
mdITS <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) 
mdITS <- mdITS[order(row.names(mdITS)),] 
mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS)) 
colnames(asvITS) == row.names(mdITS) # sanity check:TRUE

#TAXONOMY
tax.ITS <- read.csv("~/Documents/Orchard_Common_Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #5983 obs of 2 variables
row.names(asvITS) == row.names(tax.ITS) #TRUE

#Cleaning
##Remove entries with insufficient sequences
asvITS[asvITS < 10] <- 0
asvITS <- asvITS[rowSums(asvITS) > 0,] # each observation needs at least 10 seqs

summary(rowSums(asvITS)) #10
summary(colSums(asvITS)) #0

asvITS <- asvITS[,colSums(asvITS) > 999] # each sample needs at least 1000 seqs

mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS)) # remove trimmed samples from metadata

mdITS.OCG <- subset(mdITS, mdITS$Project == "OCG")   
mdITS.OCG <- subset(mdITS.OCG, mdITS.OCG$Location != "NEG")  
asvITS.t <- t(asvITS) 
asvITS.OCG_t <- subset(asvITS.t, row.names(asvITS.t) %in% row.names(mdITS.OCG)) 
asvITS.OCG <- t(asvITS.OCG_t)
asvITS.OCG <- as.data.frame(asvITS.OCG) 

asvITS.OCG[asvITS.OCG < 10] <- 0  # repeat cleaning after trimming
asvITS.OCG <- asvITS.OCG[rowSums(asvITS.OCG) > 0,] # each observation needs at least 10 seqs

summary(rowSums(asvITS.OCG)) #10
summary(colSums(asvITS.OCG)) #1018

asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 999] # each sample needs at least 1000 seqs

mdITS.OCG <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% colnames(asvITS.OCG)) 
colnames(asvITS.OCG) == row.names(mdITS.OCG)# sanity check :TRUE

taxITS.OCG <- subset(tax.ITS, row.names(tax.ITS) %in% row.names(asvITS.OCG))
row.names(taxITS.OCG) == row.names(asvITS.OCG)  # sanity check: TRUE

#Create ANCOM Function ####

ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL,check.names=FALSE)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T & repeated==F ){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  if(adjusted==T & repeated==T ){
    adjusted.formula = paste0(base.formula," + ", adj.formula," | ", repeat.var)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
}



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}


#####End of functions

#Create "Sample.ID" column for all data tables
#ANCOM requires that data be formatted so that first *column* is named "Sample.ID"
#Sample IDs as row names does not count!
asvITS.OCG_t <- t(asvITS.OCG)                                                                
taxITS.OCG_t <- t(taxITS.OCG) #transpose ASV and taxonomy data frames so sample IDs become row names
asvITS.OCG_sbst <- data.frame("Sample.ID" = row.names(asvITS.OCG_t), asvITS.OCG_t, check.names = F)#create new column in each of ASV, metadata, and taxonomy data frames labelled "Sample.ID" containing sample IDs
mdITS.OCG_sbst <- data.frame("Sample.ID" = row.names(mdITS.OCG), mdITS.OCG)
taxITS.OCG_sbst <- data.frame("Sample.ID" = row.names(taxITS.OCG_t), taxITS.OCG_t, check.names = F)

row.names(asvITS.OCG_sbst) == row.names(mdITS.OCG_sbst) #TRUE

#SUBSPECIES 
ANCOM_subspecies <- ANCOM.main(asvITS.OCG_sbst,mdITS.OCG_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)
#Create objects of significant ASVs
sigASVs_subspecies <- subset(ANCOM_subspecies$W.taxa, ANCOM_subspecies$W.taxa$W_stat > 0)[,1]
sigASVs_subspecies <- as.data.frame(sigASVs_subspecies) #convert significant ASV return to data frame
row.names(sigASVs_subspecies) <- sigASVs_subspecies[1:9,1] #change row names to ASV sequence

#By default ANCOM sorts ASVs in order of decreasing abundance. This section instead aims to sort by decreasing abundance, which can be cross-referenced with the .csv generated after the break
#If you want to look at the default ANCOM output, the full code will be included at the end of this section

sigASVs_subspecies[,1] <- c(1:9) #append numerical ranking of significance in additional column
sigASVs_subspecies_t <- t(sigASVs_subspecies) #some transposition and fiddling to get ASV sequences as usable column names to mimic format of original data frames
sigASVs_subspecies_t <- as.data.frame(sigASVs_subspecies_t)
colnames(sigASVs_subspecies_t) <- as.character(colnames(sigASVs_subspecies_t))                   
print(colnames(sigASVs_subspecies_t))                                                            
sigASVs_subspecies_t <- sigASVs_subspecies_t[,order(colnames(sigASVs_subspecies_t))]#sort alphabetically by ASV sequence, this is how we will align with original ASV data frame subsets
rownames(sigASVs_subspecies_t) <- c("sig_rank")#rename row to something more informative

write.csv(ANCOM_subspecies$W.taxa, file = "data_csv/ANCOM/ANCOM_subspecies.csv") #write ANCOM output to .csv to observe significance cutoffs

#Subset original data by significant ASVs
asvITS.OCG_sigsbst_subspecies <-  t(subset(t(asvITS.OCG_sbst), colnames(asvITS.OCG_sbst) %in% row.names(sigASVs_subspecies))) #subset of ASV counts to only significant taxa
taxITS.OCG_sigsbst_subspecies <-  t(subset(t(taxITS.OCG_sbst), colnames(taxITS.OCG_sbst) %in% row.names(sigASVs_subspecies))) #subset of taxonomy to only significant taxa

#This next section is a continuation of the above work to sort ASVs by decreasing significance
asvITS.OCG_sigsbst_subspecies <- asvITS.OCG_sigsbst_subspecies[,order(colnames(asvITS.OCG_sigsbst_subspecies))] #sort alphabetically for alignment with significance ranking
colnames(sigASVs_subspecies_t) == colnames(asvITS.OCG_sigsbst_subspecies) #sanity check:TRUE

asvITS.OCG_sigsbst_subspecies <- rbind(asvITS.OCG_sigsbst_subspecies, sigASVs_subspecies_t)#attach significance rankings
asvITS.OCG_sigsbst_subspecies_t <- as.data.frame(t(asvITS.OCG_sigsbst_subspecies))#transpose data frame to prepare for sorting
asvITS.OCG_sigsbst_subspecies_t$sig_rank <- as.numeric(asvITS.OCG_sigsbst_subspecies_t$sig_rank) #transform significance rankings to numeric (else 10 comes before 2)
asvITS.OCG_sigsbst_subspecies_t <- asvITS.OCG_sigsbst_subspecies_t[order(asvITS.OCG_sigsbst_subspecies_t$sig_rank),] #sort by significance ranking
asvITS.OCG_sigsbst_subspecies_t <- subset(asvITS.OCG_sigsbst_subspecies_t, select=-c(sig_rank)) #remove rank column
asvITS.OCG_sigsbst_subspecies_t <- as.data.frame(t(asvITS.OCG_sigsbst_subspecies_t))#transpose back and convert to data frame

#Build objects for plotting
sbsplot <- data.frame(asvITS.OCG_sigsbst_subspecies_t[,1:3], "subspecies" = mdITS.OCG_sbst$Subspecies, check.names = FALSE)
sbsplot[,1:3] <- lapply(sbsplot[,1:3], function(x) as.numeric(as.character(x)))
sbsplot[,1:3] <- log(sbsplot[,1:3]+1)
sbsplot_sub <- data.frame(sample=rownames(sbsplot),sbsplot, check.names = F)
sbsplot_sublong <- melt(sbsplot_sub)

ANCOM_subspecies$W.taxa

sbsplot_sublong$subspecies <-  factor(sbsplot_sublong$subspecies, levels = c("T", "V", "W"))

#SUBSPECIES FIGURE TOP ANCOM ASVs
ggplot(sbsplot_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_colour_manual('Species',labels=c('Dothidea sambuci', 'Endoconidioma populi','Parastagonospora novozelandica'),values=c("darkkhaki", "palevioletred", "darkcyan"))+
  ylab("Log  rel. abundance") + xlab("Subspecies") + theme_classic() #Highest in tridentata and then two are very low in vaseyana and one low in wyomingensis. Look it up in the taxonomy in csv and then BLAST it. Dothidea sambuci=1bfe6883c827ae962725d53d567780fc (97.95% identity), 89728a516cec30e6a6c3b1e428ca5f40= Endoconidioma populi (99.60%), 3518c430b3f2d2c3844aee6b7d76990a=	Parastagonospora novozelandica (95.62%)

print(taxITS.OCG_sigsbst_subspecies)

###SIGNIFICANCE BY SUBSPECIES - Default ANCOM sorting 
#Run ANCOM, specify variable
ANCOM_subspecies_default <- ANCOM.main(asvITS.OCG_sbst,mdITS.OCG_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)
#Create objects of significant ASVs
sigASVs_subspecies_default <- subset(ANCOM_subspecies_default$W.taxa, ANCOM_subspecies_default$W.taxa$W_stat > 0)[,1]
sigASVs_subspecies_default <- as.data.frame(sigASVs_subspecies_default) #convert significant ASV return to data frame
row.names(sigASVs_subspecies_default) <- sigASVs_subspecies_default[1:9,1] #change row names to ASV sequence

#Subset original data by significant ASVs
asvITS.OCG_sigsbst_subspecies_default <-  t(subset(t(asvITS.OCG_sbst), colnames(asvITS.OCG_sbst) %in% row.names(sigASVs_subspecies_default))) #subset of ASV counts to only significant taxa
taxITS.OCG_sigsbst_subspecies_default <-  t(subset(t(taxITS.OCG_sbst), colnames(taxITS.OCG_sbst) %in% row.names(sigASVs_subspecies_default))) #subset of taxonomy to only significant taxa

#Build objects for plotting
sbsdfplot <- data.frame(asvITS.OCG_sigsbst_subspecies_default[,1:9], "subspecies" = mdITS.OCG_sbst$Subspecies, check.names = FALSE)
sbsdfplot[,1:9] <- lapply(sbsdfplot[,1:9], function(x) as.numeric(as.character(x)))
sbsdfplot[,1:9] <- log(sbsdfplot[,1:9]+1)

sbsdfplot_sub <- data.frame(sample=rownames(sbsdfplot),sbsdfplot, check.names = F)
sbsdfplot_sublong <- melt(sbsdfplot_sub)

sbsdfplot_sublong$subspecies <-  factor(sbsdfplot_sublong$subspecies, levels = c("T", "V", "W"))

### Figure - OCG top ANCOM ASVs by Subspecies with default sorting 
ggplot(sbsdfplot_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable, shape = variable), alpha =.4) +
  ylab("Log  rel. abundance") + xlab("Subspecies") + theme_classic()

print(taxITS.OCG_sigsbst_subspecies_default)
print(taxITS.OCG_sisubspeciesgsbst_subspecies_default)

#PLOIDY - no sig differences found
#Run ANCOM, specify variable
ANCOM_ploidy <- ANCOM.main(asvITS.OCG_sbst,mdITS.OCG_sbst,F,F,"Ploidy",NULL,NULL,F,NULL,2,.05,.9)
#Create objects of significant ASVs
sigASVs_ploidy <- subset(ANCOM_ploidy$W.taxa, ANCOM_ploidy$W.taxa$W_stat > 0)[,1] #this is empty, as there are no significant differences found
# sigASVs_ploidy <- as.data.frame(sigASVs_ploidy)
# row.names(sigASVs_ploidy) <- sigASVs_ploidy[1:10,1]
#sigASVs_ploidy[,1] <- c(1:10)
#sigASVs_ploidy_t <- t(sigASVs_ploidy)
#sigASVs_ploidy_t <- as.data.frame(sigASVs_ploidy_t)
#colnames(sigASVs_ploidy_t) <- as.character(colnames(sigASVs_ploidy_t))
#print(colnames(sigASVs_ploidy_t))
#sigASVs_ploidy_t <- sigASVs_ploidy_t[,order(colnames(sigASVs_ploidy_t))]
#rownames(sigASVs_ploidy_t) <- c("sig_rank")

write.csv(ANCOM_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_ploidy.csv")

#Subset data by sig ASVs
#asvA_sigsbst_ploidy <-  t(subset(t(asvA_sbst), colnames(asvA_sbst) %in% row.names(sigASVs_ploidy)))
#taxA_sigsbst_ploidy <-  t(subset(t(taxA_sbst), colnames(taxA_sbst) %in% row.names(sigASVs_ploidy)))

#asvA_sigsbst_ploidy <- asvA_sigsbst_ploidy[,order(colnames(asvA_sigsbst_ploidy))]
#colnames(sigASVs_ploidy_t) == colnames(asvA_sigsbst_ploidy)
#asvA_sigsbst_ploidy <- rbind(asvA_sigsbst_ploidy, sigASVs_ploidy_t)
#asvA_sigsbst_ploidy_t <- as.data.frame(t(asvA_sigsbst_ploidy))
#asvA_sigsbst_ploidy_t$sig_rank <- as.numeric(asvA_sigsbst_ploidy_t$sig_rank)
#asvA_sigsbst_ploidy_t <- asvA_sigsbst_ploidy_t[order(asvA_sigsbst_ploidy_t$sig_rank),]
#asvA_sigsbst_ploidy_t <- subset(asvA_sigsbst_ploidy_t, select=-c(sig_rank))
#asvA_sigsbst_ploidy <- as.data.frame(t(asvA_sigsbst_ploidy_t))

#plplot <- data.frame(asvA_sigsbst_ploidy[,1:10], "ploidy" = metaA_sbst$Ploidy, check.names = FALSE)
#plplot[,1:10] <- lapply(plplot[,1:10], function(x) as.numeric(as.character(x)))
#plplot[,1:10] <- log(plplot[,1:10]+1)

#plplot_sub <- data.frame(sample=rownames(plplot),plplot, check.names = F)
#plplot_sublong <- melt(plplot_sub)

#plplot_sublong$ploidy <-  factor(plplot_sublong$ploidy, levels = c("2n", "4n"))

### (DNE) Figure - OCG top ANCOM ASVs by ploidy 
# ggplot(plplot_sublong, aes(y = value, x = ploidy, color=variable))+
#   geom_boxplot(outlier.shape = NA) + 
#   geom_point(position=position_dodge(width=0.75), aes(group=variable, shape = variable), alpha =.4) +
#   ylab("Log  rel. abundance") + xlab("Ploidy") + theme_classic()
# 
# #print(taxA_sigsbst_ploidy)

##SIGNIFICANCE BY SUBSPECIES PLOIDY
#Run ANCOM, specify variable
ANCOM_subsp_ploidy <- ANCOM.main(asvITS.OCG_sbst,mdITS.OCG_sbst,F,F,"Subsp_ploidy",NULL,NULL,F,NULL,2,.05,.9)
#Create objects of significant ASVs
sigASVs_subsp_ploidy <- subset(ANCOM_subsp_ploidy$W.taxa, ANCOM_subsp_ploidy$W.taxa$W_stat > 0)[,1]
sigASVs_subsp_ploidy <- as.data.frame(sigASVs_subsp_ploidy)
row.names(sigASVs_subsp_ploidy) <- sigASVs_subsp_ploidy[1:9,1]
sigASVs_subsp_ploidy[,1] <- c(1:9)
sigASVs_subsp_ploidy_t <- t(sigASVs_subsp_ploidy)
sigASVs_subsp_ploidy_t <- as.data.frame(sigASVs_subsp_ploidy_t)
colnames(sigASVs_subsp_ploidy_t) <- as.character(colnames(sigASVs_subsp_ploidy_t))
print(colnames(sigASVs_subsp_ploidy_t))
sigASVs_subsp_ploidy_t <- sigASVs_subsp_ploidy_t[,order(colnames(sigASVs_subsp_ploidy_t))]
rownames(sigASVs_subsp_ploidy_t) <- c("sig_rank")

write.csv(ANCOM_subsp_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_subsp_ploidy.csv")

#Subset data by sig ASVs
asvITS.OCG_sigsbst_subsp_ploidy <-  t(subset(t(asvITS.OCG_sbst), colnames(asvITS.OCG_sbst) %in% row.names(sigASVs_subsp_ploidy)))
taxITS.OCG_sigsbst_subsp_ploidy <-  t(subset(t(taxITS.OCG_sbst), colnames(taxITS.OCG_sbst) %in% row.names(sigASVs_subsp_ploidy)))

asvITS.OCG_sigsbst_subsp_ploidy <- asvITS.OCG_sigsbst_subsp_ploidy[,order(colnames(asvITS.OCG_sigsbst_subsp_ploidy))]
colnames(sigASVs_subsp_ploidy_t) == colnames(asvITS.OCG_sigsbst_subsp_ploidy) #TRUE
asvITS.OCG_sigsbst_subsp_ploidy <- rbind(asvITS.OCG_sigsbst_subsp_ploidy, sigASVs_subsp_ploidy_t)
asvITS.OCG_sigsbst_subsp_ploidy_t <- as.data.frame(t(asvITS.OCG_sigsbst_subsp_ploidy))
asvITS.OCG_sigsbst_subsp_ploidy_t$sig_rank <- as.numeric(asvITS.OCG_sigsbst_subsp_ploidy_t$sig_rank)
asvITS.OCG_sigsbst_subsp_ploidy_t <- asvITS.OCG_sigsbst_subsp_ploidy_t[order(asvITS.OCG_sigsbst_subsp_ploidy_t$sig_rank),]
asvITS.OCG_sigsbst_subsp_ploidy_t <- subset(asvITS.OCG_sigsbst_subsp_ploidy_t, select=-c(sig_rank))
asvITS.OCG_sigsbst_subsp_ploidy <- as.data.frame(t(asvITS.OCG_sigsbst_subsp_ploidy_t))

sbsplplot <- data.frame(asvITS.OCG_sigsbst_subsp_ploidy[,1:9], "subsp_ploidy" = mdITS.OCG_sbst$Subsp_ploidy, check.names = FALSE) 

sbsplplot[,1:9] <- lapply(sbsplplot[,1:9], function(x) as.numeric(as.character(x)))
sbsplplot[,1:9] <- log(sbsplplot[,1:9]+1)

sbsplplot_sub <- data.frame(sample=rownames(sbsplplot),sbsplplot, check.names = F)
sbsplplot_sublong <- melt(sbsplplot_sub)

sbsplplot_sublong$subsp_ploidy <-  factor(sbsplplot_sublong$subsp_ploidy, levels = c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))

### FIGURE FOR SUBSPECIES PLOIDY
ggplot(sbsplplot_sublong, aes(y = value, x = subsp_ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable, shape = variable), alpha =.4) +
  ylab("Log  rel. abundance") + xlab("Subsp_ploidy") + theme_classic()

print(taxA_sigsbst_subsp_ploidy)

##SIGNIFICANCE BY YEAR
#Run ANCOM, specify variable
ANCOM_year <- ANCOM.main(asvITS.OCG_sbst,mdITS.OCG_sbst,F,F,"Year",NULL,NULL,F,NULL,2,.05,.9)
#Create objects of significant ASVs
sigASVs_year <- subset(ANCOM_year$W.taxa, ANCOM_year$W.taxa$W_stat > 0)[,1]
sigASVs_year <- as.data.frame(sigASVs_year)
row.names(sigASVs_year) <- sigASVs_year[1:10,1]
sigASVs_year[,1] <- c(1:10)
sigASVs_year_t <- t(sigASVs_year)
sigASVs_year_t <- as.data.frame(sigASVs_year_t)
colnames(sigASVs_year_t) <- as.character(colnames(sigASVs_year_t))
print(colnames(sigASVs_year_t))
sigASVs_year_t <- sigASVs_year_t[,order(colnames(sigASVs_year_t))]
rownames(sigASVs_year_t) <- c("sig_rank")

write.csv(ANCOM_year$W.taxa, file = "data_csv/ANCOM/ANCOM_Year.csv")

#Subset data by sig ASVs
asvITS.OCG_sigsbst_year <-  t(subset(t(asvITS.OCG_sbst), colnames(asvITS.OCG_sbst) %in% row.names(sigASVs_year)))
taxITS.OCG_sigsbst_year <-  t(subset(t(taxITS.OCG_sbst), colnames(taxITS.OCG_sbst) %in% row.names(sigASVs_year)))

asvITS.OCG_sigsbst_year <- asvITS.OCG_sigsbst_year[,order(colnames(asvITS.OCG_sigsbst_year))]
colnames(sigASVs_year_t) == colnames(asvITS.OCG_sigsbst_year) #TRUE
asvITS.OCG_sigsbst_year <- rbind(asvITS.OCG_sigsbst_year, sigASVs_year_t)
asvITS.OCG_sigsbst_year_t <- as.data.frame(t(asvITS.OCG_sigsbst_year))
asvITS.OCG_sigsbst_year_t$sig_rank <- as.numeric(asvITS.OCG_sigsbst_year_t$sig_rank)
asvITS.OCG_sigsbst_year_t <- asvITS.OCG_sigsbst_year_t[order(asvITS.OCG_sigsbst_year_t$sig_rank),]
asvITS.OCG_sigsbst_year_t <- subset(asvITS.OCG_sigsbst_year_t, select=-c(sig_rank))
asvITS.OCG_sigsbst_year <- as.data.frame(t(asvITS.OCG_sigsbst_year_t))

yrplot <- data.frame(asvITS.OCG_sigsbst_year[,1:10], "year" = mdITS.OCG_sbst$Year, check.names = FALSE)
yrplot[,1:10] <- lapply(yrplot[,1:10], function(x) as.numeric(as.character(x)))
yrplot[,1:10] <- log(yrplot[,1:10]+1)

yrplot_sub <- data.frame(sample=rownames(yrplot),yrplot, check.names = F)
yrplot_sublong <- melt(yrplot_sub)

yrplot_sublong$year <-  factor(yrplot_sublong$year, levels = c("2012", "2021"))

### FIGURE FOR YEAR
ggplot(yrplot_sublong, aes(y = value, x = year, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable, shape = variable), alpha =.4) +
  ylab("Log  rel. abundance") + xlab("Year") + theme_classic()

print(taxITS.OCG_sigsbst_year)

#CLEAR GLOBAL ENVIRONMENT
rm(list = ls())

#METACODER ####
#Read in project data 
#ASV
asvITS<- read.csv("data_csv/asv-table-dada2-ITS-sagebrush.csv",head=T,row.names=1, check.names = F) #5983 obs of 463 variable
asvITS<- asvITS[,order(colnames(asvITS))]

##METADATA
mdITS <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #505 obs of 16 variables.
mdITS <- mdITS[order(row.names(mdITS)),]
mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS)) #361 of 16 variables
colnames(asvITS) == row.names(mdITS) # sanity check true.

##TAXONOMY
tax.ITS <- read.csv("~/Documents/Orchard_Common_Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #5983 obs of 2 variables
row.names(asvITS) == row.names(tax.ITS) #TRUE

#Clean up data for analysis
asvITS[asvITS < 10] <- 0
asvITS <- asvITS[rowSums(asvITS) > 0,]
summary(rowSums(asvITS)) #10
summary(colSums(asvITS)) #0
asvITS <- asvITS[,colSums(asvITS) > 999]
mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS)) # remove trimmed samples from metadata

##Subset data frame to only include samples from the Orchard Common Garden project #
mdITS.OCG <- subset(mdITS, mdITS$Project == "OCG")  # trim metadata to only OCG entries
mdITS.OCG <- subset(mdITS.OCG, mdITS.OCG$Location != "NEG")# remove negative controls
print(mdITS.OCG)
asvITS.t <- t(asvITS) # transpose asv data frame before subsetting from metadata (subsetting didn't play well when comparing rows to columns or vice versa)
asvITS.OCG.t <- subset(asvITS.t, row.names(asvITS.t) %in% row.names(mdITS.OCG))  # trim transposed asv data frame based on trimmed metadata
asvITS.OCG <- t(asvITS.OCG.t)  # transpose asv table back
asvITS.OCG <- as.data.frame(asvITS.OCG)# return asv table to data frame format after transposing

##Trimming of improper samples from OCG subset #
asvITS.OCG[asvITS.OCG < 10] <- 0                  
asvITS.OCG <- asvITS.OCG[rowSums(asvITS.OCG) > 0,] # each observation needs at least 10 seqs

summary(rowSums(asvITS.OCG)) #10
summary(colSums(asvITS.OCG)) #1018

asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 999] # each sample needs at least 1000 seqs

mdITS.OCG <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% colnames(asvITS.OCG)) # remove trimmed samples from metadata

#Create taxmap object for use with Metacoder functions 
##Create data frame matching taxonomy information with ASV sequences ###
taxITS.OCG <- subset(tax.ITS, row.names(tax.ITS) %in% row.names(asvITS.OCG))# remove entries in taxonomy not present in trimmed ASV data
row.names(asvITS.OCG) == row.names(taxITS.OCG)# sanity check:TRUE
totMC_OCG <- asvITS.OCG # new data frame to hold combined asv + taxonomy information
totMC_OCG$Taxon <- taxITS.OCG$Taxon # append taxonomy information to ASV data frame
row.names(taxITS.OCG) == row.names(totMC_OCG)# sanity check for ASVs: TRUE
taxITS.OCG$Taxon == totMC_OCG$Taxon # sanity check for taxonomy: TRUE

##Create taxmap object ####
obj <- parse_tax_data(totMC_OCG,
                      class_cols = "Taxon",                   # name of column that contains input taxon data
                      class_sep = ";",                        # character that separates taxon data
                      class_regex = "^(.+)__(.+)$",           # regex to identify taxon entries
                      class_key = c(tax_rank = "info",        # this is the key that labels each column pulled from the regex, since we had two sections for each identifier we need two columns
                                    tax_name = "taxon_name"))

print(obj) #ERROR                   
print(obj$data$tax_data)
print(obj$data$class_data)
obj$data$class_data <- NULL     # class_data is repetitive/unnecessary
names(obj$data) <- "asv_counts" # rename "data" to something more relevant
print(obj)                      

heat_tree(obj, # very basic heat tree showing the overall composition of Orchard Common Garden samples
          node_label = taxon_names, 
          node_size = n_obs,        
          node_color = n_obs)

#DIFFERENTIAL HEAT TREES FOR OCG
##Adjust taxmap so we can investigate taxa instead of ASVs #
#Now we need to calculate abundances based on taxon not ASV
obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "asv_counts",
                                                  cols = row.names(mdITS.OCG))
print(obj$data$tax_abund)
##**Notice** ###
# Metacoder did not play well with separate naming conventions. As such, all heatmaps are each written into the same places: "diff_table" and "diff_heattree_color"
# To look at different maps, run code starting from the beginning of each section: "Heatmap matrix by _____"

##HEATMAP MATRIX BY SUBSPECIES
obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund",
                                                 cols = row.names(mdITS.OCG), # What columns of sample data to use (sample ID is stored in row names of metadata)
                                                 groups = mdITS.OCG$Subspecies) # What category each sample is assigned to from metadata

print(obj$data$diff_table)
obj <- mutate_obs(obj, "diff_table",                                                                
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

#look at the p-values and see if there is any significance
range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
# [1] 0.000186618 0.99738256
# lower range is significant

## Focus only on significant taxa
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

set.seed(1)
diff_heattree_color <- metacoder::heat_tree_matrix(obj, data = "diff_table",
                                                   node_size = n_obs, 
                                                   node_label = taxon_names,
                                                   node_color = log2_median_ratio,
                                                   node_color_range = diverging_palette(),
                                                   node_color_trans = "linear", 
                                                   node_color_interval = c(-3, 3), 
                                                   edge_color_interval = c(-3, 3), 
                                                   node_size_axis_label = "Number of ASVs",
                                                   node_color_axis_label = "Log2 ratio median proportions",
                                                   layout = "davidson-harel", 
                                                   initial_layout = "reingold-tilford")

## This plot takes a while to load
print(diff_heattree_color) ## Show taxonomic heat tree 



##HEATMAP MATRIX BY PLOIDY
obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund",
                                                 cols = row.names(mdITS.OCG), # What columns of sample data to use (sample ID is stored in row names of metadata)
                                                 groups = mdITS.OCG$Ploidy) # What category each sample is assigned to from metadata

print(obj$data$diff_table)
obj <- mutate_obs(obj, "diff_table",                                                                
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

#lets look at the p-values and see if there is any significance
range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
# [1] 0.5037569 0.9958621
# not significant

## Focus only on significant taxa
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0


set.seed(1)
diff_heattree_color <- metacoder::heat_tree_matrix(obj, data = "diff_table",
                                                   node_size = n_obs, 
                                                   node_label = taxon_names,
                                                   node_color = log2_median_ratio,
                                                   node_color_range = diverging_palette(),
                                                   node_color_trans = "linear", 
                                                   node_color_interval = c(-3, 3), 
                                                   edge_color_interval = c(-3, 3), 
                                                   node_size_axis_label = "Number of ASVs",
                                                   node_color_axis_label = "Log2 ratio median proportions",
                                                   layout = "davidson-harel", 
                                                   initial_layout = "reingold-tilford")

## This plot takes a while to load
print(diff_heattree_color) ## Show taxonomic heat tree 

#HEATMAP MATRIX BY SUBSPECIES PLOIDY
obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund",
                                                 cols = row.names(mdITS.OCG), # What columns of sample data to use (sample ID is stored in row names of metadata)
                                                 groups = mdITS.OCG$Subsp_ploidy) # What category each sample is assigned to from metadata

print(obj$data$diff_table)
obj <- mutate_obs(obj, "diff_table",                                                               
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

#lets look at the p-values and see if there is any significance
range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
# [1] 0.04266539 1.00000000
# the lower range is significant

## Focus only on significant taxa
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

print(obj$data$diff_table)

set.seed(1)
diff_heattree_color <- metacoder::heat_tree_matrix(obj, data = "diff_table",
                                                   node_size = n_obs, 
                                                   node_label = taxon_names,
                                                   node_color = log2_median_ratio,
                                                   node_color_range = diverging_palette(),
                                                   node_color_trans = "linear", 
                                                   node_color_interval = c(-3, 3), 
                                                   edge_color_interval = c(-3, 3), 
                                                   node_size_axis_label = "Number of ASVs",
                                                   node_color_axis_label = "Log2 ratio median proportions",
                                                   layout = "davidson-harel", 
                                                   initial_layout = "reingold-tilford")

## This plot takes a while to load
print(diff_heattree_color) ## Show taxonomic heat tree

##HEATMAP BY YEAR
obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund",
                                                 cols = row.names(mdITS.OCG), # What columns of sample data to use (sample ID is stored in row names of metadata)
                                                 groups = mdITS.OCG$Year) # What category each sample is assigned to from metadata

print(obj$data$diff_table)
obj <- mutate_obs(obj, "diff_table",                                                                
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

#lets look at the p-values and see if there is any significance
range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
# [1] 2.288654e-06 9.026667e-01
# the lower range is significant

## Focus only on significant taxa
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0


set.seed(1)
diff_heattree_color <- metacoder::heat_tree(obj,  # heat tree code slightly different here on account of only having two variables, use heat_tree, exclude data argument
                                            node_size = n_obs, 
                                            node_label = taxon_names,
                                            node_color = log2_median_ratio,
                                            node_color_range = diverging_palette(),
                                            node_color_trans = "linear", 
                                            node_color_interval = c(-3, 3), 
                                            edge_color_interval = c(-3, 3), 
                                            node_size_axis_label = "Number of ASVs",
                                            node_color_axis_label = "Log2 ratio median proportions",
                                            layout = "davidson-harel", 
                                            initial_layout = "reingold-tilford")

## This plot takes a while to load
print(diff_heattree_color) ## Show taxonomic heat tree
#after cross-referencing with ANCOM results, this heat map shows changes from 2012 -> 2021; i.e. positive ratio indicates more sequences in 2021

#Thermal data testing#### 

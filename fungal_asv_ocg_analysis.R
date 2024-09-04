# Orchard Common Garden - Microbial & Chemistry data cleaning, analysis, and visualization
# Set working directory and load necessary packages####
setwd("/Users/ellehorwath/Documents/Orchard_Common_Garden/commongarden")
# List of CRAN packages
cran_packages <- c("dplyr", "effects", "exactRankTests", "factoextra", 
                   "ggcorrplot", "ggfortify", "ggplot2", "glmm", "gridExtra", 
                   "iNEXT", "lme4", "MASS", "metacoder", "nlme", 
                   "picante", "raster", "readr", "reshape2", 
                   "tidyr", "tidyverse", "vegan")

# List of Bioconductor packages
bioc_packages <- c("phyloseq", "qiime2R")

# List of GitHub packages
github_packages <- c("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# Function to install and load CRAN packages
install_and_load_cran <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Function to install and load Bioconductor packages
install_and_load_bioc <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Function to install and load GitHub packages
install_and_load_github <- function(repo) {
  pkg <- basename(repo)
  if (!require(pkg, character.only = TRUE)) {
    devtools::install_github(repo)
    library(pkg, character.only = TRUE)
  }
}

# Check for and install BiocManager if needed
if (!require("BiocManager", character.only = TRUE)) {
  install.packages("BiocManager")
  library(BiocManager, character.only = TRUE)
}

# Check for and install devtools if needed
if (!require("devtools", character.only = TRUE)) {
  install.packages("devtools")
  library(devtools, character.only = TRUE)
}

# Install and load all CRAN packages
sapply(cran_packages, install_and_load_cran)

# Install and load all Bioconductor packages
sapply(bioc_packages, install_and_load_bioc)

# Install and load all GitHub packages
sapply(github_packages, install_and_load_github)

# ASV analysis: includes alpha diversity analysis with glm; beta diversity analysis with NMDS plots, PERMANOVA tests, and pairwise adonis; barchart plots at asv level for just the common garden; betadispersion
# Read in cleaned data : START HERE ####
##METADATA
#Data includes just the orchard common garden plants. The data will include plant ID number, location of origin, year of sampling (2012 or 2021), subspecies, ploidy, and subspecies ploidy, status (dead or alive) in 2020.
mdITS.OCG <- read.csv("data_csv/metadataITS_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #206 of 21 variables
#make variables factor to plot
mdITS.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant")] <- lapply(mdITS.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant")], as.factor)

str(mdITS.OCG)

## ASV DATA
#Data has not been filtered and is not yet clean to include just observations with a certain number of sequences and a certain number of sequences per sample. The data is ordered alphabetically.
asvITS.OCG <- read.csv("data_csv/asvITS.OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #206 obs of 2377 variables
summary(rowSums(asvITS.OCG)) #0
summary(colSums(asvITS.OCG)) #1

### Remove duplicates from ASV
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1')
asvITS.OCG <- asvITS.OCG[!rownames(asvITS.OCG) %in% rows_to_remove, ] #199

## Remove negative control
asvITS.OCG <- asvITS.OCG[!(row.names(asvITS.OCG) == "NEG_8-28-21"),]

## Remove MTW.3.7.R_2012
asvITS.OCG <- asvITS.OCG[!(row.names(asvITS.OCG) == "MTW.3.7.R_2012"),] #197

asvITS.OCG[asvITS.OCG < 10] <- 0 # each observation needs at least 10 seqs.
asvITS.OCG <- asvITS.OCG[rowSums(asvITS.OCG) > 0,] #the values that are greater than zero
summary(rowSums(asvITS.OCG)) #19
summary(colSums(asvITS.OCG)) #0

asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 999] # each sample needs at least 100 seqs. #166 of 228
asvITS.OCG <- asvITS.OCG[rowSums(asvITS.OCG) > 0,] #the values that are greater than zero 

summary(colSums(asvITS.OCG)) #401
summary(rowSums(asvITS.OCG)) #16

mdITS.OCG <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(asvITS.OCG)) ##166 of 21 var

row.names(asvITS.OCG) == row.names(mdITS.OCG) # sanity check:TRUE

mdITS.OCG <- droplevels(mdITS.OCG)
levels(mdITS.OCG$Year)

##TAXONOMY
#taxonomy table is used to match to amplicon sequence variant table to fungal ID.
tax.ITS <- read.csv("~/Documents/Orchard_Common_Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #taxonomy read in 5983 obs of 2 variables

# Alpha diversity asv level####
## Rarefying
asvITS.OCG.r <- rrarefy(asvITS.OCG, 16) ## rarefy: Warning message
asvITS.OCG.shannon <- diversity(asvITS.OCG.r)
asvITS.OCG.ef <- exp(asvITS.OCG.shannon)
asvITS.OCG.efr <- round(asvITS.OCG.ef)
mdITS.OCG <- cbind(mdITS.OCG, effective_species = asvITS.OCG.efr)

glm.OCG <- glm(effective_species ~ Subspecies + Year, family = poisson, data=mdITS.OCG)
summary(glm.OCG) 

glm.OCG.gamma <- glm(effective_species ~ Subspecies + Year, family = Gamma, data=mdITS.OCG)
summary(glm.OCG.gamma) #year is significant

plot(allEffects(glm.OCG))

plot(mdITS.OCG$Year,asvITS.OCG.ef)
plot(mdITS.OCG$Subspecies,asvITS.OCG.ef)

ggplot(mdITS.OCG, aes(Year, effective_species))+
  geom_boxplot(aes(group = Year, fill = Year))+
  theme_classic()

ggplot(data = mdITS.OCG, mapping = aes(x = Subspecies, y = effective_species, fill = Subspecies)) +
  geom_boxplot() +
  theme_classic()

#Beta diversity - NMDS plots#### 
set.seed(1)
asvITS.OCG.nmds <- metaMDS(asvITS.OCG.r, trymax=500) ### Solution reached
# save(asvITS.OCG.nmds, file = "nmds/asvITS_OCG_nmds.rda") #save the nmds so you won't need to run it again
load("nmds/asvITS_OCG_nmds.rda") #load it to use in code anytime after the initial run

ordiplot(asvITS.OCG.nmds, type = "t",display = "sites",cex = .6)
rownames(asvITS.OCG.nmds$points) == rownames(mdITS.OCG)

#SUBSPECIES
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="Sagebrush fungal community by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.OCG$Subspecies],
     pch=c(19,17)[mdITS.OCG$Year])
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
asvITS.OCG.subsp.pw #T vs V= 0.003, T vs W= 0.039, and W vs V= 0.540.

#PLOIDY
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by ploidy", 
     col= c("red","blue")[mdITS.OCG$Ploidy],
     pch=c(19,17)[mdITS.OCG$Year])
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
     main="Fungal community by subspecies, ploidy, & year", 
     col= c("pink","brown","darkgreen",'tan','lightblue')[mdITS.OCG$Subsp_ploidy],
     pch=c(17,19)[mdITS.OCG$Year])
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
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

### PERMANOVA and adonis for subspecies ploidy ##
asvITS.OCG.subsp_ploi <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Subsp_ploidy) 
asvITS.OCG.subsp_ploi #subspecies ploidy is significant 0.002

#pairwiseadonis
asvITS.OCG.subsp.pw <- pairwise.adonis(asvITS.OCG.r, mdITS.OCG$Subsp_ploidy)
asvITS.OCG.subsp.pw #p-adjusted:T_4n vs V_2n = 0.04, T_4n vs V_4n = 0.03, T_4n vs W_4n = 0.05

#YEAR
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush fungal community by year", 
     col= c("maroon","cyan")[mdITS.OCG$Year],
     pch=16)
legend("bottomleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=16,
       cex=0.8,
       bty = "n")
# ordiarrows(asvITS.OCG.nmds,mdITS.OCG$Plant)
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Year, show.groups = "2012", col = "maroon")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Year, show.groups = "2021", col = "cyan")

### PERMANOVAs for year ##
asvITS.OCG.subsp_loc_yr <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Subspecies + mdITS.OCG$Year + mdITS.OCG$Location, by = "margin") 
asvITS.OCG.subsp_loc_yr #Year is significant

asvITS.OCG.subsp_yr <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Subspecies*mdITS.OCG$Year) 
asvITS.OCG.subsp_yr #year and subspecies are sig

asvITS.OCG.yr <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Year) 
asvITS.OCG.yr #year is significant'

# k- means clustering #####
#read in lcms with cluster data saved
md_LCMS_cluster <- read.csv("data_csv/md.OCG.LCMS_cluster.csv", row.names = 1)
md_LCMS_cluster <- md_LCMS_cluster[order(row.names(md_LCMS_cluster)),]
md_LCMS_cluster <- subset(md_LCMS_cluster, row.names(md_LCMS_cluster) %in% row.names(asvITS.OCG))
asvITS.OCG_lcms <- subset(asvITS.OCG, row.names(asvITS.OCG) %in% row.names(md_LCMS_cluster))
row.names(asvITS.OCG_lcms) == row.names(md_LCMS_cluster) # sanity check:TRUE

# rarefy the asv data
asvITS.OCG_lcms <- asvITS.OCG_lcms[,colSums(asvITS.OCG_lcms) > 999]
summary(colSums(asvITS.OCG_lcms)) #1002
summary(rowSums(asvITS.OCG_lcms)) #17
asvITS.OCG_lcms.r <- rrarefy(asvITS.OCG_lcms, 17)

# model
k_means_asv_fit <- adonis2(asvITS.OCG_lcms.r ~ md.OCG.LCMS.3.asv$cluster_assignments, by = "margin")
k_means_asv_fit #subspecies ploidy is significant 0.001

asv_subsp_LCMS <- adonis2(asvITS.OCG_lcms.r ~ md.OCG.LCMS.3.asv$Subsp_ploidy, by = "margin")

# read in GC clustering data
md_gc_2012_cluster <- read.csv("data_csv/md.OCG.GC.2012_cluster.csv")
md_gc_2021_cluster <- read.csv("data_csv/md.OCG.GC.2021_cluster.csv")
md_gc_cluster <- merge(md_gc_2012_cluster, md_gc_2021_cluster, all = TRUE)
rownames(md_gc_cluster) <- md_gc_cluster[, 1]
md_gc_cluster <- md_gc_cluster[, -1]

md_gc_cluster <- md_gc_cluster[order(row.names(md_gc_cluster)),]
md_gc_cluster <- subset(md_gc_cluster, row.names(md_gc_cluster) %in% row.names(asvITS.OCG)) ##166 of 21 var
asvITS.OCG_gc <- subset(asvITS.OCG, row.names(asvITS.OCG) %in% row.names(md_gc_cluster)) ##166 of 21 var
row.names(asvITS.OCG_gc) == row.names(md_gc_cluster) # sanity check:TRUE

# rarefy the asv data
asvITS.OCG_gc <- asvITS.OCG_gc[,colSums(asvITS.OCG_gc) > 999]
summary(colSums(asvITS.OCG_gc)) #1002
summary(rowSums(asvITS.OCG_gc)) #16
asvITS.OCG_gc.r <- rrarefy(asvITS.OCG_gc, 16)

# model
k_means_asv.gc_fit <- adonis2(asvITS.OCG_gc.r ~ md_gc_cluster$cluster_assignments, by = "margin")
k_means_asv.gc_fit #subspecies ploidy is significant 0.001

asv_subsp_GC <- adonis2(asvITS.OCG_gc.r ~ md_gc_cluster$Subspecies, by = "margin")

## 2012 asv NMDS ####
#NMDS
asvITS.2012 <- subset(asvITS.OCG, mdITS.OCG$Year=="2012") #105 of 2135 variables
asvITS.2012 <- asvITS.2012[,colSums(asvITS.2012) > 0]
summary(rowSums(asvITS.2012)) #16
summary(colSums(asvITS.2012)) #305

mdITS.2012 <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(asvITS.2012)) #112 of 21
str(mdITS.2012)

#rarefy
asvITS.2012.r <- rrarefy(asvITS.2012,16) ## rarefy. warning message

set.seed(12)
asvITS.2012.nmds <- metaMDS(asvITS.2012.r, trymax=500) ###solution reached! warning message
save(asvITS.2012.nmds, file = "nmds/asvITS.2012.nmds.rda")
load("nmds/asvITS.2012.nmds.rda")

ordiplot(asvITS.2012.nmds, type = "t",display = "sites",cex = .6)
rownames(asvITS.2012.nmds$points) == rownames(mdITS.2012)
#SUBSPECIES
plot(asvITS.2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2012 fungal community by subspecies", 
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.2012$Subspecies],
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
asvITS.2012.ploidy <- adonis2(asvITS.2012.r ~ mdITS.2012$Ploidy,by="margin") # Bray-Curtis is the default metric
asvITS.2012.ploidy #ploidy is not significant

#SUBSP PLOIDY
plot(asvITS.2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2012 fungal community by subspecies and ploidy", 
     col= c("pink","brown","darkgreen",'tan','lightblue')[mdITS.2012$Subsp_ploidy],
     pch=c(19))
legend("topright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
# ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "T_2n", col = "red")
# ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "T_4n", col = "orange")
# ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "V_2n", col = "green")
# ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
# ordispider(asvITS.2012.nmds,groups = mdITS.2012$Subsp_ploidy, show.groups = "W_4n", col = "purple")

### PERMANOVA and adonis for subspecies ploidy ##
asvITS.2012.subsp_ploi <- adonis2(asvITS.2012.r ~ mdITS.2012$Subsp_ploidy) 
asvITS.2012.subsp_ploi #subspecies ploidy is not sig

#pairwiseadonis
asvITS.2012.subsp.pw <- pairwise.adonis(asvITS.2012.r, mdITS.2012$Subsp_ploidy)
asvITS.2012.subsp.pw #none sig

## 2021 asv NMDS ####
asvITS.2021 <- subset(asvITS.OCG, mdITS.OCG$Year=="2021") #49 of 769 var 
asvITS.2021 <- asvITS.2021[,colSums(asvITS.2021) > 0]

summary(rowSums(asvITS.2021)) #19
summary(colSums(asvITS.2021)) #73

mdITS.2021 <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(asvITS.2021)) #54 obs of 22

#rarefying
asvITS.2021.r <- rrarefy(asvITS.2021,19) ## rarefy. warning message

set.seed(78)
# asvITS.2021.nmds <- metaMDS(asvITS.2021.r, trymax=500) ###solution reached!
# save(asvITS.2021.nmds, file = "nmds/asvITS.2021.nmds.rda")
load("nmds/asvITS.2021.nmds.rda")

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
asvITS.2021.ad <- adonis2(asvITS.2021.r ~ mdITS.2021$Subspecies) # Bray-Curtis is the default metric
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
asvITS.2021.ploidy <- adonis2(asvITS.2021.r ~ mdITS.2021$Ploidy,by="margin") # Bray-Curtis is the default metric
asvITS.2021.ploidy #ploidy is significant 0.024

#SUBSPECIES PLOIDY
plot(asvITS.2021.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush 2021 fungal community by subspecies and ploidy", 
     col= c("pink","brown","darkgreen",'tan',"lightblue")[mdITS.2021$Subsp_ploidy],
     pch=19)
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan',"lightblue"),
       pch=19,
       cex=0.8,
       bty = "n")

#there is no V_2n in 2021

### PERMANOVA and adonis for subspecies ploidy ##
asvITS.2021.subsp_ploi <- adonis2(asvITS.2021.r ~ mdITS.2021$Subsp_ploidy) 
asvITS.2021.subsp_ploi #subspecies ploidy is not sig

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

##### CLASS LEVEL OCG
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

# Betadispersion ####
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

#ANCOM: analysis of composition of microbiomes. Differential abundance analysis for common garden####
asvITS.OCG.t <- t(asvITS.OCG)
asvITS.OCG.t <- asvITS.OCG.t[,order(colnames(asvITS.OCG.t))] # order samples alphabetically

mdITS.OCG <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(asvITS.OCG)) 
row.names(asvITS.OCG) == row.names(mdITS.OCG)# sanity check :TRUE

taxITS.OCG <- subset(tax.ITS, row.names(tax.ITS) %in% rownames(asvITS.OCG.t))
row.names(taxITS.OCG) == rownames(asvITS.OCG.t)  # sanity check: TRUE

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
  
  ######################## Function #################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  ######################## Function #################
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
taxITS.OCG_t <- t(taxITS.OCG) #transpose ASV and taxonomy data frames so sample IDs become row names
asvITS.OCG_sbst <- data.frame("Sample.ID" = row.names(asvITS.OCG), asvITS.OCG, check.names = F)#create new column in each of ASV, metadata, and taxonomy data frames labelled "Sample.ID" containing sample IDs
mdITS.OCG_sbst <- data.frame("Sample.ID" = row.names(mdITS.OCG), mdITS.OCG)
taxITS.OCG_sbst <- data.frame("Sample.ID" = row.names(taxITS.OCG_t), taxITS.OCG_t, check.names = F)

row.names(asvITS.OCG_sbst) == row.names(mdITS.OCG_sbst) #TRUE

#SUBSPECIES 
ANCOM_subspecies <- ANCOM.main(asvITS.OCG_sbst,mdITS.OCG_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)
#Create objects of significant ASVs
sigASVs_subspecies <- subset(ANCOM_subspecies$W.taxa, ANCOM_subspecies$W.taxa$W_stat > 0)[,1]
sigASVs_subspecies <- as.data.frame(sigASVs_subspecies) #convert significant ASV return to data frame
row.names(sigASVs_subspecies) <- sigASVs_subspecies[1:8,1] #change row names to ASV sequence

#By default ANCOM sorts ASVs in order of decreasing abundance. This section instead aims to sort by decreasing abundance, which can be cross-referenced with the .csv generated after the break
#If you want to look at the default ANCOM output, the full code will be included at the end of this section

sigASVs_subspecies[,1] <- c(1:8) #append numerical ranking of significance in additional column
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
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  ylab("Log  rel. abundance") + xlab("Subspecies") + theme_classic()

print(taxITS.OCG_sigsbst_subspecies_default)
print(taxITS.OCG_sigsbst_subspecies_default)

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
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  ylab("Log  rel. abundance") + xlab("Subsp_ploidy") + theme_classic()

print(taxITS.OCG_sigsbst_subsp_ploidy)

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
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  ylab("Log  rel. abundance") + xlab("Year") + theme_classic()

print(taxITS.OCG_sigsbst_year)

## ANCOM BC 2 ####
# Create phyloseq object
asvITS.OCG_t <- t(asvITS.OCG)
asvITS.OCG_t <- as.data.frame(asvITS.OCG_t)
taxITS.OCG <- subset(taxITS.OCG, row.names(taxITS.OCG) %in% rownames(asvITS.OCG_t))
taxITS.OCG <- taxITS.OCG[order(row.names(taxITS.OCG)),] # order samples alphabetically
row.names(taxITS.OCG) == row.names(asvITS.OCG_t) #TRUE
colnames(asvITS.OCG_t) == row.names(mdITS.OCG) #TRUE
otu_asv_table <- otu_table(as.matrix(asvITS.OCG_t), taxa_are_rows = TRUE)
tax_ITS_table <- tax_table(as.matrix(taxITS.OCG))
md_ITS_phy <- sample_data(mdITS.OCG) #make metadata into sample data to create phyloseq
physeqITS <- merge_phyloseq(phyloseq(otu_asv_table),md_ITS_phy,tax_ITS_table) #create phyloseq object for ANCOM BC
physeqITS #86 taxa and 164 samples

#MAKE TSE
tse.ITS = mia::makeTreeSummarizedExperimentFromPhyloseq(physeqITS)
tse.ITS$Subspecies <- factor(tse.ITS$Subspecies, levels=c("T", "W", "V"))
tse.ITS$Ploidy <- factor(tse.ITS$Ploidy, levels=c("2n", "4n"))
tse.ITS$Subsp_ploidy <- factor(tse.ITS$Subsp_ploidy, levels=c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))
tse.ITS$Year <- factor(tse.ITS$Year, levels =c("2012", "2021"))

#Model for Subspecies
result_ITS_sub <- ancombc2( 
  data = tse.ITS, assay_name = "counts", tax_level = NULL,
  fix_formula = "Subspecies + Ploidy + Year",
  p_adj_method = "fdr", pseudo_sens = TRUE,
  group = "Subspecies", #change variable of interest
  alpha = 0.05, verbose = TRUE,
  global = TRUE, prv_cut = 0.02
)

summary(tse.ITS$Year)

ancomITSres_df.sub <- result_ITS_sub$res

#Model for year
result_ITS_yr <- ancombc2( 
  data = tse.ITS, assay_name = "counts", tax_level = NULL,
  fix_formula = "Subsp_ploidy + Year",
  p_adj_method = "fdr", pseudo_sens = TRUE,
  group = "Year", #change variable of interest
  alpha = 0.05, verbose = TRUE,
  global = TRUE, prv_cut = 0.1
)

ancomITSres_df.yr <- result_ITS_yr$res

#iter_control = list(tol = 1e-5, max_iter = 200, 
                    #verbose = FALSE)

#METACODER ####
#Create taxmap object for use with Metacoder functions 
##Create data frame matching taxonomy information with ASV sequences ###
row.names(asvITS.OCG) == row.names(taxITS.OCG)# sanity check:TRUE
totMC_OCG <- asvITS.OCG # new data frame to hold combined asv + taxonomy information
totMC_OCG$Taxon <- taxITS.OCG$Taxon # append taxonomy information to ASV data frame
row.names(taxITS.OCG) == row.names(totMC_OCG)# sanity check for ASVs: TRUE
taxITS.OCG$Taxon == totMC_OCG$Taxon # sanity check for taxonomy: TRUE

## Create taxmap object ####
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

###### Clear Global Environment ####
rm(list = ls())

# Procrustes analysis here ####
#Read in cleaned data ####
#METADATA
md.OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #metadata read in #154 of 16 variables
head(md.OCG)

#ASV
asvITS.OCG <- read.csv("data_csv/asvITS.OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #asv table read in 154 obs of 2135 variables

summary(rowSums(asvITS.OCG)) #507
summary(colSums(asvITS.OCG)) #2

row.names(asvITS.OCG) == row.names(md.OCG) # TRUE

#LCMS 3UL
#LCMS 3uL read in and subset with asv
OCG_LCMS_3uL <- read.csv("data_csv/OCG_LCMS_3uL_cleaned.csv", row.names = 1) #110 obs of 309 var
OCG_LCMS_3uL <- OCG_LCMS_3uL[order(row.names(OCG_LCMS_3uL)),] # order samples alphabetically

OCG_LCMS_3uLp <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.OCG)) #59 of 309 variables

asvITS.OCG.LC3 <- subset(asvITS.OCG, row.names(asvITS.OCG) %in% row.names(OCG_LCMS_3uLp)) #59 of 2135 variables

#GC
OCG_GC <- read.csv("data_csv/OCG_GC_full_clean.csv", row.names = 1) #217
OCG_GCp <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(md.OCG)) #140
OCG_GCp <- OCG_GCp[order(row.names(OCG_GCp)),]
asvITS.OCG.GC <- subset(asvITS.OCG, row.names(asvITS.OCG) %in% row.names(OCG_GCp)) #140 of 2135 var
rownames(OCG_GCp) == rownames(asvITS.OCG.GC) #TRUE

#LCMS to GC subset
OCG_GCp2 <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(OCG_LCMS_3uL)) #107
OCG_GCp2 <- OCG_GCp2[order(row.names(OCG_GCp2)),]
OCG_LCMS_3uLp2 <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(OCG_GC)) #107 var
rownames(OCG_GCp2) == rownames(OCG_LCMS_3uLp2) #TRUE

## Run NMDS (not needed if already done) ####
# Replace NA with 0
OCG_LCMS_3uLp[is.na(OCG_LCMS_3uLp)] <- 0
OCG_LCMS_3uLp2[is.na(OCG_LCMS_3uLp2)] <- 0
OCG_GCp[is.na(OCG_GCp)] <- 0
OCG_GCp2[is.na(OCG_GCp2)] <- 0

#rarefy

#ASV FOR LCMS 3
summary(rowSums(asvITS.OCG.LC3)) #507
summary(colSums(asvITS.OCG.LC3)) #0.00 Need to get rid of ASVs that are 0!
asvITS.OCG.LC3 <- asvITS.OCG.LC3[,colSums(asvITS.OCG.LC3) > 0]
summary(colSums(asvITS.OCG.LC3)) #2
asvITS.LCMS_3.r <- rrarefy(asvITS.OCG.LC3,507) ## rarefy.

#LCMS 3UL
summary(rowSums(OCG_LCMS_3uLp)) #47200283
summary(colSums(OCG_LCMS_3uLp)) #12471
OCG_LCMS_3uLp.r <- rrarefy(round(OCG_LCMS_3uLp),47200283) ## rarefy.

#LCMS FOR ALL GC
summary(rowSums(OCG_LCMS_3uLp2)) #31016795
summary(colSums(OCG_LCMS_3uLp2)) #254831
OCG_LCMS_3uLp2.r <- rrarefy(round(OCG_LCMS_3uLp2),31016795) ## rarefy. 

#ASV FOR ALL GC
summary(rowSums(asvITS.OCG.GC)) #507
summary(colSums(asvITS.OCG.GC)) #0.00
asvITS.OCG.GC <- asvITS.OCG.GC[,colSums(asvITS.OCG.GC) > 0]
summary(colSums(asvITS.OCG.GC)) #2
asvITS.OCG.GC.r <- rrarefy(asvITS.OCG.GC,507) ## rarefy

#ALL GC
summary(rowSums(OCG_GCp)) #179
summary(colSums(OCG_GCp)) #151
OCG_GCp.r <- rrarefy(round(OCG_GCp),179) ## rarefy. 

#ALL GC FOR LCMS
summary(rowSums(OCG_GCp2)) #179
summary(colSums(OCG_GCp2)) #0
OCG_GCp2.r <- rrarefy(round(OCG_GCp2),179) ## rarefy. 

#ASV NMDS TO MATCH LCMS 3
set.seed(87)
asvITS_OCG_LCMS3.nmds <- metaMDS(asvITS.LCMS_3.r, trymax=500) ###solution reached! warning message
save(asvITS_OCG_LCMS3.nmds, file = "nmds/asvITS_OCG_LCMS3.nmds.rda")

#LCMS 3 TO MATCH ASV NMDS
set.seed(56)
OCG_LCMS3_pro.nmds <- metaMDS(OCG_LCMS_3uLp.r, trymax=500) ###solution reached!
save(OCG_LCMS3_pro.nmds, file = "nmds/OCG_LCMS3_pro.nmds.rda")

#LCMS 3 TO MATCH ALL GC NMDS
set.seed(55)
OCG_LCMS3_prop2.nmds <- metaMDS(OCG_LCMS_3uLp2.r, trymax=500) ###solution reached!
save(OCG_LCMS3_prop2.nmds, file = "nmds/OCG_LCMS3_prop2.nmds.rda")

#ASV NMDS TO MATCH ALL GC
set.seed(98)
asvITS_OCG_GCpp.nmds <- metaMDS(asvITS.OCG.GC.r, trymax=500) ###solution reached! warning message
save(asvITS_OCG_GCpp.nmds, file = "nmds/asvITS_OCG_GCpp.nmds.rda")

#ALL GC TO MATCH ASV
set.seed(5)
OCG_GCp_pro.nmds <- metaMDS(OCG_GCp.r, trymax=500) ###solution reached!
save(OCG_GCp_pro.nmds, file = "nmds/OCG_GCp_pro.nmds.rda")

#ALL GC TO MATCH ASV
set.seed(5)
OCG_GCp2_pro.nmds <- metaMDS(OCG_GCp2.r, trymax=500) ###solution reached!
save(OCG_GCp2_pro.nmds, file = "nmds/OCG_GCp2_pro.nmds.rda")

# Load all nmds files ####
load("nmds/asvITS_OCG_LCMS3.nmds.rda") #asv file for LCMS 3uL data = asvITS.OCG
load("nmds/OCG_LCMS3_prop2.nmds.rda") #LCMS 3uL file for GC data = OCG_LCMS_3uLp2
load("nmds/OCG_LCMS3_pro.nmds.rda") #LCMS 3uL file for asv data = OCG_LCMS_3uLp
load("nmds/asvITS_OCG_GCpp.nmds.rda") #asv file for all GC data = asvITS.OCG
load("nmds/OCG_GCp_pro.nmds.rda") #all GC file for asv data = OCG_GCp
load("nmds/OCG_GCp2_pro.nmds.rda") #all GC file for LCMS 3uL data = OCG_GCp2

# Procrustes plots ####
# LCMS 3uL procrustes ####
asv_LCMS_3.pro <- protest(asvITS_OCG_LCMS3.nmds, OCG_LCMS3_pro.nmds, symmetric=T) 
asv_LCMS_3.pro ## Correlation in a symmetric Procrustes rotation: 0.1575, Significance: 0.434
summary(asv_LCMS_3.pro) 
plot(asv_LCMS_3.pro)

## Re-plotting the procrustes 
asv_LCMS_3_prodat <- as.data.frame(asv_LCMS_3.pro$X)
asv_LCMS_3_prodat <- cbind(asv_LCMS_3_prodat,asv_LCMS_3.pro$Yrot)
colnames(asv_LCMS_3_prodat)[colnames(asv_LCMS_3_prodat)=="1"] <- "Xend"
colnames(asv_LCMS_3_prodat)[colnames(asv_LCMS_3_prodat)=="2"] <- "Yend"

ggplot() + 
  geom_segment(data=asv_LCMS_3_prodat, mapping=aes(x=NMDS1, y=NMDS2, xend=Xend, yend=Yend), size=0.8, color="gray") + 
  geom_point(data=asv_LCMS_3_prodat, mapping=aes(x=NMDS1, y=NMDS2), size=2, shape=19, color = "maroon") +
  geom_point(data=asv_LCMS_3_prodat, mapping=aes(x=Xend, y=Yend), size=2, shape=17, color = "lightseagreen") +
  labs(x="Procrustes axis 1", y="Procrustes axis 2",title = "LCMS3 vs asv procrustes plot") +
  theme_classic() 

# All GC procrustes ####
asv_GC.pro <- protest(asvITS_OCG_GCpp.nmds, OCG_GCp_pro.nmds, symmetric=T) 
asv_GC.pro ## Correlation in a symmetric Procrustes rotation: 0.348, Significance: 0.001
summary(asv_GC.pro) 
plot(asv_GC.pro)

## Re-plotting the procrustes 
asv_GC.pro_prodat <- as.data.frame(asv_GC.pro$X)
asv_GC.pro_prodat <- cbind(asv_GC.pro_prodat,asv_GC.pro$Yrot)
colnames(asv_GC.pro_prodat)[colnames(asv_GC.pro_prodat)=="1"] <- "Xend"
colnames(asv_GC.pro_prodat)[colnames(asv_GC.pro_prodat)=="2"] <- "Yend"

ggplot() + 
  geom_segment(data=asv_GC.pro_prodat, mapping=aes(x=NMDS1, y=NMDS2, xend=Xend, yend=Yend), size=0.8, color="gray") + 
  geom_point(data=asv_GC.pro_prodat, mapping=aes(x=NMDS1, y=NMDS2), size=2, shape=19, color = "maroon") +
  geom_point(data=asv_GC.pro_prodat, mapping=aes(x=Xend, y=Yend), size=2, shape=17, color = "lightseagreen") +
  labs(x="Procrustes axis 1", y="Procrustes axis 2",title = "GC vs ASV Procrustes Plot") +
  theme_classic() 

# GC vs LCMS procrustes ####
LCMS_GC.pro <- protest(OCG_LCMS3_prop2.nmds, OCG_GCp2_pro.nmds, symmetric=T) 
LCMS_GC.pro ## Correlation in a symmetric Procrustes rotation: 0.3824, Significance: 0.001
summary(LCMS_GC.pro) 
plot(LCMS_GC.pro)

## Re-plotting the procrustes 
LCMS_GC.pro_prodat <- as.data.frame(LCMS_GC.pro$X)
LCMS_GC.pro_prodat <- cbind(LCMS_GC.pro_prodat,LCMS_GC.pro$Yrot)
colnames(LCMS_GC.pro_prodat)[colnames(LCMS_GC.pro_prodat)=="1"] <- "Xend"
colnames(LCMS_GC.pro_prodat)[colnames(LCMS_GC.pro_prodat)=="2"] <- "Yend"

ggplot() + 
  geom_segment(data=LCMS_GC.pro_prodat, mapping=aes(x=NMDS1, y=NMDS2, xend=Xend, yend=Yend), size=0.6, color="lightgrey") + 
  geom_point(data=LCMS_GC.pro_prodat, mapping=aes(x=NMDS1, y=NMDS2), size=2, shape=19, color = "maroon") +
  geom_point(data=LCMS_GC.pro_prodat, mapping=aes(x=Xend, y=Yend), size=2, shape=17, color = "lightseagreen") +
  labs(x="Procrustes axis 1", y="Procrustes axis 2",title = "GC vs LCMS Procrustes Plot") +
  theme_classic() 

# #LCMS against 1ul and 3ul procrustes to check #### 
# OCG_LCMS_3uLp. <- subset(OCG_LCMS_3uLp, row.names(OCG_LCMS_3uLp) %in% row.names(OCG_LCMS_1uL)) #58
# OCG_LCMS_1uL. <- subset(OCG_LCMS_1uL, row.names(OCG_LCMS_1uL) %in% row.names(OCG_LCMS_3uLp)) #58
# 
# set.seed(72)
# OCG_LCMS3_prot.nmds <- metaMDS(OCG_LCMS_3uLp., trymax=500) ###solution reached!
# 
# 
# set.seed(92)
# OCG_LCMS1_prot.nmds <- metaMDS(OCG_LCMS_1uL., trymax=500) ###solution reached!
# 
# LCMS_3_1.pro <- protest(OCG_LCMS1_prot.nmds ,OCG_LCMS3_prot.nmds, symmetric=T) 
# LCMS_3_1.pro ## Correlation in a symmetric Procrustes rotation: 0.7926, Significance: 0.001
# summary(LCMS_3_1.pro) 
# plot(LCMS_3_1.pro)

#Mantel tests measuring correlation between the distance matrices####
#GC 
OCG_GCp.dist <- vegdist(OCG_GCp)
asvITS.OCG.GC.dist <- vegdist(asvITS.OCG.GC)
GC.asv.mant <- mantel(OCG_GCp.dist,asvITS.OCG.GC.dist, permutations = 9999) 
GC.asv.mant ## r = 0.008694, P = 0.4007

#LCMS3
OCG_LCMS_3uLp.dist <- vegdist(OCG_LCMS_3uLp)
asvITS.OCG.LC3.dist <- vegdist(asvITS.OCG.LC3)
LCMS3.asv.mant <- mantel(OCG_LCMS_3uLp.dist,asvITS.OCG.LC3.dist, permutations = 9999) 
LCMS3.asv.mant ## r = -0.1345, P =0.9879

#Trash ####
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
# asvITS.OCG[asvITS < 10] <- 0
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
# mdITS.OCG <- subset(mdITS, mdITS$Project == "OCG") #146 obs of 16 var   
# mdITS.OCG <- subset(mdITS.OCG, mdITS.OCG$Location != "NEG")  
# asvITS.OCG.s <- subset(asvITS.OCG, row.names(asvITS.OCG) %in% row.names(mdITS.OCG)) 
# asvITS.OCG <- as.data.frame(asvITS.OCG.s) 
# 
# asvITS.OCG[asvITS.OCG < 10] <- 0  # repeat cleaning after trimming
# asvITS.OCG <- asvITS.OCG[rowSums(asvITS.OCG) > 0,] # each observation needs at least 10 seqs
# 
# summary(rowSums(asvITS.OCG)) #16
# summary(colSums(asvITS.OCG)) #1002
# 
# asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 999] # each sample needs at least 1000 seqs

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
# mdITS.OCG <- subset(mdITS, mdITS$Project == "OCG")  # trim metadata to only OCG entries
# mdITS.OCG <- subset(mdITS.OCG, mdITS.OCG$Location != "NEG")# remove negative controls
# print(mdITS.OCG)
# asvITS.t <- t(asvITS) # transpose asv data frame before subsetting from metadata (subsetting didn't play well when comparing rows to columns or vice versa)
# asvITS.OCG.t <- subset(asvITS.t, row.names(asvITS.t) %in% row.names(mdITS.OCG))  # trim transposed asv data frame based on trimmed metadata
# asvITS.OCG <- t(asvITS.OCG.t)  # transpose asv table back
# asvITS.OCG <- as.data.frame(asvITS.OCG)# return asv table to data frame format after transposing
# 
# ##Trimming of improper samples from OCG subset #
# asvITS.OCG[asvITS.OCG < 10] <- 0                  
# asvITS.OCG <- asvITS.OCG[rowSums(asvITS.OCG) > 0,] # each observation needs at least 10 seqs
# 
# summary(rowSums(asvITS.OCG)) #10
# summary(colSums(asvITS.OCG)) #1018
# 
# asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 999] # each sample needs at least 1000 seqs
# 
# mdITS.OCG <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% colnames(asvITS.OCG)) # remove trimmed samples from metadata



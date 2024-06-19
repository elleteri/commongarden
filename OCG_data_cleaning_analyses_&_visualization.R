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
str(mdITS.OCG)

## ASV DATA
#Data has not been filtered and is not yet clean to include just observations with a certain number of sequences and a certain number of sequences per sample. The data is ordered alphabetically.
asvITS.OCG <- read.csv("data_csv/asvITS.OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #206 obs of 2372 variables
summary(rowSums(asvITS.OCG)) #0
summary(colSums(asvITS.OCG)) #1

### Remove duplicates from ASV
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1')
asvITS.OCG <- asvITS.OCG[!rownames(asvITS.OCG) %in% rows_to_remove, ]

## Remove negative control
asvITS.OCG <- asvITS.OCG[!(row.names(asvITS.OCG) == "NEG_8-28-21"),]

## Remove MTW.3.7.R_2012
asvITS.OCG <- asvITS.OCG[!(row.names(asvITS.OCG) == "MTW.3.7.R_2012"),] #197

asvITS.OCG[asvITS.OCG < 4] <- 0 # each observation needs at least 4 seqs.
asvITS.OCG <- asvITS.OCG[rowSums(asvITS.OCG) > 0,] #the values that are greater than zero
summary(rowSums(asvITS.OCG)) #4
summary(colSums(asvITS.OCG)) #0

asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 49] # each sample needs at least 100 seqs. #189 of 842

summary(colSums(asvITS.OCG)) #50
summary(rowSums(asvITS.OCG)) #0

mdITS.OCG <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(asvITS.OCG)) ##189 of 21 var

row.names(asvITS.OCG) == row.names(mdITS.OCG) # sanity check:TRUE

##TAXONOMY
#taxonomy table is used to match to amplicon sequence variant table to fungal ID.
tax.ITS <- read.csv("~/Documents/Orchard_Common_Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #taxonomy read in 5983 obs of 2 variables

# Alpha diversity asv level####
## Rarefying
asvITS.OCG.r <- rrarefy(asvITS.OCG,50) ## rarefy: Warning message
asvITS.OCG.shannon <- diversity(asvITS.OCG.r)
asvITS.OCG.ef <- exp(asvITS.OCG.shannon)
asvITS.OCG.efr <- round(asvITS.OCG.ef)
mdITS.OCG <- cbind(mdITS.OCG, effective_species = asvITS.OCG.efr)

glm.OCG <- glm(effective_species ~ Subspecies + Year, family = poisson, data=mdITS.OCG)
summary(glm.OCG) # year p-value = 0.000460

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
set.seed(41)
#asvITS.OCG.nmds <- metaMDS(asvITS.OCG.r, trymax=500) ### Solution reached! gives me a warning
# save(asvITS.OCG.nmds, file = "nmds/asvITS_OCG_nmds.rda") #save the nmds so you won't need to run it again
load("nmds/asvITS_OCG_nmds.rda") #load it to use in code anytime after the initial run

ordiplot(asvITS.OCG.nmds, type = "t",display = "sites",cex = .6)
rownames(asvITS.OCG.nmds$points) == rownames(mdITS.OCG)

#make variables factor to plot
mdITS.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant")] <- lapply(mdITS.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant")], as.factor)
str(mdITS.OCG)

#SUBSPECIES
plot(asvITS.OCG.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="Sagebrush fungal community by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.OCG$Subspecies],
     pch=c(16,17)[mdITS.OCG$Year])
legend("bottomleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=16,
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
legend("bottomleft", 
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
     pch=c(17,16)[mdITS.OCG$Year])
legend("bottomleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=16,
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
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Year, show.groups = "2012", col = "maroon")
ordispider(asvITS.OCG.nmds,groups = mdITS.OCG$Year, show.groups = "2021", col = "cyan")

### PERMANOVAs for year ##
asvITS.OCG.subsp_loc_yr <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Subspecies + mdITS.OCG$Year + mdITS.OCG$Location, by = "margin") 
asvITS.OCG.subsp_loc_yr #Year is significant

asvITS.OCG.subsp_yr <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Subspecies*mdITS.OCG$Year) 
asvITS.OCG.subsp_yr #year and subspecies are sig

asvITS.OCG.yr <- adonis2(asvITS.OCG.r ~ mdITS.OCG$Year) 
asvITS.OCG.yr #year is significant

## 2012 asv NMDS ####
#NMDS
asvITS.2012 <- subset(asvITS.OCG, mdITS.OCG$Year=="2012") #105 of 2135 variables
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

summary(rowSums(asvITS.2021)) #648
summary(colSums(asvITS.2021)) #2

mdITS.2021 <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(asvITS.2021)) #49 obs of 16

#rarefying
asvITS.2021.r <- rrarefy(asvITS.2021,648) ## rarefy. warning message

set.seed(78)
#asvITS.2021.nmds <- metaMDS(asvITS.2021.r, trymax=500) ###solution reached!
#save(asvITS.2021.nmds, file = "nmds/asvITS.2021.nmds.rda")
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
     col= c("pink","brown","darkgreen",'tan')[mdITS.2021$Subsp_ploidy],
     pch=19)
legend("topleft", 
       legend=c("T_2n","T_4n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan'),
       pch=19,
       cex=0.8,
       bty = "n")

#there is no V_2n in 2021

### PERMANOVA and adonis for subspecies ploidy ##
asvITS.2021.subsp_ploi <- adonis2(asvITS.2021.r ~ mdITS.2021$Subsp_ploidy) 
asvITS.2021.subsp_ploi #subspecies ploidy is sig

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

###### Clear Global Environment ####
rm(list = ls())

#ANCOM: analysis of composition of microbiomes. Differential abundance analysis for common garden####
#Read in project data 
#ASV
asvITS <- read.csv("data_csv/asv-table-dada2-ITS-sagebrush.csv",head=T,row.names=1, check.names = F) #5983 obs of 463 var
asvITS <- asvITS[,order(colnames(asvITS))] # order samples alphabetically

#METADATA
mdITS <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) # 505 obs of 21 var
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
sort(colSums(asvITS)) 


asvITS <- asvITS[,colSums(asvITS) > 499] # each sample needs at least 1000 seqs 

mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS)) # remove trimmed samples from metadata #361 obs o f 16 var

mdITS.OCG <- subset(mdITS, mdITS$Project == "OCG") #146 obs of 16 var   
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
row.names(taxITS.OCG) == row.names(asvITS.OCG) #TRUE
colnames(asvITS.OCG) == row.names(mdITS.OCG) #TRUE
otu_asv_table <- otu_table(as.matrix(asvITS.OCG), taxa_are_rows = TRUE)
tax_ITS_table <- tax_table(as.matrix(taxITS.OCG))
md_ITS_phy <- sample_data(mdITS.OCG) #make metadata into sample data to create phyloseq
physeqITS <- merge_phyloseq(phyloseq(otu_asv_table),md_ITS_phy,tax_ITS_table) #create phyloseq object for ANCOM BC
physeqITS #1719 taxa and 146 samples

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

#It will not run with subspecies ploidy in the model structure

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

###### Clear Global Environment ####
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

# Chemistry analysis: GC and LCMS data from both 2012 and 2021. This analysis has PCAs, NMDS, PERMANOVAS, pairwise adonis, alpha diversity glms, binary jaccard plots, PCoA, and procrustes####

# Read data in, DONT RUN#### 
##METADATA READ IN
md <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #505 obs of 16 variables.
md <- md[order(row.names(md)),] #alphabetical

### Subsetting to just the plant in the common garden (OCG)#
md.OCG <- subset(md, md$Project=="OCG")#246 of 16 variables

### Remove duplicates
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1', 'ORT.1.5_2012')
md.OCG <- md.OCG[!rownames(md.OCG) %in% rows_to_remove, ]

## Remove negative control
md.OCG <- md.OCG[!(row.names(md.OCG) == "NEG_8-28-21"),]
md.OCG <- md.OCG[!(row.names(md.OCG) == "NEG_10-2-20"),]

## Remove MTW.3.7.R_2012
md.OCG <- md.OCG[!(row.names(md.OCG) == "MTW.3.7.R_2012"),] 

#subset the md to only have observations from 2012 to avoid duplicates
md.OCG.2012 <- subset(md.OCG, md.OCG$Year=="2012") #159
str(md.OCG.2012)

#subset the md to only have observations from 2021 to avoid duplicates
md.OCG.2021 <- subset(md.OCG, md.OCG$Year=="2021") #76
str(md.OCG.2021)

## 2012 GC RAW DATA READ IN
#read in 2012 GC chemistry data
OCG_GC_2012 <- read.csv("data_csv/OCG_2012_GC.csv", head=T,check.names = F,stringsAsFactors = T, skip = 1) #183 obs of 223 var

## 2012 GC CLEANING 
names(OCG_GC_2012)[1] <- "Plant_ID"

sum(duplicated(OCG_GC_2012$Plant_ID)) #4 duplicates (control, cocktails, and blanks)

# Remove rows where "empty", "Empty", "Ct, "CT", "M", "w"  occur
OCG_GC_2012 <- OCG_GC_2012[!(apply(OCG_GC_2012, 1, function(row) any(grepl("^empty|^Empty|^Ct|^CT |^M|w", row)))), ] #167 of 220 variables

sum(duplicated(OCG_GC_2012$Plant_ID)) # 0 duplicates
#head(OCG_GC_2021)

# Rename duplicate columns
names(OCG_GC_2012) <- make.unique(names(OCG_GC_2012))

# Add 0.537 to non-NA values in RT columns
# OCG_GC_2012 <- OCG_GC_2012 %>%
#   mutate(across(starts_with("RT"), ~ ifelse(!is.na(.x), .x + 0.537, .x)))

#subset to only columns that contain "Peak Area" and "Plant ID". This removes "RT". comment out if using RT
OCG_GC_2012 <- OCG_GC_2012 [, grepl("Peak.Area|Plant_ID", colnames(OCG_GC_2012))] #167 obs of 149 variables. 

#subset to remove columns that contain Peak Area Percent and just keep Peak Area.
OCG_GC_2012 <- OCG_GC_2012 [, !grepl("Peak.Area.Percent", colnames(OCG_GC_2012))] #167 obs of 75 variables

#subset to have just the columns that contain "Peak Area" 1:74   
peak_area_cols <- grep("Peak.Area", colnames(OCG_GC_2012)) 

#the new column names will replace the repeating "Peak.Area" names to be "C001" through "C0073" increasing sequentially. 
new_col_names <- paste0("C", sprintf("%03d", seq_along(peak_area_cols)))
colnames(OCG_GC_2012)[peak_area_cols] <- new_col_names 

#subset to have just the columns that contain "RT" 1:74. if going to try graphing   
# RT_cols <- grep("RT", colnames(OCG_GC_2012))
# 
# new_col_names_RT <- paste0("RT", sprintf("%03d", seq_along(RT_cols)))
# colnames(OCG_GC_2012)[RT_cols] <- new_col_names_RT

names(OCG_GC_2012)[names(OCG_GC_2012) == 'Plant_ID'] <- 'Garden Plant ID' 

OCG_GC_2012 <- merge(md.OCG.2012, OCG_GC_2012, by="Garden Plant ID") #156 obs of 90

OCG_GC_2012 <- OCG_GC_2012[,-c(1:15)] #removing everything except area under the curve 156 obs of variables

# OCG_GC_2012_w_RT <- OCG_GC_2012[,-c(1:15)] #removing everything except area under the curve and RT. 156 obs of 149 variables

##Save csv with RT
# write.csv(OCG_GC_2012_w_RT, file = "data_csv/OCG_GC_2012_cleaned_w_RT.csv",row.names = FALSE)

##Save csv
write.csv(OCG_GC_2012, file = "data_csv/OCG_GC_2012_cleaned.csv",row.names = FALSE)

## 2021 GC RAW DATA READ IN
OCG_GC_2021 <- read.csv("data_csv/OCG_2021_GCData.csv", head=T, skip = 1) #100 of 225 variables

## 2021 GC CLEANING
sum(duplicated(OCG_GC_2021$Plant_ID)) #3 duplicates (100 and then cocktails and blanks)

#sample 100 was ionized during the first run but then we ran out of gas and the second run is the 100 we actually want to keep so we can remove the first row
OCG_GC_2021 <- OCG_GC_2021[-1, ] #99 of 225 variables

# Remove rows where "Blank" or "Cocktail" or 429 or 333 appear in Plant ID column 
OCG_GC_2021 <- OCG_GC_2021[!(apply(OCG_GC_2021, 1, function(row) any(grepl("^Blank_|^Cocktail_|^333|^429|^BLANK|Cocktail ", row)))), ] #87 obs of 225 variables

#removes "RT". comment out if needed
OCG_GC_2021 <- OCG_GC_2021 [, grepl("Peak.Area|Plant_ID", colnames(OCG_GC_2021))] #87 obs of 149 variables

#subset to keep Peak Area.
OCG_GC_2021 <- OCG_GC_2021 [, !grepl("Peak.Area.Percent", colnames(OCG_GC_2021))] #87 obs of 75 variables

#remove column with compound 1 to match to the 2012 GC data that starts at compound C002. fixed
#OCG_GC_2021 <- OCG_GC_2021[,-2] #87 of 74 variables

#subset to "Peak Area" 1:74
peak_area_cols <- grep("Peak.Area", colnames(OCG_GC_2021)) 

#the new column "C001" through "C0074" increasing sequentially. 
new_col_names <- paste0("C", sprintf("%03d", seq_along(peak_area_cols)))

#Rename
colnames(OCG_GC_2021)[peak_area_cols] <- new_col_names #87 obs and 74 variables

# #subset to "RT" 1:74 if keeping in RT
# RT_cols <- grep("RT", colnames(OCG_GC_2021))
# 
# #the new column "C001" through "C0074" increasing sequentially.
# new_col_names_RT <- paste0("RT", sprintf("%03d", seq_along(RT_cols)))
# 
# #Rename
# colnames(OCG_GC_2021)[RT_cols] <- new_col_names_RT #87 obs and 151 variables

names(OCG_GC_2021)[names(OCG_GC_2021) == 'Plant_ID'] <- 'Garden Plant ID' 

md.OCG.2021$`Garden Plant ID` <- as.character(md.OCG.2021$`Garden Plant ID`)

OCG_GC_2021 <- merge(md.OCG.2021, OCG_GC_2021, by="Garden Plant ID") #70 obs of 89 variables

#Going to remove everything except AUC and plant ID description and RT
#OCG_GC_2021_w_RT <- OCG_GC_2021[,-c(1:15,17:18)] #70 obs of 149 variables

#Going to remove everything except AUC and plant ID description
OCG_GC_2021 <- OCG_GC_2021[,-c(1:15)] #70 obs of 75 variables

# #Save csv with RT
# write.csv(OCG_GC_2021_w_RT, file = "data_csv/OCG_GC_2021_cleaned_w_RT.csv",row.names = FALSE)

#Save csv
write.csv(OCG_GC_2021, file = "data_csv/OCG_GC_2021.csv",row.names = FALSE)

#FULL GC COMBINED
OCG_GC <- rbind(OCG_GC_2012, OCG_GC_2021) #226 obs of 75 variables

#FULL GC COMBINED W RT
#OCG_GC_w_RT <- rbind(OCG_GC_2012_w_RT, OCG_GC_2021_w_RT) #226 obs of 149 variables

#restructure the RT GC full data. 
#remove plant ID?
#make the row names the plant id
rownames(OCG_GC_w_RT) <- OCG_GC_w_RT[,1]
#pivot longer
df_long <- OCG_GC_w_RT %>%
  pivot_longer(cols = everything(), 
               names_to = c(".value", "Compound"), 
               names_pattern = "(C|RT)0*(\\d+)")

# Rename columns for clarity
names(df_long) <- c("Compound", "RT", "PeakArea")

# ##Save csv
# write.csv(OCG_GC_w_RT, file = "data_csv/OCG_GC_w_RT_full_clean.csv")

OCG_GC_w_RT_2012 <- read.csv("data_csv/OCG_GC_2012_cleaned_w_RT.csv", head=T, check.names = F,stringsAsFactors = T, row.names = 1)
OCG_GC_w_RT_2021 <- read.csv("data_csv/OCG_GC_2021_cleaned_w_RT.csv", head=T, check.names = F,stringsAsFactors = T, row.names = 1)
OCG_GC_w_RT <- rbind(OCG_GC_w_RT_2012, OCG_GC_w_RT_2021) #226 obs of 149 variables
OCG_GC_w_RT <- read.csv("data_csv/OCG_GC_w_RT_full_clean.csv", head=T, check.names = F,stringsAsFactors = T)

##Save csv
write.csv(OCG_GC, file = "data_csv/OCG_GC_full_clean.csv",row.names = FALSE)

## 3UL LCMS RAW DATA READ IN
OCG_LCMS_3uL <- read.csv("data_csv/3uL_injection_results_LCMS.csv", head=T, check.names = F,stringsAsFactors = T, skip = 1) #120 of 929 variables

sum(duplicated(OCG_LCMS_3uL$Plant_ID)) #0

## LCMS 3UL CLEANING
OCG_LCMS_3uL <- OCG_LCMS_3uL [, grepl("Peak.Area|Plant_ID", colnames(OCG_LCMS_3uL))] #120 obs of 617 variables

OCG_LCMS_3uL <- OCG_LCMS_3uL [, !grepl("Peak.Area.Percent", colnames(OCG_LCMS_3uL))] #120 obs of 309 variables

#now i'm going to have to rename all of the peak areas by their compound name? unless there is a better way to retain those when reading in the data. The compounds are in order 10, 100-109, 11, 110-119, 12, 120-129, 13.... etc. but not all the way throughout so I think I will just create a key until I have found an alternate way to look at the compound name

#1:308
peak_area_cols <- grep("Peak.Area", colnames(OCG_LCMS_3uL)) 

#"C001" through "C0074"
new_col_names <- paste0("C", sprintf("%03d", seq_along(OCG_LCMS_3uL)))

colnames(OCG_LCMS_3uL)[peak_area_cols] <- new_col_names #120 obs of 309 variables 

#column names are now C001-C309

# Splitting the Plant_ID column and creating new columns
OCG_LCMS_3uL <- OCG_LCMS_3uL %>%
  separate(Plant_ID, into = c("Plant_ID_Number", "Year"), sep = "_")

# Convert 'Year' column to integer
OCG_LCMS_3uL$Year <- as.integer(OCG_LCMS_3uL$Year)

#dataframe now includes the plant ID number and year as their own columns to match with the metadata. 
names(OCG_LCMS_3uL)[names(OCG_LCMS_3uL) == 'Plant_ID_Number'] <- 'Garden Plant ID' 

#2012 LCMS
OCG_LCMS_3uL_2012 <- subset(OCG_LCMS_3uL, OCG_LCMS_3uL$Year=="2012") #46 observations and 310 variables

OCG_LCMS_3uL_2012 <- merge(md.OCG.2012, OCG_LCMS_3uL_2012, by="Garden Plant ID") #41 obs of 332 variables

OCG_LCMS_3uL_2012 <- OCG_LCMS_3uL_2012[,-c(1:15,17)] #41 obs of 309 variables

#2021 LCMS
OCG_LCMS_3uL_2021 <- subset(OCG_LCMS_3uL, OCG_LCMS_3uL$Year=="2021") #73 observations and 310 variables

OCG_LCMS_3uL_2021 <- merge(md.OCG.2021, OCG_LCMS_3uL_2021, by="Garden Plant ID") #71 obs of 325

OCG_LCMS_3uL_2021 <- OCG_LCMS_3uL_2021[,-c(1:15,17)] #71 obs of 309 variables

#FULL LCMS
OCG_LCMS_3uL <- data.frame(rbind(OCG_LCMS_3uL_2012,OCG_LCMS_3uL_2021)) #112 of 309 var

##Save csv
write.csv(OCG_LCMS_3uL, file = "data_csv/OCG_LCMS_3uL_cleaned.csv",row.names = FALSE)

#Clear Global Environment
rm(list = ls())

# Cleaned data read in: START HERE ####
#METADATA
md <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #505 obs of 16 variables.
md <- md[order(row.names(md)),]

### Subsetting to just the plant in the common garden (OCG)
md.OCG <- subset(md, md$Project=="OCG") #246 of 16 variables

### Remove duplicates, negatve controls, and MTW.3.7.R_2012
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1', 'ORT.1.5_2012', 'NEG_8-28-21', 'NEG_10-2-20', 'MTW.3.7.R_2012')
md.OCG <- md.OCG[!rownames(md.OCG) %in% rows_to_remove, ] #235 of 16 var

#make variables factor to plot and droplevels
md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Location","Description")] <- lapply(md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Location","Description")], as.factor)
md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Location","Description")] <- lapply(md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Location","Description")], droplevels)
str(md.OCG)

#subset the md to only have observations from 2012 to avoid duplicates
md.OCG.2012 <- subset(md.OCG, md.OCG$Year=="2012") #159
str(md.OCG.2012)

#subset the md to only have observations from 2021 to avoid duplicates
md.OCG.2021 <- subset(md.OCG, md.OCG$Year=="2021") #76
str(md.OCG.2021)

##FULL CLEAN GC
OCG_GC <- read.csv("data_csv/OCG_GC_full_clean.csv", row.names = 1)#226 obs of 74 variables
OCG_GC <- OCG_GC[order(row.names(OCG_GC)),]
OCG_GC <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(md.OCG)) #217 of 74 variables
md.OCG.GC <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC)) #217 of 16 variables
OCG_GC[is.na(OCG_GC)] <- 0


## 2012 CLEAN GC
OCG_GC_2012 <- read.csv("data_csv/OCG_GC_2012_cleaned.csv", row.names = 1) #157 obs of 74 variables
OCG_GC_2012 <- OCG_GC_2012[order(row.names(OCG_GC_2012)),]
OCG_GC_2012 <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md.OCG)) #147 of 74 variables
md.OCG.GC.2012 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC_2012)) #147 of 16 variables
OCG_GC_2012[is.na(OCG_GC_2012)] <- 0

## 2021 CLEAN GC
OCG_GC_2021 <- read.csv("data_csv/OCG_GC_2021.csv", row.names = 1)#70 obs of 74 variables
OCG_GC_2021 <- OCG_GC_2021[order(row.names(OCG_GC_2021)),] 
OCG_GC_2021 <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md.OCG)) #70 of 74 variables
md.OCG.GC.2021 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC_2021)) #70 of 16 variables
OCG_GC_2021[is.na(OCG_GC_2021)] <- 0

## CLEAN LCMS 
OCG_LCMS_3uL <- read.csv("data_csv/OCG_LCMS_3uL_cleaned.csv", row.names = 1) #112 obs of 308 var
OCG_LCMS_3uL <- OCG_LCMS_3uL[order(row.names(OCG_LCMS_3uL)),]
OCG_LCMS_3uL <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.OCG)) #111 of 308 variables
md.OCG.LCMS.3 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_LCMS_3uL)) #111 of 16 variables
OCG_LCMS_3uL[is.na(OCG_LCMS_3uL)] <- 0

# Alpha diversity ####
## GC alpha diversity ####
OCG.GC.shannon <- diversity(OCG_GC)
OCG.GC.ef <- exp(OCG.GC.shannon)
OCG.GC.ef.r <- round(OCG.GC.ef)

md.OCG.GC <- cbind(md.OCG.GC, effective_species = OCG.GC.ef.r)

glm.OCG.GC <- glm(effective_species ~ Year + Subspecies + Ploidy, family = poisson, data = md.OCG.GC)
summary(glm.OCG.GC) #year and subspecies sig

glm.OCG.GC.gamma <- glm(effective_species ~ Subspecies + Year + Ploidy, family = Gamma, data=md.OCG.GC)
summary(glm.OCG.GC.gamma) #year and subspecies is significant

plot(allEffects(glm.OCG.GC))

plot(md.OCG.GC$Year,OCG.GC.ef)
plot(md.OCG.GC$Subspecies,OCG.GC.ef)

ggplot(md.OCG.GC, aes(Year, effective_species))+
  geom_boxplot(aes(group = Year, fill = Year))+
  theme_classic()

ggplot(data = md.OCG.GC, mapping = aes(x = Subspecies, y = effective_species, fill = Subspecies)) +
  geom_boxplot() +
  theme_classic()

## LCMS alpha diversity ####
OCG.LCMS.shannon <- diversity(OCG_LCMS_3uL)
OCG.LCMS.ef <- exp(OCG.LCMS.shannon)
OCG.LCMS.ef.r <- round(OCG.LCMS.ef)

md.OCG.LCMS.3 <- cbind(md.OCG.LCMS.3, effective_species = OCG.LCMS.ef.r)
glm.OCG.LCMS <- glm(effective_species ~ Subspecies + Year, family = poisson, data = md.OCG.LCMS.3)
summary(glm.OCG.LCMS)

glm.OCG.LCMS.gamma <- glm(effective_species ~ Subspecies + Year + Ploidy + Location, family = Gamma, data=md.OCG.LCMS.3)
summary(glm.OCG.LCMS.gamma) #location, subspecies and ploidy is significant

plot(allEffects(glm.OCG.LCMS))

plot(md.OCG.LCMS.3$Year,OCG.LCMS.ef)
plot(md.OCG.LCMS.3$Subspecies,OCG.LCMS.ef)

ggplot(md.OCG.LCMS.3, aes(Year, effective_species))+
  geom_boxplot(aes(group = Year, fill = Year))+
  theme_classic()

ggplot(data = md.OCG.LCMS.3, mapping = aes(x = Subspecies, y = effective_species, fill = Subspecies)) +
  geom_boxplot() +
  theme_classic()


# Cleaning for PCAs####
## 2012 GC threshold defined ####
OCG_GC_2012[OCG_GC_2012 == 0] <- NA
colSums(is.na(OCG_GC_2012)) #checking for null values since they are used for filtering

#So now the data needs to be cleaned to only contain compounds that occur in more than 20% of the samples (plants). AKA columns need to be removed that don't have at least 20% of the rows containing a number=NA

#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_GC_2012))
print(na_proportion)

#Now define the threshold of 10% - there are 73 compounds that remain after this
threshold <- 0.90

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2012_subset <- OCG_GC_2012[, columns_to_keep] #147 obs of 74 variables

# Replace NA values with zeroes
OCG_GC_2012_subset[is.na(OCG_GC_2012_subset)] <- 0
colSums(is.na(OCG_GC_2012_subset))

## 2021 GC threshold defined ####
OCG_GC_2021[OCG_GC_2021 == 0] <- NA
colSums(is.na(OCG_GC_2021))

na_proportion <- colMeans(is.na(OCG_GC_2021))
print(na_proportion)

threshold <- 0.90

columns_to_keep <- na_proportion <= threshold

OCG_GC_2021_subset <- OCG_GC_2021[, columns_to_keep] #70 obs of 43 var variables

OCG_GC_2021_subset[is.na(OCG_GC_2021_subset)] <- 0
colSums(is.na(OCG_GC_2021_subset)) 

## Full GC threshold defined #### 
OCG_GC[OCG_GC == 0] <- NA
colSums(is.na(OCG_GC)) 

na_proportion <- colMeans(is.na(OCG_GC))
print(na_proportion)

threshold <- 0.90

columns_to_keep <- na_proportion <= threshold

OCG_GC_subset <- OCG_GC[, columns_to_keep] #217 obs of 54 variables

OCG_GC_subset[is.na(OCG_GC_subset)] <- 0

colSums(is.na(OCG_GC_subset)) 

## LCMS threshold defined ####
OCG_LCMS_3uL[OCG_LCMS_3uL == 0] <- NA
colSums(is.na(OCG_LCMS_3uL))

na_proportion <- colMeans(is.na(OCG_LCMS_3uL))
print(na_proportion)

threshold <- 0.90

columns_to_keep <- na_proportion <= threshold

OCG_LCMS_3uL_subset <- OCG_LCMS_3uL[, columns_to_keep] #111 of 302 variables

OCG_LCMS_3uL_subset[is.na(OCG_LCMS_3uL_subset)] <- 0

colSums(is.na(OCG_LCMS_3uL_subset)) 

# Scaling####
## 2012 GC scaling ####
#Compound
data_normalized_2012 <- scale(OCG_GC_2012_subset) #this will center the data 1:147, 1:47

#2012 plant ID
OCG_GC_2012_subset.t <- t(OCG_GC_2012_subset) #transpose the data first to put plant ID as columns and compound as row
colSums(is.na(OCG_GC_2012_subset.t)) #checking for null values
data_normalized_2012.ID <- scale(OCG_GC_2012_subset.t) #1:47, 1:147

## 2021 GC scaling ####
#Compound 
data_normalized_2021 <- scale(OCG_GC_2021_subset) #1:70, 1:42

#2021 plant ID
OCG_GC_2021_subset.t <- t(OCG_GC_2021_subset) 
colSums(is.na(OCG_GC_2021_subset.t))
data_normalized_2021.ID <- scale(OCG_GC_2021_subset.t) #1:42, 1:70

## Full GC scaling ####
#Compound
data_normalized_GC <- scale(OCG_GC_subset) #1:217, 1:54

#Plant ID
OCG_GC_subset.t <- t(OCG_GC_subset)
colSums(is.na(OCG_GC_subset.t))
data_normalized_GC.ID <- scale(OCG_GC_subset.t) #1:54, 1:217

## GC scaling for both years ####
# OCG_GC_yr12<- subset(OCG_GC_2012_subset, row.names(OCG_GC_2012_subset) %in% row.names(md.OCG.2012)) #147 of 47
# 
# OCG_GC_yr21 <- subset(OCG_GC_2021_subset, row.names(OCG_GC_2021_subset) %in% row.names(md.OCG.2021)) #70 of 42 variables

# OCG_GC_yr12 <- merge(OCG_GC_yr12, md.OCG.2012, by = "row.names", all.x = TRUE) #147 of 64
# OCG_GC_yr21 <- merge(OCG_GC_yr21, md.OCG.2021, by = "row.names", all.x = TRUE) #147 of 59

OCG_GC_btyr <- subset(OCG_GC_subset, row.names(OCG_GC_subset) %in% row.names(md.OCG.GC)) #217 of 54. 
OCG_GC_btyr <- merge(OCG_GC_btyr, md.OCG.GC, by = "row.names", all.x = TRUE) #217 of 71

# Identify plant IDs with duplicates in both years
duplicated_plant_ids <- OCG_GC_btyr[duplicated(OCG_GC_btyr$Plant),] #65 of 71 var

#subset to include only 'TRUE' for paired
OCG_GC_btyr <- OCG_GC_btyr[OCG_GC_btyr$Paired == TRUE, ] #139 of 71 var

# Plant names to remove since they didnt meet the threshold parameters
plants_to_remove <- c("WAT.2.8", "WAT.2.4", "UTWV.2.10", "UTW.1.10", "UTV.3.5", "UTT.1.1", "MTT.1.6", "IDW.1.6", "CAT.2.9")

# Filter to exclude specified plant names
OCG_GC_btyr <- OCG_GC_btyr %>%
  filter(!Plant %in% plants_to_remove) #130 of 71 variables

#cleaning to remove everything except plant ID and compounds
#5671
OCG_GC_btyr <- OCG_GC_btyr[, -c(56:71)] #130 of 55

#make the row names the plant id
rownames(OCG_GC_btyr) <- OCG_GC_btyr[,1]

# Remove the column with the rownames from the dataframe
OCG_GC_btyr <- OCG_GC_btyr[,-1 ] #130 obs of 54 var

#Compound
data_normalized_GC_btyr <- scale(OCG_GC_btyr) #1:130, 1:54
OCG_GC_btyr.t <- t(OCG_GC_btyr)
data_normalized_GC_btyr.ID <- scale(OCG_GC_btyr.t) #1:54, 1:130

## LCMS scaling ####
#Compound
data_LCMS3_normalized <- scale(OCG_LCMS_3uL_subset) #1:111, 1:302

#Plant ID
OCG_LCMS_3uL.t <- t(OCG_LCMS_3uL_subset)
colSums(is.na(OCG_LCMS_3uL.t)) 
data_LCMS3_normalized.ID <- scale(OCG_LCMS_3uL.t) #1:302, 1:111 

## Correlation ####
### 2012 GC CORRELATION
#Compound
corr_matrix_2012 <- cor(data_normalized_2012) #47 compounds
#ggcorrplot(corr_matrix_2012)
#Plant ID
corr_matrix_2012.ID <- cor(data_normalized_2012.ID) #147 plants

### 2021 GC CORRELATION 
#Compound
corr_matrix_2021 <- cor(data_normalized_2021) #42 compounds
#ggcorrplot(corr_matrix_2021)
#Plant ID
corr_matrix_2021.ID <- cor(data_normalized_2021.ID) # 70 plants

### GC CORRELATION
#Compound
corr_matrix_GC <- cor(data_normalized_GC) #54 compounds
#ggcorrplot(corr_matrix_GC)
#Plant ID
corr_matrix_GC.ID <- cor(data_normalized_GC.ID) # 217 plants

###GC CORRELATION BOTH YEAR
corr_matrix_GC_btyr <- cor(data_normalized_GC_btyr) #55 compounds
#ggcorrplot(corr_matrix_GC)
#Plant ID
corr_matrix_GC_btyr.ID <- cor(data_normalized_GC_btyr.ID) # 217 plants

### LCMS CORRELATION 
#Compound
corr_matrix_LCMS3 <- cor(data_LCMS3_normalized)#large matrix 91204 elements
#ggcorrplot(corr_matrix_LCMS3)
#Plant ID
corr_matrix_LCMS3.ID <- cor(data_LCMS3_normalized.ID) #111 plants

# PCA ####
## 2012 GC PCA ####
data.pca_2012 <- princomp(corr_matrix_2012)
summary(data.pca_2012)
data.pca_2012$loadings[, 1:2]

fviz_eig(data.pca_2012, addlabels = TRUE) #scree plot is used to visualize the importance of each principal component and can be used to determine the number of principal components to retain. (52%, 17%)

#With the biplot, it is possible to visualize the similarities and dissimilarities between the samples, and further shows the impact of each attribute on each of the principal components.

# Graph of the compounds
fviz_pca_var(data.pca_2012, col.var = "black")

#Contribution of each compound
fviz_cos2(data.pca_2012, choice = "var", axes = 1:2, xtickslab.rt = 90, top = 20)

#Biplot combined with cos2
fviz_pca_var(data.pca_2012, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

data.pca.GC_2012 <- princomp(data_normalized_2012)
biplot(data.pca.GC_2012, cex = 0.4, pc.biplot = TRUE)

## 2012 GC PCA BY PLANT ID
#prcomp() has improved numerical accuracy, so is preferable to use this function.
data.pca_2012_ID <- prcomp(data_normalized_2012)
summary(data.pca_2012_ID)
fviz_eig(data.pca_2012_ID, addlabels = TRUE) #scree plot. 19.8%, 12.5% 
fviz_cos2(data.pca_2012_ID, choice = "ind", axes = 1:2) #Contribution of each plant
autoplot(data.pca_2012_ID)
autoplot(data.pca_2012_ID, label = TRUE)

rownames(data_normalized_2012) == rownames(md.OCG.GC.2012)

#BY PLOIDY
plot(data.pca_2012_ID$x[, 1], data.pca_2012_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC comp of 2012 plant by ploidy", 
     col= c("red","blue")[md.OCG.GC.2012$Ploidy],
     pch=c(19),
     xlim = range(data.pca_2012_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2012_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_2012_ID$x[, 1], data.pca_2012_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC of 2012 plant by subspecies", 
     col= c("pink","brown",'darkgreen')[md.OCG.GC.2012$Subspecies],
     pch=c(19),
     xlim = range(data.pca_2012_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2012_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("pink","brown","darkgreen"),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_2012_ID$x[, 1], data.pca_2012_ID$x[, 2],
     xlab="PC 1 (19.81%)", ylab="PC 2 (12.46%)", 
     main="PCA of 2012 GC data by subspecies and ploidy", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC.2012$Subsp_ploidy],
     pch=c(19),
     xlim = range(data.pca_2012_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2012_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

### ANOVA for 2012 GC subspecies ploidy ####
pca_scores_GC12 <- data.pca_2012_ID$x
GC12_aov_df <- as.data.frame(pca_scores_GC12)
pca_model_GC_12_subsppl <- aov(cbind(PC1, PC2) ~ md.OCG.GC.2012$Subsp_ploidy, data = GC12_aov_df)
summary(pca_model_GC_12_subsppl) #PC 1 1.067e-15, PC2 0.0006886 of 2 dof. F(4) = 32.893, p < 0.001 for PC1. F(4) = 5.5143, p < 0.001

#### k- means clustering on 2012 GC####
set.seed(92)
pca_scores_GC12 <- data.pca_2012_ID$x 

#ELBOW METHOD
# plot the within cluster sum of squares against the number of clusters
wcss <- numeric(10)
for (i in 1:10) {
  kmeans_model <- kmeans(pca_scores_GC12, centers = i)
  wcss[i] <- kmeans_model$tot.withinss
}
plot(1:10, wcss, type = "b", xlab = "Number of Clusters (k)", ylab = "WCSS")

#Look for the "elbow" point where the rate of decrease in WCSS slows down significantly.
#he elbow point is a good estimate for the optimal k.

k <- 5  # Adjust this number based on your analysis

# Perform k-means clustering on the PCA scores
kmeans_result <- kmeans(pca_scores_GC12, centers = k)

# Get cluster assignments for each sample
cluster_assignments <- kmeans_result$cluster

# Visualize the clusters (optional)
plot(pca_scores_GC12, col = cluster_assignments, main = "PCA of 2012 GC data with k-means Clustering")

## Adding k means clusters to md 
md.OCG.GC.2012$cluster_assignments <- cluster_assignments

plot(data.pca_2012_ID$x[, 1], data.pca_2012_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="PCA of 2012 GC data with k-means cluster (5 clusters)", 
     pch = 19,
     col= c("lightcoral","rosybrown",'darkseagreen','peachpuff',"darkturquoise")[md.OCG.GC.2012$cluster_assignments],
     xlim = range(data.pca_2012_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2012_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("1","2","3","4","5"),
       col= c("lightcoral","rosybrown",'darkseagreen','peachpuff',"darkturquoise"),
       pch=19,
       cex=0.8,
       bty = "n")

#SPECTRAL CLUSTERING 2012 GC 
set.seed(75)
scgc12 <- specc(pca_scores_GC12, centers = 4)
scgc12
centers(scgc12)
size(scgc12)
withinss(scgc12)
plot(pca_scores_GC12, col = scgc12, main = "PCA of 2012 GC data with spectral clustering")

## 2021 GC PCA ####
data.pca_2021 <- princomp(data_normalized_2021)
summary(data.pca_2021)
data.pca_2021$loadings[, 1:2]

fviz_eig(data.pca_2021, addlabels = TRUE) #(49%, 16%)

fviz_pca_var(data.pca_2021, col.var = "black")

fviz_cos2(data.pca_2021, choice = "var", axes = 1:2, xtickslab.rt = 90, top = 20)

fviz_pca_var(data.pca_2021, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

data.pca.GC_2021 <- princomp(data_normalized_2021)
biplot(data.pca.GC_2021, cex = 0.4, pc.biplot = TRUE)

## 2021 GC PCA BY PLANT ID
data.pca_2021_ID <- prcomp(data_normalized_2021)
summary(data.pca_2021_ID)
fviz_eig(data.pca_2021_ID, addlabels = TRUE) #20%, 13%
fviz_cos2(data.pca_2021_ID, choice = "ind", axes = 1:2)
autoplot(data.pca_2021_ID)
autoplot(data.pca_2021_ID, label = TRUE)

rownames(data_normalized_2021) == rownames(md.OCG.GC.2021)

#BY PLOIDY
plot(data.pca_2021_ID$x[, 1], data.pca_2021_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC of 2021 plant by ploidy", 
     col= c("red","blue")[md.OCG.GC.2021$Ploidy],
     pch=c(19),
     xlim = range(data.pca_2021_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2021_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES
plot(data.pca_2021_ID$x[, 1], data.pca_2021_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC of 2021 plant by subspecies", 
     col= c("pink","brown",'darkgreen')[md.OCG.GC.2021$Subspecies],
     pch=c(19),
     xlim = range(data.pca_2021_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2021_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("pink","brown","darkgreen"),
       pch=19,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES PLOIDY
plot(data.pca_2021_ID$x[, 1], data.pca_2021_ID$x[, 2],
     xlab="PC 1 (20.54%)", ylab="PC 2 (13.33%)", 
     main="PCA of 2021 GC data by subspecies and ploidy", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC.2021$Subsp_ploidy],
     pch=c(19),
     xlim = range(data.pca_2021_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2021_ID$x[, 2], na.rm = TRUE))
legend("bottomright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

### ANOVA for 2021 GC subspecies ploidy ####
pca_scores_GC21 <- data.pca_2021_ID$x
GC21_aov_df <- as.data.frame(pca_scores_GC21)
pca_model_GC_21_subsppl <- aov(cbind(PC1, PC2) ~ md.OCG.GC.2021$Subsp_ploidy, data = GC21_aov_df)
summary(pca_model_GC_21_subsppl) #F(4) = 4.5663 p < 0.005 (0.002598) for PC1. F(4) = 5.3852, p < 0.001 

#### k- means clustering attempt on 2021 GC####
pca_scores_GC21 <- data.pca_2021_ID$x

#ELBOW METHOD
# plot the within cluster sum of squares against the number of clusters
wcss <- numeric(10)
for (i in 1:10) {
  kmeans_model <- kmeans(pca_scores_GC21, centers = i)
  wcss[i] <- kmeans_model$tot.withinss
}
plot(1:10, wcss, type = "b", xlab = "Number of Clusters (k)", ylab = "WCSS")

#Look for the "elbow" point where the rate of decrease in WCSS slows down significantly.
#he elbow point is a good estimate for the optimal k.

k <- 5  # Adjust this number based on your analysis

# Perform k-means clustering on the PCA scores
kmeans_result <- kmeans(pca_scores_GC21, centers = k)

# Get cluster assignments for each sample
cluster_assignments <- kmeans_result$cluster

# Visualize the clusters (optional)
plot(pca_scores_GC21, col = cluster_assignments, main = "PCA of 2021 GC data with k-means Clustering")

## Adding k means clusters to md 
md.OCG.GC.2021$cluster_assignments <- cluster_assignments

plot(data.pca_2021_ID$x[, 1], data.pca_2021_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="PCA of 2021 GC data with k-means cluster (5 clusters)", 
     pch = 19,
     col= c("lightcoral","rosybrown",'darkseagreen','peachpuff',"darkturquoise")[md.OCG.GC.2021$cluster_assignments],
     xlim = range(data.pca_2021_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2021_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("1","2","3","4","5"),
       col= c("lightcoral","rosybrown",'darkseagreen','peachpuff',"darkturquoise"),
       pch=19,
       cex=0.8,
       bty = "n")

#SPECTRAL CLUSTERING 2021 GC 
set.seed(24)
scgc21 <- specc(pca_scores_GC21, centers = 4)
scgc21
centers(scgc21)
size(scgc21)
withinss(scgc21)
plot(pca_scores_GC21, col = scgc21, main = "PCA of 2021 GC data with spectral clustering")

## Full GC PCA ####
data.pca_GC_ID <- prcomp(data_normalized_GC)
summary(data.pca_GC_ID)
fviz_eig(data.pca_GC_ID, addlabels = TRUE) #17.6, 10.7%
fviz_cos2(data.pca_GC_ID, choice = "var", axes = 1:2, xtickslab.rt = 90, top = 20) 
autoplot(data.pca_GC_ID)
autoplot(data.pca_GC_ID, label = TRUE)

data.pca_GC <- princomp(data_normalized_GC)
biplot(data.pca_GC, cex = 0.4, pc.biplot = TRUE)

rownames(data_normalized_GC) == rownames(md.OCG.GC)

#BY PLOIDY
plot(data.pca_GC_ID$x[, 1], data.pca_GC_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="Full GC comp plant by ploidy", 
     col= c("red","blue")[md.OCG.GC$Ploidy],
     pch=c(19),
     xlim = range(data.pca_GC_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES
plot(data.pca_GC_ID$x[, 1], data.pca_GC_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="Full GC comp plant by subspecies", 
     col= c("pink","brown",'darkgreen')[md.OCG.GC$Subspecies],
     pch=c(19),
     xlim = range(data.pca_GC_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("pink","brown","darkgreen"),
       pch=19,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES PLOIDY
#### This is a nice plot to share!
plot(data.pca_GC_ID$x[, 1], data.pca_GC_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="PCA of GC data by subspecies, ploidy, and year", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC$Subsp_ploidy],
     pch=c(17,19)[md.OCG.GC$Year],
     xlim = range(data.pca_GC_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.6,
       bty = "n")
legend("topleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.6,
       bty = "n")

OCG_GC_subset[is.na(OCG_GC_subset)] <- 0
summary(rowSums(OCG_GC_subset)) #179 seqs in smallest sample
summary(colSums(OCG_GC_subset)) #2500
OCG_GC_subset <- OCG_GC_subset[,colSums(OCG_GC_subset) > 0]
summary(colSums(OCG_GC_subset)) #2500

OCG_GC_subset.r <- rrarefy(round(OCG_GC_subset),sample = 179)

#permanova by subspecies ploidy
PCA_gc_subsploi <- adonis2(OCG_GC_subset.r ~ md.OCG.GC$Subsp_ploidy, by = "margin")
PCA_gc_subsploi #subspecies ploidy is significant 0.001

#BY YEAR
plot(data.pca_GC_ID$x[, 1], data.pca_GC_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="PCA of GC data by year", 
     col= c("maroon","cyan")[md.OCG.GC$Year],
     pch=19)
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

PCA_gc_yr <- adonis2(OCG_GC_subset.r ~ md.OCG.GC$Year, by = "margin")
PCA_gc_yr #year is significant 0.001

## GC PCA for existing in both years ####
data.pca_GC_btyr_ID <- prcomp(data_normalized_GC_btyr)
summary(data.pca_GC_btyr_ID)
fviz_eig(data.pca_GC_btyr_ID, addlabels = TRUE) #22% and 9.7%
fviz_cos2(data.pca_GC_btyr_ID, choice = "var", axes = 1:2, xtickslab.rt = 90, top = 20) #Contribution of each compound
fviz_cos2(data.pca_GC_btyr_ID, choice = "ind", axes = 1:2, xtickslab.rt = 90, top = 20) #Contribution of each plant
autoplot(data.pca_GC_btyr_ID)
autoplot(data.pca_GC_btyr_ID, label = TRUE)

md.OCG.GC.btyr <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC_btyr)) #130

rownames(OCG_GC_btyr) == rownames(md.OCG.GC.btyr) 

#BY PLOIDY
plot(data.pca_GC_btyr_ID$x[, 1], data.pca_GC_btyr_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC plants in both years by ploidy", 
     col= c("red","blue")[md.OCG.GC.btyr$Ploidy],
     pch=c(19),
     xlim = range(data.pca_GC_btyr_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_btyr_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES
plot(data.pca_GC_btyr_ID$x[, 1], data.pca_GC_btyr_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC of plants from both years by subspecies", 
     col= c("pink","brown",'darkgreen')[md.OCG.GC.btyr$Subspecies],
     pch=c(19),
     xlim = range(data.pca_GC_btyr_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_btyr_ID$x[, 2], na.rm = TRUE))
legend("topleft",
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("pink","brown","darkgreen"),
       pch=19,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES PLOIDY
plot(data.pca_GC_btyr_ID$x[, 1], data.pca_GC_btyr_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="PCA of GC from plants in both years by subspecies, ploidy, and year", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC.btyr$Subsp_ploidy],
     pch=c(17,19)[md.OCG.GC.btyr$Year],
     xlim = range(data.pca_GC_btyr_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_btyr_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.6,
       bty = "n")
legend("topleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.6,
       bty = "n")

#permanova for subspecies ploidy
PCA_GCbtyr_subsploi <- adonis2(OCG_GC_btyr ~ md.OCG.GC.btyr$Subsp_ploidy, by = "margin")
PCA_GCbtyr_subsploi #subspecies ploidy is significant 0.001

#BY YEAR
plot(data.pca_GC_btyr_ID$x[, 1], data.pca_GC_btyr_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC of plants in both years by year", 
     col= c("maroon","cyan")[md.OCG.GC.btyr$Year],
     pch=c(19),
     xlim = range(data.pca_GC_btyr_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_btyr_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

## Full LCMS PCA ####
data_LCMS3_normalized.t <- t(data_LCMS3_normalized)
data.pca_LCMS3 <- princomp(data_LCMS3_normalized.t)
summary(data.pca_LCMS3)

data.pca_LCMS3 <- prcomp(data_LCMS3_normalized.ID)
summary(data.pca_LCMS3)
fviz_cos2(data.pca_LCMS3, choice = "var", axes = 1:2, xtickslab.rt = 90, top = 20) #Contribution of each plant
fviz_cos2(data.pca_LCMS3, choice = "ind", axes = 1:2, xtickslab.rt = 90, top = 20) #Contribution of each compound

data.pca_LCMS3_ID <- prcomp(data_LCMS3_normalized)
summary(data.pca_LCMS3_ID)
fviz_eig(data.pca_LCMS3_ID, addlabels = TRUE) #14.1% and 8.9%
fviz_cos2(data.pca_LCMS3_ID, choice = "var", axes = 1:2, xtickslab.rt = 90, top = 20) #Contribution of each compound
fviz_cos2(data.pca_LCMS3_ID, choice = "ind", axes = 1:2, xtickslab.rt = 90, top = 20) #Contribution of each plant
autoplot(data.pca_LCMS3_ID)
autoplot(data.pca_LCMS3_ID, label = TRUE)

biplot(data.pca_LCMS3_ID, cex = 0.3, pc.biplot = TRUE)

# data.pca_LCMS3 <- princomp(data_LCMS3_normalized)
# autoplot(data.pca_LCMS3, label = TRUE)

# merged_data <- merge(data.frame(Row.names = rownames(data.pca_LCMS3_ID$x), data.pca_LCMS3_ID$x), OCG_LCMS_3uL_subset)
# 
# ggplot(merged_data)+
#   geom_point(aes(PC1, PC2, color = "Row.names"))+
#   geom_point(aes(C001,C002))

rownames(data_LCMS3_normalized) == rownames(md.OCG.LCMS.3)

#BY PLOIDY
plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS 3uL plant by ploidy", 
     col= c("red","blue")[md.OCG.LCMS.3$Ploidy],
     pch=c(19),
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES
plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS 3uL plant by subspecies", 
     col= c("pink","brown",'darkgreen')[md.OCG.LCMS.3$Subspecies],
     pch=c(19),
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topleft",
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("pink","brown","darkgreen"),
       pch=19,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES PLOIDY
## This is a good plot to share!
plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1 (13.51%)", ylab="PC 2 (9.51%)", 
     main="PCA of LCMS data by subspecies, ploidy, and year", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.LCMS.3$Subsp_ploidy],
     pch=c(17,19)[md.OCG.LCMS.3$Year],
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
legend("topleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

summary(rowSums(OCG_LCMS_3uL_subset)) #31016795 seqs in smallest sample
summary(colSums(OCG_LCMS_3uL_subset)) #384380
OCG_LCMS_3uL_subset <- OCG_LCMS_3uL_subset[,colSums(OCG_LCMS_3uL_subset) > 0]
summary(colSums(OCG_LCMS_3uL_subset)) #384380

#OCG_LCMS_3uL_subset.r <- rrarefy(round(OCG_LCMS_3uL_subset),sample = 31016795)

#permanova for subspecies ploidy
PCA_lcms_subsploi.r <- adonis2(OCG_LCMS_3uL_subset.r ~ md.OCG.LCMS.3$Subsp_ploidy, by = "margin")
PCA_lcms_subsploi.r #subspecies ploidy is significant 0.001

### ANOVA for LCMS subspecies ploidy ####
pca_scores_LCMS <- data.pca_LCMS3_ID$x
LCMS_subspl_df <- as.data.frame(pca_scores_LCMS)
pca_model_LCMS_subsppl <- aov(cbind(PC1, PC2) ~ md.OCG.LCMS.3$Subsp_ploidy, data = LCMS_subspl_df)
summary(pca_model_LCMS_subsppl) #F(4) = 34.767 p < 0.001 for PC1. F(4) = 9.9945, p < 0.001 

#BY YEAR
plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS data by year", 
     col= c("maroon","cyan")[md.OCG.LCMS.3$Year],
     pch=c(19),
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

#BY LOCATION
plotcolor <- c("firebrick","cadetblue","sienna","skyblue","rosybrown","tomato","olivedrab","turquoise","burlywood","mediumaquamarine","darkseagreen")

plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1 (13.51%)", ylab="PC 2 (9.51%)", 
     main="PCA of LCMS data by location and subspecies", 
     col= plotcolor[md.OCG.LCMS.3$Location],
     pch=c(17,19,15)[md.OCG.LCMS.3$Subspecies],
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.8,
       bty = "n")
legend("topleft", 
       legend=c("T","V","W"),
       col= "black",
       pch=c(17,19,15),
       cex=0.8,
       bty = "n")

PCA_lcms_loc.r <- adonis2(OCG_LCMS_3uL_subset.r ~ md.OCG.LCMS.3$Location, by = "margin")
PCA_lcms_loc.r #subspecies ploidy is significant 0.001

# k- means clustering for LCMS####
set.seed(65)
pca_scores_LCMS <- data.pca_LCMS3_ID$x

#ELBOW METHOD
# plot the within cluster sum of squares against the number of clusters
wcss <- numeric(10)
for (i in 1:10) {
  kmeans_model <- kmeans(pca_scores_LCMS, centers = i)
  wcss[i] <- kmeans_model$tot.withinss
}
plot(1:10, wcss, type = "b", xlab = "Number of Clusters (k)", ylab = "WCSS")

#Look for the "elbow" point where the rate of decrease in WCSS slows down significantly.
#he elbow point is a good estimate for the optimal k.
k <- 4

# Perform k-means clustering on the PCA scores
kmeans_result <- kmeans(pca_scores_LCMS, centers = k)

# Get cluster assignments for each sample
cluster_assignments <- kmeans_result$cluster

# Visualize the clusters (optional)
plot(pca_scores_LCMS, col = cluster_assignments, main = "PCA of LCMS data with k-means Clustering")

## Adding k means clusters to md 
md.OCG.LCMS.3$cluster_assignments <- cluster_assignments

plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="PCA of LCMS data with k-means cluster (4 clusters)", 
     pch = 19,
     col= c("lightcoral","rosybrown",'darkseagreen','peachpuff')[md.OCG.LCMS.3$cluster_assignments],
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("1","2","3","4"),
       col= c("lightcoral","rosybrown",'darkseagreen','peachpuff'),
       pch=19,
       cex=0.8,
       bty = "n")

##DBSCAN ####
# Perform DBSCAN clustering
dbscan::kNNdistplot(OCG_LCMS_3uL_subset, k =  5)
abline(h = 6.0e+06, lty = 2)
set.seed(123)
# fpc package
res.fpc <- fpc::dbscan(OCG_LCMS_3uL_subset, eps = 6.0e+06 , MinPts = 2)
plot(res.fpc, pca_scores_LCMS)
# dbscan package
res.db <- dbscan::dbscan(pca_scores_LCMS, 14.5, 2)
all(res.fpc$cluster == res.db$cluster) #TRUE

fviz_cluster(res.fpc, pca_scores_LCMS, geom = "point")

set.seed(1)
dbscan_result <- dbscan(pca_scores_LCMS, eps = 14.5, 2)
dbscan_result$cluster
fviz_cluster(dbscan_result, data = pca_scores_LCMS, geom = "point", 
             outlier.pointsize = 1, main = "DBSCAN Clustering on LCMS PCA-transformed Data", palette = "RdYlGn", shape = 19)+ theme_classic() 

#HDBSCAN
hdbscan_result <- hdbscan(pca_scores_LCMS, minPts = 5)
plot(hdbscan_result, col = hdbscan_result$cluster+1, pch = 20)

#SPECTRAL CLUSTERING
set.seed(43)
sc <- specc(pca_scores_LCMS, centers = 5)
sc
centers(sc)
size(sc)
withinss(sc)
plot(pca_scores_LCMS, col = sc, main = "PCA of LCMS data with spectral clustering")

##LCMS PCA of just tridentata ####
#subset LCMS to just tridentata
OCG_LCMS_tri <- subset(OCG_LCMS_3uL, md.OCG.LCMS.3$Subspecies=="T") #80 of 308 variables
md.OCG.LCMS.tri <- subset(md.OCG.LCMS.3, row.names(md.OCG.LCMS.3) %in% row.names(OCG_LCMS_tri)) #80

#make 0 NA to redefine threshold
OCG_LCMS_tri[OCG_LCMS_tri == 0] <- NA
colSums(is.na(OCG_LCMS_tri))
na_proportion <- colMeans(is.na(OCG_LCMS_tri))
print(na_proportion)
threshold <- 0.90
columns_to_keep <- na_proportion <= threshold
OCG_LCMS_tri_subset <- OCG_LCMS_tri[, columns_to_keep] #80 0f 297 var
OCG_LCMS_tri_subset[is.na(OCG_LCMS_tri_subset)] <- 0
colSums(is.na(OCG_LCMS_tri_subset)) 

#SCALING
#Compound
data_LCMS_tri_normalized <- scale(OCG_LCMS_tri_subset) #1:80, 1:297

#Plant ID
OCG_LCMS_tri.t <- t(OCG_LCMS_tri_subset)
colSums(is.na(OCG_LCMS_tri.t)) 
data_LCMS_tri_normalized.ID <- scale(OCG_LCMS_tri.t) #1:302, 1:111

#PCA
data.pca_LCMS_tri_ID <- prcomp(data_LCMS_tri_normalized)
summary(data.pca_LCMS_tri_ID)
fviz_eig(data.pca_LCMS_tri_ID, addlabels = TRUE) #17.1% & 9.3%%
fviz_cos2(data.pca_LCMS_tri_ID, choice = "var", axes = 1:2, xtickslab.rt = 90, top = 20) #Contribution of each compound
fviz_cos2(data.pca_LCMS_tri_ID, choice = "ind", axes = 1:2, xtickslab.rt = 90, top = 20) #Contribution of each plant
autoplot(data.pca_LCMS_tri_ID)
autoplot(data.pca_LCMS_tri_ID, label = TRUE)

biplot(data.pca_LCMS_tri_ID, cex = 0.3, pc.biplot = TRUE)

rownames(data_LCMS_tri_normalized) == rownames(md.OCG.LCMS.tri)

#BY PLOIDY
plot(data.pca_LCMS_tri_ID$x[, 1], data.pca_LCMS_tri_ID$x[, 2],
     xlab="PC 1 (17.1%)", ylab="PC 2", 
     main="LCMS data of tridentata plants by ploidy", 
     col= c("red","blue")[md.OCG.LCMS.tri$Ploidy],
     pch=c(19),
     xlim = range(data.pca_LCMS_tri_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS_tri_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(data.pca_LCMS_tri_ID,groups = md.OCG.LCMS.tri$Ploidy, show.groups = "2n", col = "red")
ordispider(data.pca_LCMS_tri_ID,groups = md.OCG.LCMS.tri$Ploidy, show.groups = "4n", col = "blue")

#BY SUBSPECIES PLOIDY
## This is a good plot to share!
plot(data.pca_LCMS_tri_ID$x[, 1], data.pca_LCMS_tri_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="PCA of LCMS data of tridentata by ploidy, and year", 
     col= c("pink","lightblue")[md.OCG.LCMS.tri$Subsp_ploidy],
     pch=c(17,19)[md.OCG.LCMS.3$Year],
     xlim = range(data.pca_LCMS_tri_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS_tri_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("T_2n","T_4n"),
       col= c("pink",'lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
legend("topleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")
ordispider(data.pca_LCMS_tri_ID,groups = md.OCG.LCMS.tri$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(data.pca_LCMS_tri_ID,groups = md.OCG.LCMS.tri$Subsp_ploidy, show.groups = "T_4n", col = "lightblue")

summary(rowSums(OCG_LCMS_tri_subset)) #49580982 seqs in smallest sample
summary(colSums(OCG_LCMS_tri_subset)) #110226
OCG_LCMS_tri_subset <- OCG_LCMS_tri_subset[,colSums(OCG_LCMS_tri_subset) > 0]
summary(colSums(OCG_LCMS_tri_subset)) #110226

#OCG_LCMS_tri_subset.r <- rrarefy(round(OCG_LCMS_tri_subset),sample = 49580982) #takes forever to run

#permanova for subspecies ploidy
#PCA_lcms_tri_subsploi.r <- adonis2(OCG_LCMS_tri_subset.r ~ md.OCG.LCMS.tri$Subsp_ploidy, by = "margin")
#PCA_lcms_tri_subsploi.r #subspecies ploidy is significant 0.001

#BY YEAR
plot(data.pca_LCMS_tri_ID$x[, 1], data.pca_LCMS_tri_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS data of tridentata by year", 
     col= c("maroon","cyan")[md.OCG.LCMS.tri$Year],
     pch=c(19),
     xlim = range(data.pca_LCMS_tri_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS_tri_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

PCA_lcms_tri_yr.r <- adonis2(OCG_LCMS_tri_subset.r ~ md.OCG.LCMS.tri$Year, by = "margin")
PCA_lcms_tri_yr.r #year is significant 0.006

#BY LOCATION
md.OCG.LCMS.tri$Location <- factor(md.OCG.LCMS.tri$Location)
levels(md.OCG.LCMS.tri$Location)
md.OCG.LCMS.tri$Location <- droplevels(md.OCG.LCMS.tri$Location)
levels(md.OCG.LCMS.tri$Location)

plotcolor.tri <- c("firebrick","cadetblue","skyblue","rosybrown","tomato","olivedrab","turquoise","burlywood","mediumaquamarine")

plot(data.pca_LCMS_tri_ID$x[, 1], data.pca_LCMS_tri_ID$x[, 2],
     xlab="PC 1 (17.08%)", ylab="PC 2 (9.3%)", 
     main="PCA of LCMS data by location for tridentata", 
     col= plotcolor.tri[md.OCG.LCMS.tri$Location],
     pch=19,
     xlim = range(data.pca_LCMS_tri_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS_tri_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("AZ","CA", "ID", "MT", "NM", "NV", "OR", "UT", "WA"),
       col= plotcolor.tri,
       pch=19,
       cex=0.8,
       bty = "n")

#PCA_lcms_tri.loc.r <- adonis2(OCG_LCMS_tri_subset ~ md.OCG.LCMS.tri$Location, by = "margin")
#PCA_lcms_tri.loc.r #location and ploidy are significant 0.001

### ANOVA for LCMS loc ####
pca_scores_LCloc.tri <- data.pca_LCMS_tri_ID$x
LCMS.tri_aov_df <- as.data.frame(pca_scores_LCloc.tri)
pca_model_LCMS.tri <- aov(cbind(PC1, PC2) ~ md.OCG.LCMS.tri$Location, data = LCMS.tri_aov_df)
summary(pca_model_LCMS.tri) #F(8) = 68.684 p < 0.001 for PC1. F(8) = 2.1365, p < 0.005 

##LCMS PCA of just wyomingensis ####
#subset LCMS to just wyomingensis
OCG_LCMS_wy <- subset(OCG_LCMS_3uL, md.OCG.LCMS.3$Subspecies=="W") #80 of 308 variables
md.OCG.LCMS.wy <- subset(md.OCG.LCMS.3, row.names(md.OCG.LCMS.3) %in% row.names(OCG_LCMS_wy)) #80

#make 0 NA to redefine threshold
OCG_LCMS_wy[OCG_LCMS_wy == 0] <- NA
colSums(is.na(OCG_LCMS_wy))
na_proportion <- colMeans(is.na(OCG_LCMS_wy))
print(na_proportion)
threshold <- 0.90
columns_to_keep <- na_proportion <= threshold
OCG_LCMS_wy_subset <- OCG_LCMS_wy[, columns_to_keep] #80 0f 297 var
OCG_LCMS_wy_subset[is.na(OCG_LCMS_wy_subset)] <- 0
colSums(is.na(OCG_LCMS_wy_subset)) 

#SCALING
#Compound
data_LCMS_wy_normalized <- scale(OCG_LCMS_wy_subset) #1:26, 1:299

#Plant ID
OCG_LCMS_wy.t <- t(OCG_LCMS_wy_subset)
colSums(is.na(OCG_LCMS_wy.t)) 
data_LCMS_wy_normalized.ID <- scale(OCG_LCMS_wy.t) #1:299, 1:26

#PCA
data.pca_LCMS_wy_ID <- prcomp(data_LCMS_wy_normalized)
summary(data.pca_LCMS_wy_ID)
fviz_eig(data.pca_LCMS_wy_ID, addlabels = TRUE) #17.1% & 9.3%%
fviz_cos2(data.pca_LCMS_wy_ID, choice = "var", axes = 1:2, xtickslab.rt = 90, top = 20) #Contribution of each compound
fviz_cos2(data.pca_LCMS_wy_ID, choice = "ind", axes = 1:2, xtickslab.rt = 90, top = 20) #Contribution of each plant
autoplot(data.pca_LCMS_wy_ID)
autoplot(data.pca_LCMS_wy_ID, label = TRUE)

biplot(data.pca_LCMS_wy_ID, cex = 0.3, pc.biplot = TRUE)

rownames(data_LCMS_wy_normalized) == rownames(md.OCG.LCMS.wy)

#BY LOCATION
md.OCG.LCMS.wy$Location <- factor(md.OCG.LCMS.wy$Location)
levels(md.OCG.LCMS.wy$Location)
md.OCG.LCMS.tri$Location <- droplevels(md.OCG.LCMS.tri$Location)
levels(md.OCG.LCMS.tri$Location)

plotcolor.wy <- c("sienna","skyblue","rosybrown","turquoise","burlywood","mediumaquamarine")

plot(data.pca_LCMS_wy_ID$x[, 1], data.pca_LCMS_wy_ID$x[, 2],
     xlab="PC 1 (19.14%)", ylab="PC 2 (15.41%)", 
     main="PCA of LCMS data by location for wyomingensis", 
     col= plotcolor.wy[md.OCG.LCMS.wy$Location],
     pch=19,
     xlim = range(data.pca_LCMS_wy_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS_wy_ID$x[, 2], na.rm = TRUE))
legend("bottomright", 
       legend=c("CO", "ID", "MT", "OR", "UT", "WA"),
       col= plotcolor.wy,
       pch=19,
       cex=0.8,
       bty = "n")

PCA_lcms_wy.loc.r <- adonis2(OCG_LCMS_wy_subset ~ md.OCG.LCMS.wy$Location, by = "margin")
PCA_lcms_wy.loc.r #location is significant 0.001

### ANOVA for LCMS loc wyomingensis ####
pca_scores_LCloc.w <- data.pca_LCMS_wy_ID$x
LCMS_wy_aov_df <- as.data.frame(pca_scores_LCloc.w)
pca_model_LCloc.w <- aov(cbind(PC1, PC2) ~ md.OCG.LCMS.wy$Location, data = LCMS_wy_aov_df)
summary(pca_model_LCloc.w) #F(5) = 10.065 p < 0.001 for PC1. F(5) = 0.5863, p = 0.71


# NMDS plots####
## 2012 GC NMDS####
str(md.OCG.2012)
OCG_GC_2012[is.na(OCG_GC_2012)] <- 0
summary(rowSums(OCG_GC_2012)) #7962 seqs in smallest sample
summary(colSums(OCG_GC_2012)) #0
OCG_GC_2012 <- OCG_GC_2012[,colSums(OCG_GC_2012) > 0]
summary(colSums(OCG_GC_2012)) #77.1

OCG_GC_2012.r <- rrarefy(round(OCG_GC_2012),sample = 7962)

set.seed(37)
#OCG_AUC_2012_ID.nmds <- metaMDS(OCG_GC_2012.r, k=3, trymax=1000) #solution reached only with k=3
#save(OCG_AUC_2012_ID.nmds, file = "nmds/OCG_AUC_2012_ID.nmds.rda")
load("nmds/OCG_AUC_2012_ID.nmds.rda")

ordiplot(OCG_AUC_2012_ID.nmds, type = "t",display = "sites",cex = .5)

rownames(md.OCG.GC.2012) == rownames(OCG_AUC_2012_ID.nmds$points)

#PLOIDY
plot(OCG_AUC_2012_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of 2012 plant by ploidy", 
     col= c("red","blue")[md.OCG.GC.2012$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

#### PERMANOVA FOR 2012 PLOIDY ##
OCG_GC_2012_ploidy.r <- adonis2(OCG_GC_2012.r ~ md.OCG.GC.2012$Ploidy,by="margin",na.rm = T)
OCG_GC_2012_ploidy.r #ploidy is significant 0.001

OCG_GC_2012_ploidy_subsp.r <- adonis2(OCG_GC_2012.r ~ md.OCG.GC.2012$Ploidy + md.OCG.GC.2012$Subspecies,by="margin",na.rm = T)
OCG_GC_2012_ploidy.r #ploidy and subspecies is significant 0.001

#SUBSPECIES
plot(OCG_AUC_2012_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="2012 GC chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.GC.2012$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")

#### PERMANOVA FOR SUBSPECIES
OCG_GC_2012_subsp.r <- adonis2(OCG_GC_2012.r ~ md.OCG.GC.2012$Subspecies,by="margin", na.rm=T) 
OCG_GC_2012_subsp.r #subspecies is signficant= 0.001

#PAIRWISE ADONIS
OCG_GC_2012_subsp.pw.r <- pairwise.adonis(OCG_GC_2012.r, md.OCG.GC.2012$Subspecies)
OCG_GC_2012_subsp.pw.r # sig between all subspecies

# SUBSPECIES PLOIDY
plot(OCG_AUC_2012_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="NMDS of 2012 GC data by subspecies and ploidy", 
     col= c("pink","brown","darkgreen",'tan','lightblue')[md.OCG.GC.2012$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_AUC_2012_ID.nmds,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(OCG_AUC_2012_ID.nmds,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(OCG_AUC_2012_ID.nmds,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(OCG_AUC_2012_ID.nmds,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(OCG_AUC_2012_ID.nmds,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

# PERMANOVAS FOR SUBSPECIES PLOIDY
OCG_GC_2012_subspploidy <- adonis2(OCG_GC_2012.r ~ md.OCG.GC.2012$Subsp_ploidy,by="margin") 
OCG_GC_2012_subspploidy #subspecies ploidy is significant

#PAIRWISE ADONIS
OCG_GC_2012_subsp_ploidy.pw.r <- pairwise.adonis(OCG_GC_2012.r, md.OCG.GC.2012$Subsp_ploidy)
OCG_GC_2012_subsp_ploidy.pw.r 

## 2021 GC NMDS #### 
OCG_GC_2021[is.na(OCG_GC_2021)] <- 0
OCG_GC_2021 <- OCG_GC_2021[-which(rownames(OCG_GC_2021) == "CAT.1.1_2021"), ] #remove this outlier
summary(rowSums(OCG_GC_2021)) #2842 seqs in smallest sample
summary(colSums(OCG_GC_2021)) #0
OCG_GC_2021 <- OCG_GC_2021[,colSums(OCG_GC_2021) > 0]
summary(colSums(OCG_GC_2021)) #100.3

md.OCG.GC.2021 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC_2021)) #70 of 16 variables

OCG_GC_2021.r <- rrarefy(round(OCG_GC_2021),sample = 2842)

set.seed(47)
#OCG_GC_2021.nmds <- metaMDS(OCG_GC_2021.r, trymax=500) #solution reached
#save(OCG_GC_2021.nmds, file = "nmds/OCG_GC_2021.nmds.rda")
load("nmds/OCG_GC_2021.nmds.rda")

ordiplot(OCG_GC_2021.nmds, type = "t",display = "sites",cex = .7)

rownames(md.OCG.GC.2021) == rownames(OCG_GC_2021.nmds$points)

#PLOIDY
plot(OCG_GC_2021.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of 2021 plant by ploidy", 
     col= c("red","blue")[md.OCG.GC.2021$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

#### PERMANOVA FOR 2021 PLOIDY##
OCG_GC_2021_ploidy.r <- adonis2(OCG_GC_2021.r ~ md.OCG.GC.2021$Ploidy,by="margin")
OCG_GC_2021_ploidy.r #ploidy is significant 0.002

#SUBSPECIES
plot(OCG_GC_2021.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="2021 GC chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.GC.2021$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")

#### PERMANOVA FOR SUBSPECIES
OCG_GC_2021_subsp.r <- adonis2(OCG_GC_2021.r ~ md.OCG.GC.2021$Subspecies,by="margin") 
OCG_GC_2021_subsp.r #subspecies is signficant= 0.001

#PAIRWISE ADONIS SUBSPECIES
OCG_GC_2021_subsp.pw.r <- pairwise.adonis(OCG_GC_2021.r, md.OCG.GC.2021$Subspecies)
OCG_GC_2021_subsp.pw.r # sig between T vs V, T vs W

# SUBSPECIES PLOIDY NMDS
plot(OCG_GC_2021.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="NMDS of 2021 GC data by subspecies and ploidy", 
     col= c("pink","brown","darkgreen",'tan','lightblue')[md.OCG.GC.2021$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_GC_2021.nmds,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(OCG_GC_2021.nmds,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(OCG_GC_2021.nmds,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(OCG_GC_2021.nmds,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(OCG_GC_2021.nmds,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

# PERMANOVAS FOR SUBSPECIES PLOIDY
OCG_GC_2021_subspploidy.r <- adonis2(OCG_GC_2021.r ~ md.OCG.GC.2021$Subsp_ploidy,by="margin") 
OCG_GC_2021_subspploidy.r #subspecies ploidy is significant

#PAIRWISE ADONIS FOR SUBSPECIES PLOIDY
OCG_GC_2021_subsp_ploidy.pw <- pairwise.adonis(OCG_GC_2021.r, md.OCG.GC.2021$Subsp_ploidy)
OCG_GC_2021_subsp_ploidy.pw 

## Full GC NMDS ####
OCG_GC[is.na(OCG_GC)] <- 0
OCG_GC <- OCG_GC[-which(rownames(OCG_GC) == "CAT.1.1_2021"), ] #remove this outlier
summary(rowSums(OCG_GC)) #2842 seqs in smallest sample
summary(colSums(OCG_GC)) #0
OCG_GC <- OCG_GC[,colSums(OCG_GC) > 0]
summary(colSums(OCG_GC)) #77.1

OCG_GC.r <- rrarefy(round(OCG_GC),sample = 2842)

md.OCG.GC <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC)) #216

rownames(md.OCG.GC) == rownames(OCG_GC_ID.nmds$points)

set.seed(64)
#OCG_GC_ID.nmds <- metaMDS(OCG_GC.r, k=3, trymax=500) #solution not reached
#save(OCG_GC_ID.nmds, file = "nmds/OCG_GC_ID.nmds.rda")
load("nmds/OCG_GC_ID.nmds.rda")

ordiplot(OCG_GC_ID.nmds, type = "t", display = "sites", cex = .4)

#PLOIDY
plot(OCG_GC_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of plant by ploidy", 
     col= c("red","blue")[md.OCG.GC$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

#### PERMANOVA FOR PLOIDY ##
OCG_GC_ploidy.r <- adonis2(OCG_GC.r ~ md.OCG.GC$Ploidy,by="margin")
OCG_GC_ploidy.r #ploidy is significant 0.003

#SUBSPECIES
plot(OCG_GC_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="GC chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.GC$Subspecies],
     pch=c(19))
legend("topright", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.6,
       bty = "n")

#### PERMANOVA FOR SUBSPECIES###
OCG_GC_subsp.r <- adonis2(OCG_GC.r ~ md.OCG.GC$Subspecies,by="margin") 
OCG_GC_subsp.r #subspecies is signficant= 0.001

OCG_GC_subsp_yr.r <- adonis2(OCG_GC.r ~ md.OCG.GC$Subspecies*md.OCG.GC$Year,by="margin") 
OCG_GC_subsp_yr.r #sig 

#PAIRWISE ADONIS FOR SUBSPECIES
OCG_GC_subsp.pw.r <- pairwise.adonis(OCG_GC.r, md.OCG.GC$Subspecies)
OCG_GC_subsp.pw.r # sig between T vs V, V vs W, & T VS W

rownames(OCG_GC_ID.nmds$points) == rownames(md.OCG.GC)

# SUBSPECIES PLOIDY YEAR NMDS 
# This is a nice plot to share
plot(OCG_GC_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="NMDS of GC data by subspecies, ploidy, and year", 
     col= c("pink","brown","darkgreen",'tan','lightblue')[md.OCG.GC$Subsp_ploidy],
     pch=c(17,19)[md.OCG.GC$Year])
legend("bottomright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.6,
       bty = "n")
legend("topright", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.6,
       bty = "n")

# PERMANOVAS FOR SUBSPECIES PLOIDY
OCG_GC_subspploidy <- adonis2(OCG_GC.r ~ md.OCG.GC$Subsp_ploidy,by="margin") 
OCG_GC_subspploidy #subspecies ploidy is significant

OCG_GC_subspploidy_yr <- adonis2(OCG_GC.r ~ md.OCG.GC$Subsp_ploidy+md.OCG.GC$Year,by="margin")
OCG_GC_subspploidy_yr #sig

#PAIRWISEADONIS FOR SUBSPECIES PLOIDY
OCG_GC_subsp_ploidy.pw <- pairwise.adonis(OCG_GC.r, md.OCG.GC$Subsp_ploidy)
OCG_GC_subsp_ploidy.pw 

## Trying to find name mismatch
rownames(OCG_GC_ID.nmds$points) == rownames(md.OCG.GC)

#YEAR
plot(OCG_GC_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="NMDS of GC data by year", 
     col= c("maroon","cyan")[md.OCG.GC$Year],
     pch=c(19))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")
text(OCG_GC_ID.nmds$points[,1:2],
     labels=md.OCG.GC$Description,
     pos=1,
     cex=0.4)

OCG_GC_yr <- adonis2(OCG_GC.r ~ md.OCG.GC$Year,by="margin") 
OCG_GC_yr #yearis signficant= 0.001

## Full LCMS NMDS ####
OCG_LCMS_3uL[is.na(OCG_LCMS_3uL)] <- 0

summary(rowSums(OCG_LCMS_3uL)) #31016795 LCMS area in smallest sample
summary(colSums(OCG_LCMS_3uL)) #265605

OCG_LCMS_3uL.r <- rrarefy(round(OCG_LCMS_3uL),sample = 31016795)

set.seed(38)
#OCG_LCMS_3uL_ID.nmds <- metaMDS(OCG_LCMS_3uL.r, trymax=500) #solution reached
#save(OCG_LCMS_3uL_ID.nmds, file = "nmds/OCG_LCMS_3uL_ID.nmds.rda")
load("nmds/OCG_LCMS_3uL_ID.nmds.rda")

ordiplot(OCG_LCMS_3uL_ID.nmds, type = "t",display = "sites",cex = .5)

md.OCG.LCMS.3 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_LCMS_3uL)) #subset md to match AUC samples #111

row.names(md.OCG.LCMS.3) == row.names(OCG_LCMS_3uL_ID.nmds$points)

#PLOIDY
plot(OCG_LCMS_3uL_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 3uL plant by ploidy", 
     col= c("red","blue")[md.OCG.LCMS.3$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

#### PERMANOVA PLOIDY
OCG_LCMS3_ploidy.r <- adonis2(OCG_LCMS_3uL.r ~ md.OCG.LCMS.3$Ploidy ,by="margin") 
OCG_LCMS3_ploidy.r #ploidy is significant 0.001

#SUBSPECIES
plot(OCG_LCMS_3uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 3uL chemistry by subspecies and ploidy",
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.LCMS.3$Subspecies],
     pch=c(17,19)[md.OCG.LCMS.3$Ploidy])
legend("topright", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomleft", 
       legend=c("2n","4n"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

#### PERMANOVA FOR SUBSPECIES
OCG_LCMS3_subsp <- adonis2(OCG_LCMS_3uL.r ~ md.OCG.LCMS.3$Subspecies,by="margin") 
OCG_LCMS3_subsp #subspecies is signficant= 0.001

#PAIRWISEADONIS
OCG_LCMS3_subsp.pw <- pairwise.adonis(OCG_LCMS_3uL.r, md.OCG.LCMS.3$Subspecies)
OCG_LCMS3_subsp.pw # sig between all subspecies

# SUBSPECIES PLOIDY 
plot(OCG_LCMS_3uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="NMDS of LCMS data by subspecies, ploidy, and year", 
     col= c("pink","brown","darkgreen",'tan','lightblue')[md.OCG.LCMS.3$Subsp_ploidy],
     pch=c(17,19)[md.OCG.LCMS.3$Year])
legend("bottomleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
legend("topleft", 
       legend=c("2012","2021"),
       col="black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

# PERMANOVAS FOR SUBSPECIES PLOIDY
OCG_LCMS3_subspploidy <- adonis2(OCG_LCMS_3uL.r ~ md.OCG.LCMS.3$Subsp_ploidy, by="margin") 
OCG_LCMS3_subspploidy #subspecies ploidy is significant

#PAIRWISE ADONIS
OCG_LCMS3_subsp_ploidy.pw <- pairwise.adonis(OCG_LCMS_3uL.r, md.OCG.LCMS.3$Subsp_ploidy)
OCG_LCMS3_subsp_ploidy.pw 

#YEAR
plot(OCG_LCMS_3uL_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 3uL plant by year", 
     col= c("maroon","cyan")[md.OCG.LCMS.3$Year],
     pch=c(19))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

# PERMANOVAS
OCG_LCMS3_year <- adonis2(OCG_LCMS_3uL.r ~ md.OCG.LCMS.3$Year ,by="margin") 
OCG_LCMS3_year #year is significant

# LOCATION 
plot(OCG_LCMS_3uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="NMDS of LCMS data by location", 
     col= plotcolor [md.OCG.LCMS.3$Location],
     pch=c(17,16)[md.OCG.LCMS.3$Year])
legend("bottomleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.8,
       bty = "n")
legend("topleft", 
       legend=c("2012","2021"),
       col="black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

## Subset LCMS to tridentata only ####
md.tridentata <- md.OCG.LCMS.3[md.OCG.LCMS.3$Subspecies == "T", ] #80 plants
LCMS.tridentata <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.tridentata)) #80
md.tridentata$Location <- factor(md.tridentata$Location)
levels(md.tridentata$Location)
md.tridentata$Location <- droplevels(md.tridentata$Location)
levels(md.tridentata$Location)

summary(rowSums(LCMS.tridentata)) #49580982 
summary(colSums(LCMS.tridentata)) #10697

LCMS.tridentata.r <- rrarefy(round(LCMS.tridentata),sample = 49580982)

set.seed(8)
#LCMS.tridentata.nmds <- metaMDS(LCMS.tridentata.r, trymax=500) #solution reached
#save(LCMS.tridentata.nmds, file = "nmds/LCMS.tridentata.nmds.rda")
load("nmds/LCMS.tridentata.nmds.rda")

rownames(LCMS.tridentata.nmds$points) == rownames(md.tridentata)

ordiplot(LCMS.tridentata.nmds, type = "t",display = "sites",cex = .5)

# LOCATION 
plot(LCMS.tridentata.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="NMDS of LCMS data by location and year for tridentata", 
     col= plotcolor.tri [md.tridentata$Location],
     pch=c(17,16)[md.tridentata$Year])
legend("bottomleft", 
       legend=c("AZ","CA", "ID", "MT", "NM", "NV", "OR", "UT", "WA"),
       col= plotcolor.tri,
       pch=19,
       cex=0.8,
       bty = "n")
legend("topleft", 
       legend=c("2012","2021"),
       col="black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

tri_LCMS3_loc <- adonis2(LCMS.tridentata.r ~ md.tridentata$Location ,by="margin") 
tri_LCMS3_loc #location is significant

lm_lcms_t <- lm(LCMS.tridentata.r ~ md.tridentata$Location)
summary(lm_lcms_t)

#YEAR
plot(LCMS.tridentata.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="NMDS of LCMS data by ploidy for tridentata", 
     col= c("pink","lightblue")[md.tridentata$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n"),
       col= c("pink","lightblue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(LCMS.tridentata.nmds,groups = md.tridentata$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(LCMS.tridentata.nmds,groups = md.tridentata$Subsp_ploidy, show.groups = "T_4n", col = "lightblue")

##Subset LCMS to wyomingensis only ####
md.wyomingensis <- md.OCG.LCMS.3[md.OCG.LCMS.3$Subspecies == "W", ] #26 plants

LCMS.wyomingensis <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.wyomingensis)) #26

summary(rowSums(LCMS.wyomingensis)) #31016795 
summary(colSums(LCMS.wyomingensis)) #0

LCMS.wyomingensis.r <- rrarefy(round(LCMS.wyomingensis),sample = 31016795)

set.seed(98)
LCMS.wyomingensis.nmds <- metaMDS(LCMS.wyomingensis.r, trymax=500) #solution reached
save(LCMS.wyomingensis.nmds, file = "nmds/LCMS.wyomingensis.nmds.rda")
load("nmds/LCMS.wyomingensis.nmds.rda")

ordiplot(LCMS.wyomingensis.nmds, type = "t",display = "sites",cex = .5)

# LOCATION 
plot(LCMS.wyomingensis.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="NMDS of LCMS data by location for wyomingensis", 
     col= plotcolor [md.wyomingensis$Location],
     pch=c(17,19)[md.wyomingensis$Year])
legend("topright", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.8,
       bty = "n")
# legend("topleft",
#        legend=c("2012","2021"),
#        col="black",
#        pch=c(17,19),
#        cex=0.6,
#        bty = "n")

wy_LCMS3_loc <- adonis2(LCMS.wyomingensis.r ~ md.wyomingensis$Location ,by="margin") 
wy_LCMS3_loc #year is significant

lm_lcms_w <- lm(LCMS.wyomingensis.r ~ md.wyomingensis$Location)
summary(lm_lcms_w)

## 2012 LCMS NMDS ####
str(md.OCG.2012)
#subset to just 2012 data
OCG_LCMS_2012 <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.OCG.2012)) #40
md.OCG.LCMS.2012 <- subset(md.OCG.2012, row.names(md.OCG.2012) %in% row.names(OCG_LCMS_2012)) #40

summary(rowSums(OCG_LCMS_2012)) #50540924 area in smallest sample
summary(colSums(OCG_LCMS_2012)) #0

OCG_LCMS_2012.r <- rrarefy(round(OCG_LCMS_2012),sample = 50540924)
OCG_LCMS_2012[is.na(OCG_LCMS_2012)] <- 0

set.seed(450)
#OCG_LCMS_2012_ID.nmds <- metaMDS(t(m_OCG_LCMS_2012.t), trymax=500) #solution reached
#save(OCG_LCMS_2012_ID.nmds, file = "nmds/OCG_LCMS_2012_ID.nmds.rda")
load("nmds/OCG_LCMS_2012_ID.nmds.rda")

ordiplot(OCG_LCMS_2012_ID.nmds, type = "t",display = "sites",cex = .7)

#PLOIDY
plot(OCG_LCMS_2012_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 2012 3uL plant by ploidy", 
     col= c("red","blue")[md.OCG.LCMS.2012$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

# PERMANOVAS
OCG_LCMS_2012_ploidy <- adonis2(OCG_LCMS_2012 ~ md.OCG.LCMS.2012$Ploidy, by="margin") 
OCG_LCMS_2012_ploidy #ploidy is significant

#by subspecies- only tridentata samples in the 2012 LCMS 3uL data
plot(OCG_LCMS_2012_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 2012 3uL plant by subspecies", 
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.LCMS.2012$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("T","V","W"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")

# Subspecies ploidy NMDS there is T_2n and T_4n 
plot(OCG_LCMS_2012_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 2012 3uL chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[md.OCG.LCMS.2012$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")

# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_LCMS_2012_subspploidy <- adonis2(OCG_LCMS_2012 ~ md.OCG.LCMS.2012$Subsp_ploidy, by="margin") 
OCG_LCMS_2012_subspploidy #subspecies ploidy is significant

#BY LOCATION 
plot(OCG_LCMS_2012_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 2012 3uL chemistry by location",
     col = plotcolor [md.OCG.LCMS.2012$Location],
     pch=19)
legend("bottomleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.6,
       bty = "n")

# PERMANOVAS
OCG_LCMS_2012_loc <- adonis2(OCG_LCMS_2012 ~ md.OCG.LCMS.2012$Location, by="margin") 
OCG_LCMS_2012_loc #location is significant

## 2021 LCMS NMDS ####
str(md.OCG.2021) 
OCG_LCMS_2021 <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.OCG.2021)) #71
md.OCG.LCMS.2021 <- subset(md.OCG.2021, row.names(md.OCG.2021) %in% row.names(OCG_LCMS_2021)) #71

OCG_LCMS_2021[is.na(OCG_LCMS_2021)] <- 0
#OCG_LCMS_2021.t <- t(OCG_LCMS_2021) 
#m_OCG_LCMS_2021.t = as.matrix(OCG_LCMS_2021.t)

set.seed(21)
#OCG_LCMS_2021_ID.nmds <- metaMDS(t(m_OCG_LCMS_2021.t), trymax=500) #solution reached
#save(OCG_LCMS_2021_ID.nmds, file = "nmds/OCG_LCMS_2021_ID.nmds.rda")
load("nmds/OCG_LCMS_2021_ID.nmds.rda")

ordiplot(OCG_LCMS_2021_ID.nmds, type = "t",display = "sites",cex = .7)

#PLOIDY
plot(OCG_LCMS_2021_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 2021 3uL plant by ploidy", 
     col= c("red","blue")[md.OCG.LCMS.2021$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

# PERMANOVAS
OCG_LCMS_2021_ploidy <- adonis2(OCG_LCMS_2021 ~ md.OCG.LCMS.2021$Ploidy, by="margin") 
OCG_LCMS_2021_ploidy #ploidy is significant

#SUBSPECIES
plot(OCG_LCMS_2021_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 2021 3uL plant by subspecies", 
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.LCMS.2021$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("T","V","W"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
# PERMANOVAS
OCG_LCMS_2021_subsp <- adonis2(OCG_LCMS_2021 ~ md.OCG.LCMS.2021$Subspecies, by="margin") 
OCG_LCMS_2021_subsp #subspecies is significant

# SUBSPECIES PLOIDY
plot(OCG_LCMS_2021_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 2021 3uL chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[md.OCG.LCMS.2021$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
# PERMANOVAS
OCG_LCMS_2021_subspploidy <- adonis2(OCG_LCMS_2021 ~ md.OCG.LCMS.2021$Subsp_ploidy, by="margin") 
OCG_LCMS_2021_subspploidy #subspecies ploidy is significant

#LOCATION
plot(OCG_LCMS_2021_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 2021 3uL chemistry by location",
     col = plotcolor [md.OCG.LCMS.2021$Location],
     pch=19)
legend("bottomleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=16,
       cex=0.8,
       bty = "n")
# PERMANOVAS
OCG_LCMS_2021_loc <- adonis2(OCG_LCMS_2021 ~ md.OCG.LCMS.2021$Location, by="margin") 
OCG_LCMS_2021_loc #location is significant

# Binary jaccard plots ####
## LCMS binary jaccard ####
set.seed(8)
#OCG_LCMS_3uL_ID.jdis <- vegdist(OCG_LCMS_3uL.r, method = "jaccard", binary = TRUE)
#OCG_LCMS_3uL_ID.jnmds <- metaMDS(OCG_LCMS_3uL_ID.jdis, trymax=500) #solution reached
#save(OCG_LCMS_3uL_ID.jnmds, file = "jnmds/OCG_LCMS_3uL_ID.jnmds.rda")
load("jnmds/OCG_LCMS_3uL_ID.jnmds.rda")

ordiplot(OCG_LCMS_3uL_ID.jnmds, type = "t",display = "sites",cex = .7)

plotcolor <- c("olivedrab","cadetblue","magenta","blue","orange","green","darkgreen","firebrick","lightgoldenrod","mediumaquamarine","cornflowerblue")

#BY SUBSPECIES AND LOCATION
plot(OCG_LCMS_3uL_ID.jnmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="Binary jaccard of LCMS by subspecies and location",
     col = plotcolor [md.OCG.LCMS.3$Location],
     pch=c(17,19,15)[md.OCG.LCMS.3$Subspecies])
legend("topleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomright", 
       legend=c("T","V","W"),
       col= "black",
       pch=c(17,19,15),
       cex=0.8,
       bty = "n")

#BY PLOIDY AND LOCATION
plot(OCG_LCMS_3uL_ID.jnmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 3uL chemistry by ploidy and location",
     col = plotcolor [md.OCG.LCMS.3$Location],
     pch=c(17,19)[md.OCG.LCMS.3$Ploidy])
legend("topleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomright", 
       legend=c("2n","4n"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

#BY SUBSPECIES PLOIDY AND LOCATION
plot(OCG_LCMS_3uL_ID.jnmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 3uL chemistry by subsp_ploidy and location",
     col = plotcolor [md.OCG.LCMS.3$Location],
     pch=c(17,19,15,18,8)[md.OCG.LCMS.3$Subsp_ploidy])
legend("topleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= "black",
       pch=c(17,19,15,18,8),
       cex=0.8,
       bty = "n")

#BY YEAR AND LOCATION
plot(OCG_LCMS_3uL_ID.jnmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 3uL chemistry by year and location",
     col = plotcolor [md.OCG.LCMS.3$Location],
     pch=c(17,19)[md.OCG.LCMS.3$Year])
legend("topleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomright", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

## GC binary jaccard ####
set.seed(6)
OCG_GC <- OCG_GC[-which(colnames(OCG_GC) == "CAT.1.1_2021"), ] #remove this outlier
OCG_GC_ID.jdis <- vegdist(OCG_GC.r, method = "jaccard", binary = TRUE)
OCG_GC_ID.jnmds <- metaMDS(OCG_GC_ID.jdis, trymax=500) #solution reached
#save(OCG_GC_ID.jnmds, file = "jnmds/OCG_GC_ID.jnmds.rda")
load("jnmds/OCG_GC_ID.jnmds.rda")

ordiplot(OCG_GC_ID.jnmds, type = "t",display = "sites",cex = .7)

### Testing match for row names
rownames(OCG_GC_ID.jnmds$points) == rownames(md.OCG.GC)

#BY PLOIDY
plot(OCG_GC_ID.jnmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of plant by ploidy and year", 
     col= c("red","blue")[md.OCG.GC$Ploidy],
     pch=c(17,19)[md.OCG.GC$Year])
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomright", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n") 


#BY SUBSPECIES
plot(OCG_GC_ID.jnmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="GC chemistry by subspecies & year",
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.GC$Subspecies],
     pch=c(17,19)[md.OCG.GC$Year])
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomright", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n") 

#BY SUBSPECIES PLOIDY
plot(OCG_GC_ID.jnmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Binary jaccard of GC data by subspecies, ploidy, & year", 
     col= c("pink","brown","darkgreen",'tan','lightblue')[md.OCG.GC$Subsp_ploidy],
     pch=c(17,19)[md.OCG.GC$Year])
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomright", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

#BY YEAR
plot(OCG_GC_ID.jnmds$point [,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Binary jaccard of GC data by year", 
     col= c("maroon","cyan")[md.OCG.GC$Year],
     pch=c(19))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

# PCoA ####
## GC PCoA ####
GC_dist_matrix <- vegdist(OCG_GC.r, method = "bray")
GC_pcoa <- cmdscale(GC_dist_matrix) #classic multidimensional scaling (cmdscale)

#BY PLOIDY
plot(GC_pcoa[,1], GC_pcoa[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "GC by ploidy PCoA", 
     col= c("red","blue")[md.OCG.GC$Ploidy],
     pch = 16)
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=16,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES
plot(GC_pcoa[,1], GC_pcoa[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "GC by subspecies PCoA", 
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.GC$Subspecies],
     pch = 16)
legend("topleft", 
       legend=c("T","V","W"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=16,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES PLOIDY
plot(GC_pcoa[,1], GC_pcoa[,2], 
     main="PCoA of GC data by subspecies, ploidy, and year", 
     col= c("pink","brown","darkgreen",'tan','lightblue')[md.OCG.GC$Subsp_ploidy],
     pch=c(17,19)[md.OCG.GC$Year])
legend("bottomright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.6,
       bty = "n")
legend("topright", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.6,
       bty = "n")


#BY YEAR
plot(GC_pcoa[,1], GC_pcoa[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "PCoA of GC by year", 
     col= c("maroon","cyan")[md.OCG.GC$Year],
     pch = 16)
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=16,
       cex=0.8,
       bty = "n")

## LCMS PCoA ####
LCMS_dist_matrix <- vegdist(OCG_LCMS_3uL.r, method = "bray")
LCMS_pcoa <- cmdscale(LCMS_dist_matrix)

plot(LCMS_pcoa)
text(LCMS_pcoa, row.names(LCMS_pcoa), cex = .5)

#BY PLOIDY
plot(LCMS_pcoa[,1], LCMS_pcoa[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "LCMS PCoA by ploidy", 
     col= c("red","blue")[md.OCG.LCMS.3$Ploidy],
     pch = 16)
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=16,
       cex=0.8,
       bty = "n")

#BY YEAR
plot(LCMS_pcoa[,1], LCMS_pcoa[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "LCMS PCoA Visualization by year", 
     col= c("red","cyan")[md.OCG.LCMS.3$Year],
     pch = 16)
legend("topleft", 
       legend=c("2012","2021"),
       col= c("red","cyan"),
       pch=16,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES
plot(LCMS_pcoa[,1], LCMS_pcoa[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "LCMS PCoA Visualization by subspecies", 
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.LCMS.3$Subspecies],
     pch = 16)
legend("topleft", 
       legend=c("T","V","W"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=16,
       cex=0.8,
       bty = "n")

#BY SUBSPECIES PLOIDY
plot(LCMS_pcoa[,1], LCMS_pcoa[,2],
     xlab="PC 1", ylab="PC 2", 
     main="PCoA of LCMS data by subspecies, ploidy, & year", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.LCMS.3$Subsp_ploidy],
     pch=c(17,19)[md.OCG.LCMS.3$Year])
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=16,
       cex=0.8,
       bty = "n")
legend("bottomleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

#BY LOCATION
plot(LCMS_pcoa[,1], LCMS_pcoa[,2],
     xlab="PC 1", ylab="PC 2", 
     main="PCoA of LCMS data by location", 
     col= plotcolor[md.OCG.LCMS.3$Location],
     pch=c(16))
legend("topleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=16,
       cex=0.8,
       bty = "n")

# Clear Global Environment ####
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

# Clear Global Environment ####
rm(list = ls())
# Stable isotope analysis#### 
#2012 
stable_iso_data_2012<-read.csv("data_csv/Stable_isotope_2012.csv")
str(stable_iso_data_2012) #96
# Extracting numbers from Plant_ID column
stable_iso_data_2012$Sample.ID <- as.numeric(gsub("[^0-9]", "", stable_iso_data_2012$ID))
stable_iso_data_2012$Sample.ID <- as.character(stable_iso_data_2012$Sample.ID)
stable_iso_data_2012$Subspecies <- recode(stable_iso_data_2012$Subspecies, 
                                          "Tridentata" = "T",
                                          "Vaseyana" = "V",
                                          "Wyomingensis" = "W")
levels(stable_iso_data_2012$Subspecies)
stable_iso_data_2012$Subspecies <- as.factor(stable_iso_data_2012$Subspecies)
str(stable_iso_data_2012)

#2021
stable_iso_data_2021<-read.csv("data_csv/Stable_isotope_2021.csv") #74
str(stable_iso_data_2021)
## Remove V-134677 (not sure what subspecies since it says "T/W"). look into later
stable_iso_data_2021 <- stable_iso_data_2021[!(row.names(stable_iso_data_2021) == "31"),] #73
# need to check what subspecies 134677 is 
levels(stable_iso_data_2021$Subspecies)<- list(Tridentata="T", Vaseyana="V", Wyomingensis="W")
levels(stable_iso_data_2021$Subspecies)
stable_iso_data_2021$Subspecies <- as.factor(stable_iso_data_2021$Subspecies)

#renaming the columns to combine the data together
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "Ampl..28"] <- "Ampl28"
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "Ampl..44"] <- "Ampl44"
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "d15N"] <- "Delta15N"
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "d13C"] <- "Delta13C"

#FULL STABLE ISOTOPE DATA
stable_iso_df <- full_join(stable_iso_data_2012,stable_iso_data_2021) #169 of 15 var
stable_iso_df$Subspecies <- droplevels(stable_iso_df$Subspecies)
str(stable_iso_df)

# 2012 BETADISPERSION
# Delta 15 N
D15Ndist <- vegdist(stable_iso_data_2012$Delta15N, method = "bray")
D15Nbetadisper <- betadisper(D15Ndist, group = stable_iso_data_2012$Subspecies)
permutest(D15Nbetadisper) #0.343

D15N12fit<-lm(Delta15N~Subspecies,data = stable_iso_data_2012)
summary(D15N12fit) #sig between subspecies

#Delta 13 C
D13Cdist <- vegdist(stable_iso_data_2012$Delta13C, method = "euclidean")
D13Cbetadisper <- betadisper(D13Cdist, group = stable_iso_data_2012$Subspecies)
permutest(D13Cbetadisper) #0.036
#pairwiseadonis
D13Csubsp.pw.12 <- pairwise.adonis(D13Cdist,as.factor(stable_iso_data_2012$Subspecies))
D13Csubsp.pw.12 #T vs V sig

D13C12fit<-lm(Delta13C~Subspecies,data = stable_iso_data_2012)
summary(D13C12fit)

#2012 VISUALIZATION
# Grouped Scatter plot with marginal density plots
ggscatterhist(
  stable_iso_data_2012, x = "Delta13C", y = "Delta15N", group = "Subspecies",
  color = "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.params = list(fill = "Subspecies", color = "black", size = 0.2)
)

ggscatterhist(
  stable_iso_data_2012, x = "Delta13C", y = "Delta15N", group="Subspecies",
  color = "Subspecies", fill= "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.plot = "boxplot",
  margin.params = list(fill = "Subspecies", color = c("olivedrab", "cadetblue", "goldenrod"), size = 0.2),
  ggtheme = theme_bw()
)

#2021 BETADISPERSION
#DELTA 15 N
D15N21fit<-lm(Delta15N~Subspecies,data = stable_iso_data_2021)
summary(D15N21fit) #p -value= 0.1812
D15N_model <- adonis2(stable_iso_data_2021$Delta15N ~ stable_iso_data_2021$Subspecies)

D15Ndist <- vegdist(stable_iso_data_2021$Delta15N, method = "euclidean")
D15Nbetadisper <- betadisper(D15Ndist, group = stable_iso_data_2021$Subspecies)
permutest(D15Nbetadisper) #0.79
#pairwiseadonis
D15Nsubsp.pw.21 <- pairwise.adonis(D15Ndist,as.factor(stable_iso_data_2021$Subspecies))
D15Nsubsp.pw.21 #T vs V sig

#DELTA 13 C
D13Cdist <- vegdist(stable_iso_data_2021$Delta13C, method = "euclidean")
D13Cbetadisper <- betadisper(D13Cdist, group = stable_iso_data_2021$Subspecies)
permutest(D13Cbetadisper) #0.491
#pairwiseadonis
D13Csubsp.pw.21 <- pairwise.adonis(D13Cdist,as.factor(stable_iso_data_2021$Subspecies))
D13Csubsp.pw.21 #T vs V sig

D15N_model <- adonis2(stable_iso_data_2021$Delta15N ~ stable_iso_data_2021$Subspecies)

D13C21fit<-lm(Delta13C~Subspecies,data = stable_iso_data_2021)
summary(D13C21fit) #p-value = 0.037

D15Nsubsp.pw.21 <- pairwise.adonis(stable_iso_data_2021$Delta13C, stable_iso_data_2021$Subspecies)

#2021 VISUALIZATION
ggscatterhist(
  stable_iso_data_2021, "Delta13C", y = "Delta15N", group="Subspecies",
  color = "Subspecies", fill= "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.plot = "boxplot",
  margin.params = list(fill = "Subspecies", color = c("olivedrab", "cadetblue", "goldenrod"), size = 0.2),
  ggtheme = theme_bw()
)

#FULL VISUALIZATION
ggscatterhist(
  stable_iso_df, "Delta13C", y = "Delta15N", group="Subspecies",
  color = "Subspecies", fill= "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.plot = "boxplot",
  margin.params = list(fill = "Subspecies", color = c("olivedrab", "cadetblue", "goldenrod"), size = 0.2),
  ggtheme = theme_bw()
)

ggscatterhist(
  stable_iso_df, x = "Delta13C", y = "Delta15N", group = "Subspecies",
  color = "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.params = list(fill = "Subspecies", color = "black", size = 0.2),
  ggtheme = theme_bw()
)

D15N_model <- adonis2(stable_iso_df$Delta15N ~ stable_iso_df$Subspecies)

D13C_model <- adonis2(stable_iso_df$Delta13C ~ stable_iso_df$Subspecies, method = "bray")

p1<-ggplot(stable_iso_data_2012,aes(Subspecies,Delta15N))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+ ggtitle("2012")+ theme(legend.position = "none") + ylab("Delta15N")
p2<-ggplot(stable_iso_data_2021,aes(Subspecies,Delta15N))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+ ggtitle("2021")+theme(legend.position = "none")
p3<-ggplot(stable_iso_data_2012,aes(Subspecies,Delta13C))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+ggtitle("2012")+ theme(legend.position = "none")
p4<-ggplot(stable_iso_data_2021,aes(Subspecies,Delta13C))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+ggtitle("2021")+theme(legend.position = "none")
grid.arrange(p1, p2, p3, p4, nrow = 2)

# Clear Global Environment ####
rm(list = ls())
# Thermal analysis#### 

## Phytochemistry manuscript script #
# Set working directory and load necessary packages####
setwd("/Users/ellehorwath/Documents/Orchard_Common_Garden/commongarden")
# List of CRAN packages
cran_packages <- c("dplyr", "effects", "exactRankTests", "factoextra", 
                   "ggcorrplot", "ggfortify", "ggplot2", "glmm", "gridExtra", 
                   "iNEXT", "lme4", "MASS", "nlme", 
                   "picante", "raster", "readr", "reshape2", 
                   "tidyr", "tidyverse", "vegan")

# List of Bioconductor packages
bioc_packages <- c("ANCOMBC","phyloseq", "qiime2R")

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

# Cleaned data read in####
#METADATA
md.OCG <- read.csv("data_csv/metadata_OCG.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #505 obs of 16 variables.
md.OCG <- md.OCG[order(row.names(md.OCG)),]

### Remove duplicates, negatve controls, and MTW.3.7.R_2012
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1', 'ORT.1.5_2012', 'NEG_8-28-21', 'NEG_10-2-20', 'MTW.3.7.R_2012')
md.OCG <- md.OCG[!rownames(md.OCG) %in% rows_to_remove, ] #235 of 21 var

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
OCG_GC <- read.csv("data_csv/OCG_GC_full_clean.csv")#226 obs of 74 variables
OCG_GC <- OCG_GC[order(row.names(OCG_GC)),]
OCG_GC <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(md.OCG)) #217 of 74 variables
md.OCG.GC <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC)) #217 of 21 variables
OCG_GC[is.na(OCG_GC)] <- 0

## 2012 CLEAN GC
OCG_GC_2012 <- read.csv("data_csv/OCG_GC_2012_cleaned.csv", row.names = 1) #157 obs of 74 variables
OCG_GC_2012 <- OCG_GC_2012[order(row.names(OCG_GC_2012)),]
OCG_GC_2012 <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md.OCG)) #147 of 74 variables
md.OCG.GC.2012 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC_2012)) #147
OCG_GC_2012[is.na(OCG_GC_2012)] <- 0

## 2021 CLEAN GC
OCG_GC_2021 <- read.csv("data_csv/OCG_GC_2021.csv", row.names = 1)#70 obs of 74 variables
OCG_GC_2021 <- OCG_GC_2021[order(row.names(OCG_GC_2021)),] 
OCG_GC_2021 <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md.OCG)) #70 of 74 variables
md.OCG.GC.2021 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC_2021)) #70
OCG_GC_2021[is.na(OCG_GC_2021)] <- 0

## CLEAN LCMS 
OCG_LCMS_3uL <- read.csv("data_csv/OCG_LCMS_3uL_cleaned.csv", row.names = 1) #112
OCG_LCMS_3uL <- OCG_LCMS_3uL[order(row.names(OCG_LCMS_3uL)),]
OCG_LCMS_3uL <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.OCG)) #111
md.OCG.LCMS.3 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_LCMS_3uL)) #111 
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
     xlab="PC 1 (19.81%)", ylab="PC 2 (12.46%)", 
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
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$cluster_assignments, show.groups = "1", col = "lightcoral")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$cluster_assignments, show.groups = "2", col = "rosybrown")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$cluster_assignments, show.groups = "3", col = "darkseagreen")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$cluster_assignments, show.groups = "4", col = "peachpuff")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$cluster_assignments, show.groups = "5", col = "darkturquoise")


# #SPECTRAL CLUSTERING 2012 GC 
# set.seed(75)
# scgc12 <- specc(pca_scores_GC12, centers = 5)
# scgc12
# centers(scgc12)
# size(scgc12)
# withinss(scgc12)
# plot(pca_scores_GC12, col = scgc12, main = "PCA of 2012 GC data with spectral clustering")

## 2021 GC PCA ####
data.pca_2021 <- princomp(data_normalized_2021)
summary(data.pca_2021)
data.pca_2021$loadings[, 1:2]

fviz_eig(data.pca_2021, addlabels = TRUE) #(20.5, 13.3)

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
set.seed(64)
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
     xlab="PC 1 (20.54%)", ylab="PC 2 (13.33%)", 
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
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$cluster_assignments, show.groups = "1", col = "lightcoral")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$cluster_assignments, show.groups = "2", col = "rosybrown")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$cluster_assignments, show.groups = "3", col = "darkseagreen")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$cluster_assignments, show.groups = "4", col = "peachpuff")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$cluster_assignments, show.groups = "5", col = "darkturquoise")

# #SPECTRAL CLUSTERING 2021 GC 
# set.seed(24)
# scgc21 <- specc(pca_scores_GC21, centers = 4)
# scgc21
# centers(scgc21)
# size(scgc21)
# withinss(scgc21)
# plot(pca_scores_GC21, col = scgc21, main = "PCA of 2021 GC data with spectral clustering")

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
k <- 5

# Perform k-means clustering on the PCA scores
kmeans_result <- kmeans(pca_scores_LCMS, centers = k)

# Get cluster assignments for each sample
cluster_assignments <- kmeans_result$cluster

# Visualize the clusters (optional)
plot(pca_scores_LCMS, col = cluster_assignments, main = "PCA of LCMS data with k-means Clustering")

## Adding k means clusters to md 
md.OCG.LCMS.3$cluster_assignments <- cluster_assignments

plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1 (13.51%)", ylab="PC 2 (9.51%)", 
     main="PCA of LCMS data with k-means clusters", 
     pch = 19,
     col= c("lightcoral","rosybrown",'darkseagreen','peachpuff','darkturquoise')[md.OCG.LCMS.3$cluster_assignments],
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("1","2","3","4","5"),
       col= c("lightcoral","rosybrown",'darkseagreen','peachpuff','darkturquoise'),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$cluster_assignments, show.groups = "1", col = "lightcoral")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$cluster_assignments, show.groups = "2", col = "rosybrown")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$cluster_assignments, show.groups = "3", col = "darkseagreen")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$cluster_assignments, show.groups = "4", col = "peachpuff")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$cluster_assignments, show.groups = "5", col = "darkturquoise")

# ##DBSCAN ####
# # Perform DBSCAN clustering
# dbscan::kNNdistplot(OCG_LCMS_3uL_subset, k =  5)
# abline(h = 6.0e+06, lty = 2)
# set.seed(123)
# # fpc package
# res.fpc <- fpc::dbscan(OCG_LCMS_3uL_subset, eps = 6.0e+06 , MinPts = 2)
# plot(res.fpc, pca_scores_LCMS)
# # dbscan package
# res.db <- dbscan::dbscan(pca_scores_LCMS, 14.5, 2)
# all(res.fpc$cluster == res.db$cluster) #TRUE
# 
# fviz_cluster(res.fpc, pca_scores_LCMS, geom = "point")
# 
# set.seed(1)
# dbscan_result <- dbscan(pca_scores_LCMS, eps = 14.5, 2)
# dbscan_result$cluster
# fviz_cluster(dbscan_result, data = pca_scores_LCMS, geom = "point", 
#              outlier.pointsize = 1, main = "DBSCAN Clustering on LCMS PCA-transformed Data", palette = "RdYlGn", shape = 19)+ theme_classic() 
# 
# #HDBSCAN
# hdbscan_result <- hdbscan(pca_scores_LCMS, minPts = 5)
# plot(hdbscan_result, col = hdbscan_result$cluster+1, pch = 20)
# 
# #SPECTRAL CLUSTERING
# set.seed(43)
# sc <- specc(pca_scores_LCMS, centers = 5)
# sc
# centers(sc)
# size(sc)
# withinss(sc)
# plot(pca_scores_LCMS, col = sc, main = "PCA of LCMS data with spectral clustering")

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

# ANCOM with Chem data ####
## Create ANCOM Function ####

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
  
  ######################## Function 
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  ######################## Function 
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

#FULL GC ANCOM ####
##Remove entries with insufficient areas
OCG_GC_subset[OCG_GC_subset < 10] <- 0
OCG_GC_subset.a <- OCG_GC_subset[rowSums(OCG_GC_subset) > 0,] # each observation needs at least 10

summary(rowSums(OCG_GC_subset.a)) #179
summary(colSums(OCG_GC_subset.a)) #2500

OCG_GC_subset.a <- OCG_GC_subset.a[,colSums(OCG_GC_subset.a) > 10] # each sample needs at least 10 

md.OCG.GC.a <- subset(md.OCG.GC, row.names(md.OCG.GC) %in% row.names(OCG_GC_subset.a)) 

OCG_GC_subset.a.t <- t(OCG_GC_subset.a) 
OCG_GC_subset.a_t <- OCG_GC_subset.a.t[, colnames(OCG_GC_subset.a.t) %in% row.names(md.OCG.GC.a), drop = FALSE]
OCG_GC_subset.a <- t(OCG_GC_subset.a_t)
OCG_GC_subset.a <- as.data.frame(OCG_GC_subset.a) 

OCG_GC_subset.a[OCG_GC_subset.a < 10] <- 0  # repeat cleaning after trimming
OCG_GC_subset.a <- OCG_GC_subset.a[rowSums(OCG_GC_subset.a) > 0,] # each observation needs at least 10

summary(rowSums(OCG_GC_subset.a)) #179
summary(colSums(OCG_GC_subset.a)) #2500

OCG_GC_subset.a <- OCG_GC_subset.a[,colSums(OCG_GC_subset.a) > 10]

#ANCOM requires that data be formatted so that first *column* is named "Sample.ID"
md.OCG.GC_sbst <- data.frame("Sample.ID" = row.names(md.OCG.GC.a), md.OCG.GC.a)
OCG_GC_sbst <- data.frame("Sample.ID" = row.names(OCG_GC_subset.a), OCG_GC_subset.a, check.names = F)
row.names(OCG_GC_sbst) == row.names(md.OCG.GC_sbst) #TRUE

### ANCOM SUBSPECIES GC ####
ANCOM_subspecies_GC <- ANCOM.main(OCG_GC_sbst,md.OCG.GC_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)

#Create objects of significant ASVs
sigGC_subspecies <- subset(ANCOM_subspecies_GC$W.taxa, ANCOM_subspecies_GC$W.taxa$W_stat > 0)[,1]
sigGC_subspecies <- as.data.frame(sigGC_subspecies) 
row.names(sigGC_subspecies) <- sigGC_subspecies[54:1,1] 

sigGC_subspecies[,1] <- c(54:1) 
sigGC_subspecies_t <- t(sigGC_subspecies) 
sigGC_subspecies_t <- as.data.frame(sigGC_subspecies_t)
colnames(sigGC_subspecies_t) <- as.character(colnames(sigGC_subspecies_t))
print(colnames(sigGC_subspecies_t))   
sigGC_subspecies_t <- sigGC_subspecies_t[,order(colnames(sigGC_subspecies_t))]
rownames(sigGC_subspecies_t) <- c("sig_rank")

GC.OCG_sig_subspecies <-  t(subset(t(OCG_GC_sbst), colnames(OCG_GC_sbst) %in% row.names(sigGC_subspecies)))
GC.OCG_sig_subspecies <- GC.OCG_sig_subspecies[,order(colnames(GC.OCG_sig_subspecies))]
colnames(sigGC_subspecies_t) == colnames(GC.OCG_sig_subspecies) #TRUE

GC.OCG_sig_subspecies <- rbind(GC.OCG_sig_subspecies, sigGC_subspecies_t)
GC.OCG_sig_subspecies_t <- as.data.frame(t(GC.OCG_sig_subspecies))
GC.OCG_sig_subspecies_t$sig_rank <- as.numeric(GC.OCG_sig_subspecies_t$sig_rank) 
GC.OCG_sig_subspecies_t <- GC.OCG_sig_subspecies_t[order(GC.OCG_sig_subspecies_t$sig_rank),] 
GC.OCG_sig_subspecies_t <- subset(GC.OCG_sig_subspecies_t, select=-c(sig_rank))
GC.OCG_sig_subspecies_t <- as.data.frame(t(GC.OCG_sig_subspecies_t))

#Build objects for plotting
sbsplotGC <- data.frame(GC.OCG_sig_subspecies_t[,1:10], "subspecies" = md.OCG.GC_sbst$Subspecies, check.names = FALSE)
sbsplotGC[,1:10] <- lapply(sbsplotGC[,1:10], function(x) as.numeric(as.character(x)))
sbsplotGC[,1:10] <- log(sbsplotGC[,1:10]+1)
sbsplotGC_sub <- data.frame(sample=rownames(sbsplotGC),sbsplotGC, check.names = F)
sbsplotGC_sublong <- melt(sbsplotGC_sub)

ANCOM_subspecies_GC$W.taxa

sbsplotGC_sublong$subspecies <-  factor(sbsplotGC_sublong$subspecies, levels = c("T", "V", "W"))

##### SUBSPECIES FIGURE TOP ANCOM COMPOUNDS ALL GC ####
sp1 <- ggplot(sbsplotGC_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("Full GC ANCOM for Subspecies")+ 
  theme_classic()

## YEAR Full GC ####
#Run ANCOM, specify variable
ANCOM_yr.GC <- ANCOM.main(OCG_GC_sbst,md.OCG.GC_sbst,F,F,"Year",NULL,NULL,F,NULL,2,.05,.9)

sigGC_yr <- subset(ANCOM_yr.GC$W.taxa, ANCOM_yr.GC$W.taxa$W_stat > 0)[,1]
sigGC_yr <- as.data.frame(sigGC_yr)
row.names(sigGC_yr) <- sigGC_yr[54:1,1]
sigGC_yr[,1] <- c(54:1)

sigGC_yr_t <- t(sigGC_yr)
sigGC_yr_t <- as.data.frame(sigGC_yr_t)
colnames(sigGC_yr_t) <- as.character(colnames(sigGC_yr_t))
print(colnames(sigGC_yr_t))
sigGC_yr_t <- sigGC_yr_t[,order(colnames(sigGC_yr_t))]
rownames(sigGC_yr_t) <- c("sig_rank")

GC_sig_yr <-  t(subset(t(OCG_GC_sbst), colnames(OCG_GC_sbst) %in% row.names(sigGC_yr)))

GC_sig_yr <- GC_sig_yr[,order(colnames(GC_sig_yr))]
colnames(sigGC_yr_t) == colnames(GC_sig_yr) #sanity check:TRUE

GC_sig_yr <- rbind(GC_sig_yr, sigGC_yr_t)
GC_sig_yr_t <- as.data.frame(t(GC_sig_yr))
GC_sig_yr_t$sig_rank <- as.numeric(GC_sig_yr_t$sig_rank) 
GC_sig_yr_t <- GC_sig_yr_t[order(GC_sig_yr_t$sig_rank),] 
GC_sig_yr_t <- subset(GC_sig_yr_t, select=-c(sig_rank))
GC_sig_yr_t <- as.data.frame(t(GC_sig_yr_t))

#Build objects for plotting
plotGC.yr <- data.frame(GC_sig_yr_t[,1:10], "year" = md.OCG.GC_sbst$Year, check.names = FALSE)
plotGC.yr[,1:10] <- lapply(plotGC.yr[,1:10], function(x) as.numeric(as.character(x)))
plotGC.yr[,1:10] <- log(plotGC.yr[,1:10]+1)
plotGC.yr_sub <- data.frame(sample=rownames(plotGC.yr),plotGC.yr, check.names = F)
plotGC.yr_sublong <- melt(plotGC.yr_sub)

plotGC.yr$year <-  factor(plotGC.yr$year, levels = c("2012","2021"))

### FIGURE FOR YEAR ####
ggplot(plotGC.yr_sublong, aes(y = value, x = year, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("GC ANCOM for year")+
  labs(color = "Compounds") +
  theme_classic()

# 2012 GC ANCOM ####
summary(rowSums(OCG_GC_2012_subset)) #7962
summary(colSums(OCG_GC_2012_subset)) #2358

OCG_GC_2012_subset.a <- OCG_GC_2012_subset[,colSums(OCG_GC_2012_subset) > 10] 
OCG_GC_2012_subset.a <- OCG_GC_2012_subset.a[rowSums(OCG_GC_2012_subset.a) > 0,] 
md.OCG.GC.2012.a <- subset(md.OCG.GC.2012, row.names(md.OCG.GC.2012) %in% row.names(OCG_GC_2012_subset.a)) 

OCG_GC_2012_subset.a.t <- t(OCG_GC_2012_subset.a) 
OCG_GC_2012_subset.a.t <- OCG_GC_2012_subset.a.t[, colnames(OCG_GC_2012_subset.a.t) %in% row.names(md.OCG.GC.2012.a), drop = FALSE]
OCG_GC_2012_subset.a <- t(OCG_GC_2012_subset.a.t)
OCG_GC_2012_subset.a <- as.data.frame(OCG_GC_2012_subset.a) 

OCG_GC_2012_subset.a[OCG_GC_2012_subset.a < 10] <- 0  
OCG_GC_2012_subset.a <- OCG_GC_2012_subset.a[rowSums(OCG_GC_2012_subset.a) > 0,] 

summary(rowSums(OCG_GC_2012_subset.a)) #7962
summary(colSums(OCG_GC_2012_subset.a)) #2358

OCG_GC_2012_subset.a <- OCG_GC_2012_subset.a[,colSums(OCG_GC_2012_subset.a) > 10]

#ANCOM requires that data be formatted so that first *column* is named "Sample.ID"
md.OCG.GC.2012_sbst <- data.frame("Sample.ID" = row.names(md.OCG.GC.2012.a), md.OCG.GC.2012.a)
OCG_GC_2012_sbst <- data.frame("Sample.ID" = row.names(OCG_GC_2012_subset.a), OCG_GC_2012_subset.a, check.names = F)
row.names(OCG_GC_2012_sbst) == row.names(md.OCG.GC.2012_sbst) #TRUE

## SUBSPECIES 2012 GC ####
ANCOM_subspecies_GC.2012 <- ANCOM.main(OCG_GC_2012_sbst,md.OCG.GC.2012_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)

#Create objects of significant ASVs
sigGC.2012_subspecies <- subset(ANCOM_subspecies_GC.2012$W.taxa, ANCOM_subspecies_GC.2012$W.taxa$W_stat > 0)[,1]
sigGC.2012_subspecies <- as.data.frame(sigGC.2012_subspecies) 
row.names(sigGC.2012_subspecies) <- sigGC.2012_subspecies[47:1,1] 

sigGC.2012_subspecies[,1] <- c(47:1) 
sigGC.2012_subspecies_t <- t(sigGC.2012_subspecies) 
sigGC.2012_subspecies_t <- as.data.frame(sigGC.2012_subspecies_t)
colnames(sigGC.2012_subspecies_t) <- as.character(colnames(sigGC.2012_subspecies_t))
print(colnames(sigGC.2012_subspecies_t))

sigGC.2012_subspecies_t <- sigGC.2012_subspecies_t[,order(colnames(sigGC.2012_subspecies_t))]
rownames(sigGC.2012_subspecies_t) <- c("sig_rank")

GC.2012.OCG_sigsbst_subspecies <-  t(subset(t(OCG_GC_2012_sbst), colnames(OCG_GC_2012_sbst) %in% row.names(sigGC.2012_subspecies)))

GC.2012.OCG_sigsbst_subspecies <- GC.2012.OCG_sigsbst_subspecies[,order(colnames(GC.2012.OCG_sigsbst_subspecies))]
colnames(sigGC.2012_subspecies_t) == colnames(GC.2012.OCG_sigsbst_subspecies) #sanity check:TRUE

GC.2012.OCG_sigsbst_subspecies <- rbind(GC.2012.OCG_sigsbst_subspecies, sigGC.2012_subspecies_t)
GC.2012.OCG_sigsbst_subspecies_t <- as.data.frame(t(GC.2012.OCG_sigsbst_subspecies))
GC.2012.OCG_sigsbst_subspecies_t$sig_rank <- as.numeric(GC.2012.OCG_sigsbst_subspecies_t$sig_rank) 
GC.2012.OCG_sigsbst_subspecies_t <- GC.2012.OCG_sigsbst_subspecies_t[order(GC.2012.OCG_sigsbst_subspecies_t$sig_rank),] 
GC.2012.OCG_sigsbst_subspecies_t <- subset(GC.2012.OCG_sigsbst_subspecies_t, select=-c(sig_rank))
GC.2012.OCG_sigsbst_subspecies_t <- as.data.frame(t(GC.2012.OCG_sigsbst_subspecies_t))

#Build objects for plotting
sbsplotGC2012 <- data.frame(GC.2012.OCG_sigsbst_subspecies_t[,1:10], "subspecies" = md.OCG.GC.2012_sbst$Subspecies, check.names = FALSE)
sbsplotGC2012[,1:10] <- lapply(sbsplotGC2012[,1:10], function(x) as.numeric(as.character(x)))
sbsplotGC2012[,1:10] <- log(sbsplotGC2012[,1:10]+1)
sbsplotGC2012_sub <- data.frame(sample=rownames(sbsplotGC2012),sbsplotGC2012, check.names = F)
sbsplotGC2012_sublong <- melt(sbsplotGC2012_sub)

ANCOM_subspecies_GC.2012$W.taxa

sbsplotGC2012_sublong$subspecies <-  factor(sbsplotGC2012_sublong$subspecies, levels = c("T", "V", "W"))

#### SUBSPECIES FIGURE TOP ANCOM COMPOUNDS 2012 GC ####
sp2 <- ggplot(sbsplotGC2012_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("2012 GC ANCOM for Subspecies")+
  labs(color = "Compounds") +
  theme_classic()

## PLOIDY 2012 GC ####
#Run ANCOM, specify variable
ANCOM_ploidy_GC <- ANCOM.main(OCG_GC_2012_sbst,md.OCG.GC.2012_sbst,F,F,"Ploidy",NULL,NULL,F,NULL,2,.05,.9)
#Create objects of significant ASVs
sigGC.2012_ploidy <- subset(ANCOM_ploidy_GC$W.taxa, ANCOM_ploidy_GC$W.taxa$W_stat > 0)[,1]
sigGC.2012_ploidy <- as.data.frame(sigGC_ploidy)
row.names(sigGC.2012_ploidy) <- sigGC.2012_ploidy[47:1,1] 

sigGC.2012_ploidy[,1] <- c(47:1) 
sigGC.2012_ploidy_t <- t(sigGC.2012_ploidy) 
sigGC.2012_ploidy_t <- as.data.frame(sigGC.2012_ploidy_t)
colnames(sigGC.2012_ploidy_t) <- as.character(colnames(sigGC.2012_ploidy_t))
print(colnames(sigGC.2012_ploidy_t))

sigGC.2012_ploidy_t <- sigGC.2012_ploidy_t[,order(colnames(sigGC.2012_ploidy_t))]
rownames(sigGC.2012_ploidy_t) <- c("sig_rank")

GC.2012.OCG_sigsbst_ploidy <-  t(subset(t(OCG_GC_2012_sbst), colnames(OCG_GC_2012_sbst) %in% row.names(sigGC.2012_ploidy)))

GC.2012.OCG_sigsbst_ploidy <- GC.2012.OCG_sigsbst_ploidy[,order(colnames(GC.2012.OCG_sigsbst_ploidy))]
colnames(sigGC.2012_subspecies_t) == colnames(GC.2012.OCG_sigsbst_ploidy) #sanity check:TRUE

GC.2012.OCG_sigsbst_ploidy <- rbind(GC.2012.OCG_sigsbst_ploidy, sigGC.2012_ploidy_t)
GC.2012.OCG_sigsbst_ploidy_t <- as.data.frame(t(GC.2012.OCG_sigsbst_ploidy))
GC.2012.OCG_sigsbst_ploidy_t$sig_rank <- as.numeric(GC.2012.OCG_sigsbst_ploidy_t$sig_rank) 
GC.2012.OCG_sigsbst_ploidy_t <- GC.2012.OCG_sigsbst_ploidy_t[order(GC.2012.OCG_sigsbst_ploidy_t$sig_rank),] 
GC.2012.OCG_sigsbst_ploidy_t <- subset(GC.2012.OCG_sigsbst_ploidy_t, select=-c(sig_rank))
GC.2012.OCG_sigsbst_ploidy_t <- as.data.frame(t(GC.2012.OCG_sigsbst_ploidy_t))

#Build objects for plotting
plplotGC2012 <- data.frame(GC.2012.OCG_sigsbst_ploidy_t[,1:10], "ploidy" = md.OCG.GC.2012_sbst$Ploidy, check.names = FALSE)
plplotGC2012[,1:10] <- lapply(plplotGC2012[,1:10], function(x) as.numeric(as.character(x)))
plplotGC2012[,1:10] <- log(plplotGC2012[,1:10]+1)
plplotGC2012_sub <- data.frame(sample=rownames(plplotGC2012),plplotGC2012, check.names = F)
plplotGC2012_sublong <- melt(plplotGC2012_sub)

ANCOM_subspecies_GC.2012$W.taxa

plplotGC2012_sublong$ploidy<-  factor(plplotGC2012_sublong$ploidy, levels = c("2n", "4n"))

#### PLOIDY FIGURE TOP ANCOM COMPOUNDS 2012 GC ####
pl1 <- ggplot(sbsplotGC2012.ploidy_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("2012 GC ANCOM for Ploidy")+
  labs(color = "Compounds") +
  theme_classic()


## SUBSPECIES PLOIDY 2012 GC ####
#Run ANCOM, specify variable
ANCOM_subsp_ploidy.GC12 <- ANCOM.main(OCG_GC_2012_sbst, md.OCG.GC.2012_sbst,F,F,"Subsp_ploidy",NULL,NULL,F,NULL,2,.05,.9)

sigGC.2012_subsp_ploidy <- subset(ANCOM_subsp_ploidy.GC12$W.taxa, ANCOM_subsp_ploidy.GC12$W.taxa$W_stat > 0)[,1]
sigGC.2012_subsp_ploidy <- as.data.frame(sigGC.2012_subsp_ploidy)
row.names(sigGC.2012_subsp_ploidy) <- sigGC.2012_subsp_ploidy[47:1,1]
sigGC.2012_subsp_ploidy[,1] <- c(47:1)

sigGC.2012_subsp_ploidy_t <- t(sigGC.2012_subsp_ploidy)
sigGC.2012_subsp_ploidy_t <- as.data.frame(sigGC.2012_subsp_ploidy_t)
colnames(sigGC.2012_subsp_ploidy_t) <- as.character(colnames(sigGC.2012_subsp_ploidy_t))
print(colnames(sigGC.2012_subsp_ploidy_t))
sigGC.2012_subsp_ploidy_t <- sigGC.2012_subsp_ploidy_t[,order(colnames(sigGC.2012_subsp_ploidy_t))]
rownames(sigGC.2012_subsp_ploidy_t) <- c("sig_rank")

#write.csv(ANCOM_subsp_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_subsp_ploidy.csv")

GC.2012.OCG_sigsbst_subspploidy <-  t(subset(t(OCG_GC_2012_sbst), colnames(OCG_GC_2012_sbst) %in% row.names(sigGC.2012_subsp_ploidy)))

GC.2012.OCG_sigsbst_subspploidy <- GC.2012.OCG_sigsbst_subspploidy[,order(colnames(GC.2012.OCG_sigsbst_subspploidy))]
colnames(sigGC.2012_subsp_ploidy_t) == colnames(GC.2012.OCG_sigsbst_subspploidy) #sanity check:TRUE

GC.2012.OCG_sigsbst_subspploidy <- rbind(GC.2012.OCG_sigsbst_subspploidy, sigGC.2012_subsp_ploidy_t)
GC.2012.OCG_sigsbst_subspploidy_t <- as.data.frame(t(GC.2012.OCG_sigsbst_subspploidy))
GC.2012.OCG_sigsbst_subspploidy_t$sig_rank <- as.numeric(GC.2012.OCG_sigsbst_subspploidy_t$sig_rank) 
GC.2012.OCG_sigsbst_subspploidy_t <- GC.2012.OCG_sigsbst_subspploidy_t[order(GC.2012.OCG_sigsbst_subspploidy_t$sig_rank),] 
GC.2012.OCG_sigsbst_subspploidy_t <- subset(GC.2012.OCG_sigsbst_subspploidy_t, select=-c(sig_rank))
GC.2012.OCG_sigsbst_subspploidy_t <- as.data.frame(t(GC.2012.OCG_sigsbst_subspploidy_t))

#Build objects for plotting
sbsp.plplotGC2012 <- data.frame(GC.2012.OCG_sigsbst_subspploidy_t[,1:10], "subsp_ploidy" = md.OCG.GC.2012_sbst$Subsp_ploidy, check.names = FALSE)
sbsp.plplotGC2012[,1:10] <- lapply(sbsp.plplotGC2012[,1:10], function(x) as.numeric(as.character(x)))
sbsp.plplotGC2012[,1:10] <- log(sbsp.plplotGC2012[,1:10]+1)
sbsp.plplotGC2012_sub <- data.frame(sample=rownames(sbsp.plplotGC2012),sbsp.plplotGC2012, check.names = F)
sbsp.plplotGC2012_sublong <- melt(sbsp.plplotGC2012_sub)

ANCOM_subsp_ploidy.GC12$W.taxa

sbsp.plplotGC2012$subsp_ploidy <-  factor(sbsp.plplotGC2012$subsp_ploidy, levels = c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))

#### FIGURE FOR SUBSPECIES PLOIDY ####
spl1 <- ggplot(sbsp.plplotGC2012_sublong, aes(y = value, x = subsp_ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("2012 GC ANCOM for Subspecies and Ploidy")+
  labs(color = "Compounds") +
  theme_classic()


# 2021 GC ANCOM ####
summary(rowSums(OCG_GC_2021_subset)) #179
summary(colSums(OCG_GC_2021_subset)) #842.5

OCG_GC_2021_subset.a <- OCG_GC_2021_subset[,colSums(OCG_GC_2021_subset) > 10] 
OCG_GC_2021_subset.a <- OCG_GC_2021_subset.a[rowSums(OCG_GC_2021_subset.a) > 0,] 
md.OCG.GC.2021.a <- subset(md.OCG.GC.2021, row.names(md.OCG.GC.2021) %in% row.names(OCG_GC_2021_subset.a)) 

OCG_GC_2021_subset.a.t <- t(OCG_GC_2021_subset.a) 
OCG_GC_2021_subset.a_t <- OCG_GC_2021_subset.a.t[, colnames(OCG_GC_2021_subset.a.t) %in% row.names(md.OCG.GC.2021.a), drop = FALSE]
OCG_GC_2021_subset.a <- t(OCG_GC_2021_subset.a_t)
OCG_GC_2021_subset.a <- as.data.frame(OCG_GC_2021_subset.a) 

OCG_GC_2021_subset.a[OCG_GC_2021_subset.a < 10] <- 0  
OCG_GC_2021_subset.a <- OCG_GC_2021_subset.a[rowSums(OCG_GC_2021_subset.a) > 0,] 

summary(rowSums(OCG_GC_2021_subset.a)) #179
summary(colSums(OCG_GC_2021_subset.a)) #842.5

OCG_GC_2021_subset.a <- OCG_GC_2021_subset.a[,colSums(OCG_GC_2021_subset.a) > 10]

md.OCG.GC.2021_sbst <- data.frame("Sample.ID" = row.names(md.OCG.GC.2021.a), md.OCG.GC.2021.a)
OCG_GC_2021_sbst <- data.frame("Sample.ID" = row.names(OCG_GC_2021_subset.a), OCG_GC_2021_subset.a, check.names = F)
row.names(OCG_GC_2021_sbst) == row.names(md.OCG.GC.2021_sbst) #TRUE

## SUBSPECIES 2021 GC ####
ANCOM_subspecies_GC.2021 <- ANCOM.main(OCG_GC_2021_sbst,md.OCG.GC.2021_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)

#Create objects of significant compounds
sigGC.2021_subspecies <- subset(ANCOM_subspecies_GC.2021$W.taxa, ANCOM_subspecies_GC.2021$W.taxa$W_stat > 0)[,1]
sigGC.2021_subspecies <- as.data.frame(sigGC.2021_subspecies) 
row.names(sigGC.2021_subspecies) <- sigGC.2021_subspecies[29:1,1] 

sigGC.2021_subspecies[,1] <- c(29:1) 
sigGC.2021_subspecies_t <- t(sigGC.2021_subspecies) 
sigGC.2021_subspecies_t <- as.data.frame(sigGC.2021_subspecies_t)
colnames(sigGC.2021_subspecies_t) <- as.character(colnames(sigGC.2021_subspecies_t))
print(colnames(sigGC.2021_subspecies_t))

sigGC.2021_subspecies_t <- sigGC.2021_subspecies_t[,order(colnames(sigGC.2021_subspecies_t))]
rownames(sigGC.2021_subspecies_t) <- c("sig_rank")

GC.2021.OCG_sigsbst_subspecies <-  t(subset(t(OCG_GC_2021_sbst), colnames(OCG_GC_2021_sbst) %in% row.names(sigGC.2021_subspecies)))

GC.2021.OCG_sigsbst_subspecies <- GC.2021.OCG_sigsbst_subspecies[,order(colnames(GC.2021.OCG_sigsbst_subspecies))]
colnames(sigGC.2021_subspecies_t) == colnames(GC.2021.OCG_sigsbst_subspecies) #sanity check:TRUE

GC.2021.OCG_sigsbst_subspecies <- rbind(GC.2021.OCG_sigsbst_subspecies, sigGC.2021_subspecies_t)
GC.2021.OCG_sigsbst_subspecies_t <- as.data.frame(t(GC.2021.OCG_sigsbst_subspecies))
GC.2021.OCG_sigsbst_subspecies_t$sig_rank <- as.numeric(GC.2021.OCG_sigsbst_subspecies_t$sig_rank) 
GC.2021.OCG_sigsbst_subspecies_t <- GC.2021.OCG_sigsbst_subspecies_t[order(GC.2021.OCG_sigsbst_subspecies_t$sig_rank),] 
GC.2021.OCG_sigsbst_subspecies_t <- subset(GC.2021.OCG_sigsbst_subspecies_t, select=-c(sig_rank))
GC.2021.OCG_sigsbst_subspecies_t <- as.data.frame(t(GC.2021.OCG_sigsbst_subspecies_t))

#Build objects for plotting
sbsplotGC2021 <- data.frame(GC.2021.OCG_sigsbst_subspecies_t[,1:10], "subspecies" = md.OCG.GC.2021_sbst$Subspecies, check.names = FALSE)
sbsplotGC2021[,1:10] <- lapply(sbsplotGC2021[,1:10], function(x) as.numeric(as.character(x)))
sbsplotGC2021[,1:10] <- log(sbsplotGC2021[,1:10]+1)
sbsplotGC2021_sub <- data.frame(sample=rownames(sbsplotGC2021),sbsplotGC2021, check.names = F)
sbsplotGC2021_sublong <- melt(sbsplotGC2021_sub)

ANCOM_subspecies_GC.2021$W.taxa

sbsplotGC2021_sublong$subspecies <-  factor(sbsplotGC2021_sublong$subspecies, levels = c("T", "V", "W"))

#### SUBSPECIES FIGURE TOP ANCOM COMPOUNDS 2021 GC ####
sp3 <- ggplot(sbsplotGC2021_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("2021 GC ANCOM for Subspecies")+theme_classic()

## PLOIDY 2021 GC ####
#Run ANCOM, specify variable
ANCOM_ploidy_GC.21 <- ANCOM.main(OCG_GC_2021_sbst,md.OCG.GC.2021_sbst,F,F,"Ploidy",NULL,NULL,F,NULL,2,.05,.9)

sigGC.2021_ploidy <- subset(ANCOM_ploidy_GC.21$W.taxa, ANCOM_ploidy_GC.21$W.taxa$W_stat > 0)[,1]
sigGC.2021_ploidy <- as.data.frame(sigGC.2021_ploidy)
row.names(sigGC.2021_ploidy) <- sigGC.2021_ploidy[37:1,1] 

sigGC.2021_ploidy[,1] <- c(37:1) 
sigGC.2021_ploidy_t <- t(sigGC.2021_ploidy) 
sigGC.2021_ploidy_t <- as.data.frame(sigGC.2021_ploidy_t)
colnames(sigGC.2021_ploidy_t) <- as.character(colnames(sigGC.2021_ploidy_t))
print(colnames(sigGC.2021_ploidy_t))

sigGC.2021_ploidy_t <- sigGC.2021_ploidy_t[,order(colnames(sigGC.2021_ploidy_t))]
rownames(sigGC.2021_ploidy_t) <- c("sig_rank")

GC.2021.OCG_sigsbst_ploidy <-  t(subset(t(OCG_GC_2021_sbst), colnames(OCG_GC_2021_sbst) %in% row.names(sigGC.2021_ploidy)))

GC.2021.OCG_sigsbst_ploidy <- GC.2021.OCG_sigsbst_ploidy[,order(colnames(GC.2021.OCG_sigsbst_ploidy))]
colnames(sigGC.2021_ploidy_t) == colnames(GC.2021.OCG_sigsbst_ploidy) #sanity check:TRUE

GC.2021.OCG_sigsbst_ploidy <- rbind(GC.2021.OCG_sigsbst_ploidy, sigGC.2021_ploidy_t)
GC.2021.OCG_sigsbst_ploidy_t <- as.data.frame(t(GC.2021.OCG_sigsbst_ploidy))
GC.2021.OCG_sigsbst_ploidy_t$sig_rank <- as.numeric(GC.2021.OCG_sigsbst_ploidy_t$sig_rank) 
GC.2021.OCG_sigsbst_ploidy_t <- GC.2021.OCG_sigsbst_ploidy_t[order(GC.2021.OCG_sigsbst_ploidy_t$sig_rank),] 
GC.2021.OCG_sigsbst_ploidy_t <- subset(GC.2021.OCG_sigsbst_ploidy_t, select=-c(sig_rank))
GC.2021.OCG_sigsbst_ploidy_t <- as.data.frame(t(GC.2021.OCG_sigsbst_ploidy_t))

#Build objects for plotting
plplotGC2021.ploidy <- data.frame(GC.2021.OCG_sigsbst_ploidy_t[,1:10], "ploidy" = md.OCG.GC.2021_sbst$Ploidy, check.names = FALSE)
plplotGC2021.ploidy[,1:10] <- lapply(plplotGC2021.ploidy[,1:10], function(x) as.numeric(as.character(x)))
plplotGC2021.ploidy[,1:10] <- log(plplotGC2021.ploidy[,1:10]+1)
plplotGC2021.ploidy_sub <- data.frame(sample=rownames(plplotGC2021.ploidy),plplotGC2021.ploidy, check.names = F)
plplotGC2021.ploidy_sublong <- melt(plplotGC2021.ploidy_sub)

ANCOM_subspecies_GC.2021$W.taxa

plplotGC2021.ploidy_sublong$ploidy<-  factor(plplotGC2021.ploidy_sublong$ploidy, levels = c("2n", "4n"))

#### PLOIDY FIGURE TOP ANCOM COMPOUNDS 2021 GC ####
pl2 <- ggplot(sbsplotGC2021.ploidy_sublong, aes(y = value, x = ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Ploidy") + 
  ggtitle("2021 GC ANCOM for Ploidy")+
  labs(color = "Compounds") +
  theme_classic()

## SIGNIFICANCE BY SUBSPECIES PLOIDY 2021 GC ####
#Run ANCOM, specify variable
ANCOM_subsp_ploidy.GC21 <- ANCOM.main(OCG_GC_2021_sbst, md.OCG.GC.2021_sbst,F,F,"Subsp_ploidy",NULL,NULL,F,NULL,2,.05,.9)

sigGC.2021_subsp_ploidy <- subset(ANCOM_subsp_ploidy.GC21$W.taxa, ANCOM_subsp_ploidy.GC21$W.taxa$W_stat > 0)[,1]
sigGC.2021_subsp_ploidy <- as.data.frame(sigGC.2021_subsp_ploidy)
row.names(sigGC.2021_subsp_ploidy) <- sigGC.2021_subsp_ploidy[36:1,1]
sigGC.2021_subsp_ploidy[,1] <- c(36:1)

sigGC.2021_subsp_ploidy_t <- t(sigGC.2021_subsp_ploidy)
sigGC.2021_subsp_ploidy_t <- as.data.frame(sigGC.2021_subsp_ploidy_t)
colnames(sigGC.2021_subsp_ploidy_t) <- as.character(colnames(sigGC.2021_subsp_ploidy_t))
print(colnames(sigGC.2021_subsp_ploidy_t))
sigGC.2021_subsp_ploidy_t <- sigGC.2021_subsp_ploidy_t[,order(colnames(sigGC.2021_subsp_ploidy_t))]
rownames(sigGC.2021_subsp_ploidy_t) <- c("sig_rank")

#write.csv(ANCOM_subsp_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_subsp_ploidy.csv")

GC.2021.OCG_sigsbst_subspploidy <-  t(subset(t(OCG_GC_2021_sbst), colnames(OCG_GC_2021_sbst) %in% row.names(sigGC.2021_subsp_ploidy)))

GC.2021.OCG_sigsbst_subspploidy <- GC.2021.OCG_sigsbst_subspploidy[,order(colnames(GC.2021.OCG_sigsbst_subspploidy))]
colnames(sigGC.2021_subsp_ploidy_t) == colnames(GC.2021.OCG_sigsbst_subspploidy) #sanity check:TRUE

GC.2021.OCG_sigsbst_subspploidy <- rbind(GC.2021.OCG_sigsbst_subspploidy, sigGC.2021_subsp_ploidy_t)
GC.2021.OCG_sigsbst_subspploidy_t <- as.data.frame(t(GC.2021.OCG_sigsbst_subspploidy))
GC.2021.OCG_sigsbst_subspploidy_t$sig_rank <- as.numeric(GC.2021.OCG_sigsbst_subspploidy_t$sig_rank) 
GC.2021.OCG_sigsbst_subspploidy_t <- GC.2021.OCG_sigsbst_subspploidy_t[order(GC.2021.OCG_sigsbst_subspploidy_t$sig_rank),] 
GC.2021.OCG_sigsbst_subspploidy_t <- subset(GC.2021.OCG_sigsbst_subspploidy_t, select=-c(sig_rank))
GC.2021.OCG_sigsbst_subspploidy_t <- as.data.frame(t(GC.2021.OCG_sigsbst_subspploidy_t))

#Build objects for plotting
sbsp.plplotGC2021 <- data.frame(GC.2021.OCG_sigsbst_subspploidy_t[,1:10], "subsp_ploidy" = md.OCG.GC.2021_sbst$Subsp_ploidy, check.names = FALSE)
sbsp.plplotGC2021[,1:10] <- lapply(sbsp.plplotGC2021[,1:10], function(x) as.numeric(as.character(x)))
sbsp.plplotGC2021[,1:10] <- log(sbsp.plplotGC2021[,1:10]+1)
sbsp.plplotGC2021_sub <- data.frame(sample=rownames(sbsp.plplotGC2021),sbsp.plplotGC2021, check.names = F)
sbsp.plplotGC2021_sublong <- melt(sbsp.plplotGC2021_sub)

ANCOM_subsp_ploidy.GC21$W.taxa

sbsp.plplotGC2021$subsp_ploidy <-  factor(sbsp.plplotGC2021$subsp_ploidy, levels = c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))

#### FIGURE FOR SUBSPECIES PLOIDY 2021 GC ####
spl2 <- ggplot(sbsp.plplotGC2021_sublong, aes(y = value, x = subsp_ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("2021 GC ANCOM for Subspecies and Ploidy")+
  labs(color = "Compounds") +
  theme_classic()

# LCMS ANCOM ####
summary(rowSums(OCG_LCMS_3uL_subset)) #31016795
summary(colSums(OCG_LCMS_3uL_subset)) #384380

OCG_LCMS_3uL_subset.a <- OCG_LCMS_3uL_subset[,colSums(OCG_LCMS_3uL_subset) > 10]
OCG_LCMS_3uL_subset.a <- OCG_LCMS_3uL_subset.a[rowSums(OCG_LCMS_3uL_subset.a) > 0,] 
md.OCG.LCMS.3.a <- subset(md.OCG.LCMS.3, row.names(md.OCG.LCMS.3) %in% row.names(OCG_LCMS_3uL_subset.a)) 

OCG_LCMS_3uL_subset.a.t <- t(OCG_LCMS_3uL_subset.a) 
OCG_LCMS_3uL_subset.a_t <- OCG_LCMS_3uL_subset.a.t[, colnames(OCG_LCMS_3uL_subset.a.t) %in% row.names(md.OCG.LCMS.3.a), drop = FALSE]
OCG_LCMS_3uL_subset.a <- t(OCG_LCMS_3uL_subset.a_t)
OCG_LCMS_3uL_subset.a <- as.data.frame(OCG_LCMS_3uL_subset.a) 

OCG_LCMS_3uL_subset.a[OCG_LCMS_3uL_subset.a < 10] <- 0  
OCG_LCMS_3uL_subset.a <- OCG_LCMS_3uL_subset.a[rowSums(OCG_LCMS_3uL_subset.a) > 0,] 

summary(rowSums(OCG_LCMS_3uL_subset.a)) #31016795
summary(colSums(OCG_LCMS_3uL_subset.a)) #384380

OCG_LCMS_3uL_subset.a <- OCG_LCMS_3uL_subset.a[,colSums(OCG_LCMS_3uL_subset.a) > 10]

md.OCG.LCMS.3_sbst <- data.frame("Sample.ID" = row.names(md.OCG.LCMS.3.a), md.OCG.LCMS.3.a)
OCG_LCMS_sbst <- data.frame("Sample.ID" = row.names(OCG_LCMS_3uL_subset.a), OCG_LCMS_3uL_subset.a, check.names = F)
row.names(OCG_LCMS_sbst) == row.names(md.OCG.LCMS.3_sbst) #TRUE

## SUBSPECIES LCMS ####
ANCOM_subspecies_LCMS <- ANCOM.main(OCG_LCMS_sbst,md.OCG.LCMS.3_sbst,F,F,"Subspecies",NULL,NULL,F,NULL,2,.05,.9)

#Create objects of significant compounds
sigLCMS_subspecies <- subset(ANCOM_subspecies_LCMS$W.taxa, ANCOM_subspecies_LCMS$W.taxa$W_stat > 0)[,1]
sigLCMS_subspecies <- as.data.frame(sigLCMS_subspecies) 
row.names(sigLCMS_subspecies) <- sigLCMS_subspecies[302:1,1] 

sigLCMS_subspecies[,1] <- c(302:1) 
sigLCMS_subspecies_t <- t(sigLCMS_subspecies) 
sigLCMS_subspecies_t <- as.data.frame(sigLCMS_subspecies_t)
colnames(sigLCMS_subspecies_t) <- as.character(colnames(sigLCMS_subspecies_t))
print(colnames(sigLCMS_subspecies_t))

sigLCMS_subspecies_t <- sigLCMS_subspecies_t[,order(colnames(sigLCMS_subspecies_t))]
rownames(sigLCMS_subspecies_t) <- c("sig_rank")

LCMS.OCG_sigsbst_subspecies <-  t(subset(t(OCG_LCMS_sbst), colnames(OCG_LCMS_sbst) %in% row.names(sigLCMS_subspecies)))

LCMS.OCG_sigsbst_subspecies <- LCMS.OCG_sigsbst_subspecies[,order(colnames(LCMS.OCG_sigsbst_subspecies))]
colnames(sigLCMS_subspecies_t) == colnames(LCMS.OCG_sigsbst_subspecies) #sanity check:TRUE

LCMS.OCG_sigsbst_subspecies <- rbind(LCMS.OCG_sigsbst_subspecies, sigLCMS_subspecies_t)
LCMS.OCG_sigsbst_subspecies_t <- as.data.frame(t(LCMS.OCG_sigsbst_subspecies))
LCMS.OCG_sigsbst_subspecies_t$sig_rank <- as.numeric(LCMS.OCG_sigsbst_subspecies_t$sig_rank) 
LCMS.OCG_sigsbst_subspecies_t <- LCMS.OCG_sigsbst_subspecies_t[order(LCMS.OCG_sigsbst_subspecies_t$sig_rank),] 
LCMS.OCG_sigsbst_subspecies_t <- subset(LCMS.OCG_sigsbst_subspecies_t, select=-c(sig_rank))
LCMS.OCG_sigsbst_subspecies_t <- as.data.frame(t(LCMS.OCG_sigsbst_subspecies_t))

#Build objects for plotting
sbsplotLCMS <- data.frame(LCMS.OCG_sigsbst_subspecies_t[,1:10], "subspecies" = md.OCG.LCMS.3_sbst$Subspecies, check.names = FALSE)
sbsplotLCMS[,1:10] <- lapply(sbsplotLCMS[,1:10], function(x) as.numeric(as.character(x)))
sbsplotLCMS[,1:10] <- log(sbsplotLCMS[,1:10]+1)
sbsplotLCMS_sub <- data.frame(sample=rownames(sbsplotLCMS),sbsplotLCMS, check.names = F)
sbsplotLCMS_sublong <- melt(sbsplotLCMS_sub)

ANCOM_subspecies_LCMS$W.taxa

sbsplotLCMS_sublong$subspecies <-  factor(sbsplotLCMS_sublong$subspecies, levels = c("T", "V", "W"))

#### SUBSPECIES FIGURE TOP ANCOM COMPOUNDS ALL LCMS ####
ggplot(sbsplotLCMS_sublong, aes(y = value, x = subspecies, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("LCMS ANCOM for Subspecies")+ 
  theme_classic()

## PLOIDY LCMS ####
#Run ANCOM, specify variable
ANCOM_ploidy_LCMS <- ANCOM.main(OCG_LCMS_sbst,md.OCG.LCMS.3_sbst,F,F,"Ploidy",NULL,NULL,F,NULL,2,.05,.9)

sigLCMS_ploidy <- subset(ANCOM_ploidy_LCMS$W.taxa, ANCOM_ploidy_LCMS$W.taxa$W_stat > 0)[,1]
sigLCMS_ploidy <- as.data.frame(sigLCMS_ploidy)
row.names(sigLCMS_ploidy) <- sigLCMS_ploidy[302:1,1] 

sigLCMS_ploidy[,1] <- c(302:1) 
sigLCMS_ploidy_t <- t(sigLCMS_ploidy) 
sigLCMS_ploidy_t <- as.data.frame(sigLCMS_ploidy_t)
colnames(sigLCMS_ploidy_t) <- as.character(colnames(sigLCMS_ploidy_t))
print(colnames(sigLCMS_ploidy_t))

sigLCMS_ploidy_t <- sigLCMS_ploidy_t[,order(colnames(sigLCMS_ploidy_t))]
rownames(sigLCMS_ploidy_t) <- c("sig_rank")

LCMS.OCG_sigsbst_ploidy <-  t(subset(t(OCG_LCMS_sbst), colnames(OCG_LCMS_sbst) %in% row.names(sigLCMS_ploidy)))

LCMS.OCG_sigsbst_ploidy <- LCMS.OCG_sigsbst_ploidy[,order(colnames(LCMS.OCG_sigsbst_ploidy))]
colnames(sigLCMS_ploidy_t) == colnames(LCMS.OCG_sigsbst_ploidy) #sanity check:TRUE

LCMS.OCG_sigsbst_ploidy <- rbind(LCMS.OCG_sigsbst_ploidy, sigLCMS_ploidy_t)
LCMS.OCG_sigsbst_ploidy_t <- as.data.frame(t(LCMS.OCG_sigsbst_ploidy))
LCMS.OCG_sigsbst_ploidy_t$sig_rank <- as.numeric(LCMS.OCG_sigsbst_ploidy_t$sig_rank) 
LCMS.OCG_sigsbst_ploidy_t <- LCMS.OCG_sigsbst_ploidy_t[order(LCMS.OCG_sigsbst_ploidy_t$sig_rank),] 
LCMS.OCG_sigsbst_ploidy_t <- subset(LCMS.OCG_sigsbst_ploidy_t, select=-c(sig_rank))
LCMS.OCG_sigsbst_ploidy_t <- as.data.frame(t(LCMS.OCG_sigsbst_ploidy_t))

#Build objects for plotting
plplotLCMS.ploidy <- data.frame(LCMS.OCG_sigsbst_ploidy_t[,1:10], "ploidy" = md.OCG.LCMS.3_sbst$Ploidy, check.names = FALSE)
plplotLCMS.ploidy[,1:10] <- lapply(plplotLCMS.ploidy[,1:10], function(x) as.numeric(as.character(x)))
plplotLCMS.ploidy[,1:10] <- log(plplotLCMS.ploidy[,1:10]+1)
plplotLCMS.ploidy_sub <- data.frame(sample=rownames(plplotLCMS.ploidy),plplotLCMS.ploidy, check.names = F)
plplotLCMS.ploidy_sublong <- melt(plplotLCMS.ploidy_sub)

ANCOM_ploidy_LCMS$W.taxa

plplotLCMS.ploidy_sublong$ploidy<-  factor(plplotLCMS.ploidy_sublong$ploidy, levels = c("2n", "4n"))

#### PLOIDY FIGURE TOP ANCOM COMPOUNDS LCMS ####
ggplot(plplotLCMS.ploidy_sublong, aes(y = value, x = ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("LCMS ANCOM for Ploidy")+
  labs(color = "Compounds") +
  theme_classic()


## SUBSPECIES PLOIDY LCMS ####
#Run ANCOM, specify variable
ANCOM_subsp_ploidy.LCMS <- ANCOM.main(OCG_LCMS_sbst,md.OCG.LCMS.3_sbst,F,F,"Subsp_ploidy",NULL,NULL,F,NULL,2,.05,.9)

sigLCMS_subsp_ploidy <- subset(ANCOM_subsp_ploidy.LCMS$W.taxa, ANCOM_subsp_ploidy.LCMS$W.taxa$W_stat > 0)[,1]
sigLCMS_subsp_ploidy <- as.data.frame(sigLCMS_subsp_ploidy)
row.names(sigLCMS_subsp_ploidy) <- sigLCMS_subsp_ploidy[302:1,1]
sigLCMS_subsp_ploidy[,1] <- c(302:1)

sigLCMS_subsp_ploidy_t <- t(sigLCMS_subsp_ploidy)
sigLCMS_subsp_ploidy_t <- as.data.frame(sigLCMS_subsp_ploidy_t)
colnames(sigLCMS_subsp_ploidy_t) <- as.character(colnames(sigLCMS_subsp_ploidy_t))
print(colnames(sigLCMS_subsp_ploidy_t))
sigLCMS_subsp_ploidy_t <- sigLCMS_subsp_ploidy_t[,order(colnames(sigLCMS_subsp_ploidy_t))]
rownames(sigLCMS_subsp_ploidy_t) <- c("sig_rank")

#write.csv(ANCOM_subsp_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_subsp_ploidy.csv")

LCMS_sigsbst_subspploidy <-  t(subset(t(OCG_LCMS_sbst), colnames(OCG_LCMS_sbst) %in% row.names(sigLCMS_subsp_ploidy)))

LCMS_sigsbst_subspploidy <- LCMS_sigsbst_subspploidy[,order(colnames(LCMS_sigsbst_subspploidy))]
colnames(sigLCMS_subsp_ploidy_t) == colnames(LCMS_sigsbst_subspploidy) #sanity check:TRUE

LCMS_sigsbst_subspploidy <- rbind(LCMS_sigsbst_subspploidy, sigLCMS_subsp_ploidy_t)
LCMS_sigsbst_subspploidy_t <- as.data.frame(t(LCMS_sigsbst_subspploidy))
LCMS_sigsbst_subspploidy_t$sig_rank <- as.numeric(LCMS_sigsbst_subspploidy_t$sig_rank) 
LCMS_sigsbst_subspploidy_t <- LCMS_sigsbst_subspploidy_t[order(LCMS_sigsbst_subspploidy_t$sig_rank),] 
LCMS_sigsbst_subspploidy_t <- subset(LCMS_sigsbst_subspploidy_t, select=-c(sig_rank))
LCMS_sigsbst_subspploidy_t <- as.data.frame(t(LCMS_sigsbst_subspploidy_t))

#Build objects for plotting
plotLCMS.subspploidy <- data.frame(LCMS_sigsbst_subspploidy_t[,1:10], "subsp_ploidy" = md.OCG.LCMS.3_sbst$Subsp_ploidy, check.names = FALSE)
plotLCMS.subspploidy[,1:10] <- lapply(plotLCMS.subspploidy[,1:10], function(x) as.numeric(as.character(x)))
plotLCMS.subspploidy[,1:10] <- log(plotLCMS.subspploidy[,1:10]+1)
plotLCMS.subspploidy_sub <- data.frame(sample=rownames(plotLCMS.subspploidy),plotLCMS.subspploidy, check.names = F)
plotLCMS.subspploidy_sublong <- melt(plotLCMS.subspploidy_sub)

ANCOM_subsp_ploidy$W.taxa

plotLCMS.subspploidy$subsp_ploidy <-  factor(plotLCMS.subspploidy$subsp_ploidy, levels = c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))

### FIGURE FOR SUBSPECIES PLOIDY ####
ggplot(plotLCMS.subspploidy_sublong, aes(y = value, x = subsp_ploidy, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("LCMS ANCOM for Subspecies and Ploidy")+
  labs(color = "Compounds") +
  theme_classic()

## YEAR LCMS ####
#Run ANCOM, specify variable
ANCOM_yr.LCMS <- ANCOM.main(OCG_LCMS_sbst,md.OCG.LCMS.3_sbst,F,F,"Year",NULL,NULL,F,NULL,2,.05,.9)

sigLCMS_yr <- subset(ANCOM_yr.LCMS$W.taxa, ANCOM_yr.LCMS$W.taxa$W_stat > 0)[,1]
sigLCMS_yr <- as.data.frame(sigLCMS_yr)
row.names(sigLCMS_yr) <- sigLCMS_yr[302:1,1]
sigLCMS_yr[,1] <- c(302:1)

sigLCMS_yr_t <- t(sigLCMS_yr)
sigLCMS_yr_t <- as.data.frame(sigLCMS_yr_t)
colnames(sigLCMS_yr_t) <- as.character(colnames(sigLCMS_yr_t))
print(colnames(sigLCMS_yr_t))
sigLCMS_yr_t <- sigLCMS_yr_t[,order(colnames(sigLCMS_yr_t))]
rownames(sigLCMS_yr_t) <- c("sig_rank")

#write.csv(ANCOM_subsp_ploidy$W.taxa, file = "data_csv/ANCOM/ANCOM_subsp_ploidy.csv")

LCMS_sig_yr <-  t(subset(t(OCG_LCMS_sbst), colnames(OCG_LCMS_sbst) %in% row.names(sigLCMS_yr)))

LCMS_sig_yr <- LCMS_sig_yr[,order(colnames(LCMS_sig_yr))]
colnames(sigLCMS_yr_t) == colnames(LCMS_sig_yr) #sanity check:TRUE

LCMS_sig_yr <- rbind(LCMS_sig_yr, sigLCMS_yr_t)
LCMS_sig_yr_t <- as.data.frame(t(LCMS_sig_yr))
LCMS_sig_yr_t$sig_rank <- as.numeric(LCMS_sig_yr_t$sig_rank) 
LCMS_sig_yr_t <- LCMS_sig_yr_t[order(LCMS_sig_yr_t$sig_rank),] 
LCMS_sig_yr_t <- subset(LCMS_sig_yr_t, select=-c(sig_rank))
LCMS_sig_yr_t <- as.data.frame(t(LCMS_sig_yr_t))

#Build objects for plotting
plotLCMS.yr <- data.frame(LCMS_sig_yr_t[,1:10], "year" = md.OCG.LCMS.3_sbst$Year, check.names = FALSE)
plotLCMS.yr[,1:10] <- lapply(plotLCMS.yr[,1:10], function(x) as.numeric(as.character(x)))
plotLCMS.yr[,1:10] <- log(plotLCMS.yr[,1:10]+1)
plotLCMS.yr_sub <- data.frame(sample=rownames(plotLCMS.yr),plotLCMS.yr, check.names = F)
plotLCMS.yr_sublong <- melt(plotLCMS.yr_sub)

plotLCMS.yr$year <-  factor(plotLCMS.yr$year, levels = c("2012","2021"))

### FIGURE FOR SUBSPECIES PLOIDY ####
ggplot(plotLCMS.yr_sublong, aes(y = value, x = year, color=variable))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_dodge(width=0.75), aes(group=variable), alpha =.4) +
  scale_color_brewer(palette = "Spectral")+
  ylab("Log  rel. abundance") + xlab("Subspecies") + 
  ggtitle("LCMS ANCOM for year")+
  labs(color = "Compounds") +
  theme_classic()

# ANCOM BC #### 
## Make phyloseq objects ####
row.names(OCG_GC_2012_subset) == row.names(md.OCG.GC.2012) #TRUE
compound_table_GC12 <- otu_table(as.matrix(OCG_GC_2012_subset), taxa_are_rows = FALSE) #create the feature table
md_GC12_phy <- sample_data(md.OCG.2012) #make metadata into sample data to create phyloseq
physeqGC12 <- phyloseq(compound_table, md_phy) #create phyloseq object for ANCOM BC

result_GC12 <- ancombc2( 
  data = physeqGC12, assay_name = "counts", tax_level = NULL,
  fix_formula = "Subspecies + Ploidy",
  rand_formula = "(1|Plant)",
  p_adj_method = "fdr",
  group = "Subspecies", #change variable of interest
  alpha = 0.05
)

#iter_control = list(tol = 1e-5, max_iter = 200, 
                    #verbose = FALSE)

ancomGC12res_df <- result_GC12$res

# LCMS SUBSPECIES ANCOM ####
# Create phyloseq object
row.names(OCG_LCMS_3uL_subset) == row.names(md.OCG.LCMS.3) #TRUE
compound_table_LCMS <- otu_table(as.matrix(OCG_LCMS_3uL_subset), taxa_are_rows = FALSE) #create the feature table
md_LCMS_phy <- sample_data(md.OCG.LCMS.3) #make metadata into sample data to create phyloseq
physeqLCMS <- phyloseq(compound_table_LCMS, md_LCMS_phy) #create phyloseq object for ANCOM BC

# Create TreeSummarizedExperimentFromPhyloseq with mia 
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(physeqLCMS)

result_LCMS_sub <- ancombc2( 
  data = tse, assay_name = "counts", tax_level = NULL,
  fix_formula = "Subspecies",
  rand_formula = "(1|Plant)",
  p_adj_method = "fdr",
  group = "Subspecies", #change variable of interest
  alpha = 0.05,
  iter_control = list(tol = 1e-5, max_iter = 200, 
                                   verbose = FALSE)
)

ancomLCMSres_df <- result_LCMS_sub$res








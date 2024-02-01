# Orchard Common Garden GC 2012 and 2021 Chemistry analysis
# Set working directory and load necessary packages for reading in and tidying the data.####
setwd("~/Documents/Orchard_Common _Garden/commongarden")
if (!require("readr")) {install.packages("readr"); require("readr")}
if (!require("dplyr")) {install.packages("dplyr"); require("dplyr")}
if (!require("tidyr")) {install.packages("tidyr"); require("tidyr")}
if (!require("factoextra")) {install.packages("factoextra"); require("factoextra")}
if (!require("vegan")) {install.packages("vegan"); require("vegan")}
if (!require("ggcorrplot")) {install.packages("ggcorrplot"); require("ggcorrplot")}
if (!require("ggpubr")) {install.packages("ggpubr"); require("ggpubr")}
if (!require("devtools")) {install.packages("devtools"); require("devtools")}
if (!require("pairwiseAdonis")) {devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"); require("pairwiseAdonis")}
if (!require("ggplot2")) {install.packages("ggplot2"); require("ggplot2")}

# Read data in#### 
##Metadata read in#####
mdITS_OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T)

## 2012 GC data read in and cleaning####
#read in 2012 GC chemistry data
OCG_AUC_2012 <- read.csv("data_csv/OCG_AUC_2012.csv", head=T, row.names = 1,check.names = F,stringsAsFactors = T) #This data has AUC. 167 pf 53 variables

## Subset to have the plants from metadata in chemistry data
#subset the md to only have observations from 2012 to avoid duplicates
mdITS_OCG_2012 <- subset(mdITS_OCG, mdITS_OCG$Year=="2012") #96 observations and 22 variables

names(OCG_AUC_2012)[names(OCG_AUC_2012) == 'garden_Plant_ID'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

OCG_AUC_2012$`Garden Plant ID` <- as.integer(OCG_AUC_2012$`Garden Plant ID`) #as integer so I can combine them

OCG_AUC_2012_subset <- left_join(mdITS_OCG_2012, OCG_AUC_2012, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 96 of 73 variables

OCG_AUC_2012_subset <- OCG_AUC_2012_subset[,-c(1,3:24)] #removing everything except area under the curve for each compound and garden plant ID number

##Save this csv so I can re-read in the data with row 1 being plant ID
#write.csv(OCG_AUC_2012_subset, file = "data_csv/OCG_AUC_subset_2012.csv",row.names = FALSE)
#write.csv(mdITS_OCG_2012, file = "data_csv/md_OCG_2012.csv",row.names = FALSE)

#read in metadata for just 2012 with row 1 being plant ID
mdITS_OCG_2012 <- read.csv("data_csv/md_OCG_2012.csv",head=T, check.names = F,stringsAsFactors = T) #96 of 22 variables
#read in subsetted 2012 GC data with row 1 being plant ID
OCG_AUC_2012_subset <- read.csv("data_csv/OCG_AUC_subset_2012.csv", row.names = 1) #96 of 49 variables

## 2021 GC data read in and cleaning####
#read in 2021 GC chemistry data#
OCG_AUC_2021 <- read.csv("data_csv/OCG_2021_GC_data.csv", head=T, check.names = F,stringsAsFactors = T) #143 of 223 variables

#subset the md to only have observations from 2021 to avoid duplicates
mdITS_OCG_2021 <- subset(mdITS_OCG, mdITS_OCG$Year=="2021") #43 observations and 22 variables

sum(duplicated(OCG_AUC_2021)) #no duplicates in the 2021 dataset

head(OCG_AUC_2021)

OCG_AUC_2021_subset <- OCG_AUC_2021[,c(1,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,147,150,153,156,159,162,165,168,171,174,177,180,183,186,189,192,195,198,201,204,207,210,213,216,219,222)] #removing everything except area under the curve for each compound and garden plant ID number

names(OCG_AUC_2021_subset)[names(OCG_AUC_2021_subset) == 'Plant_ID'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

OCG_AUC_2021_subset$`Garden Plant ID` <- as.integer(OCG_AUC_2021_subset$`Garden Plant ID`) #as integer so I can combine them

OCG_AUC_2021_subset <- left_join(mdITS_OCG_2021, OCG_AUC_2021_subset, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 43 of 96 variables

#Going to remove everything except chem data and plant ID
OCG_AUC_2021_subset <- OCG_AUC_2021_subset[,-c(1,3:22)] #43 of 75 variables

##Save this csv so I can re-read in the data with row 1 being plant ID
#write.csv(OCG_AUC_2021_subset, file = "data_csv/OCG_AUC_2021.csv",row.names = FALSE)
#write.csv(mdITS_OCG_2021, file = "data_csv/md_OCG_2021.csv",row.names = FALSE)

#read in metadata
mdITS_OCG_2021 <- read.csv("data_csv/md_OCG_2021.csv",head=T, check.names = F,stringsAsFactors = T) #43 of 22 variables

#read in 2021 chemistry data
OCG_AUC_2021_subset <- read.csv("data_csv/OCG_AUC_2021.csv", row.names = 1) #43 of 74 variables

#Cleaning####
## 2012 data cleaning ####
colSums(is.na(OCG_AUC_2012_subset)) #checking for null values since pca wont run with NAs. All zeroes for each column which indicates there are no NA values

## 2021 data cleaning ####
#So now the data needs to be cleaned to only contain compounds that occur in more than 20% of the samples (plants). AKA columns need to be removed that don't have at least 20% of the rows containing a number=NA

#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_AUC_2021_subset))

#Now define the threshold of 20% - update there are no compounds that occur in 20% of samples. 90% had 2 compounds left and 95% had 19 compounds
threshold <- 0.95

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_AUC_2021_subset_filtered <- OCG_AUC_2021_subset[, columns_to_keep] #43 of 19 variables

# Replace NA values with zeroes
OCG_AUC_2021_subset_filtered[is.na(OCG_AUC_2021_subset_filtered)] <- 0

colSums(is.na(OCG_AUC_2021_subset_filtered)) #all zeroes for each column which indicates there are no NA values

#Scaling####
## 2012 scaling ####
#Compound
data_normalized_2012 <- scale(OCG_AUC_2012_subset) #so this will center the data.... able to code this instead of using excel to do this
head(data_normalized_2012) #looks good

#Plant ID
#transpose the data first to put plant ID as columns and compound as row
OCG_AUC_2012_subset.t <- t(OCG_AUC_2012_subset) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_AUC_2012_subset.t)) #checking for null values

data_normalized_2012.ID <- scale(OCG_AUC_2012_subset.t) 
head(data_normalized_2012.ID) #looks good

## 2021 scaling ####
#Compound correlation plot
data_normalized_2021 <- scale(OCG_AUC_2021_subset_filtered) #so this will center the data.... able to code this instead of using excel to do this
head(data_normalized_2021) #mostly negative which is due to the lack of values

#Plant correlation plot
#transpose the data first to put plant ID as columns and compound as row 
OCG_AUC_2021_subset.t <- t(OCG_AUC_2021_subset_filtered) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_AUC_2021_subset.t)) #checking for null values

data_normalized_2021.ID <- scale(OCG_AUC_2021_subset.t) 
head(data_normalized_2021.ID) 

# Data visualization####
## Correlation plots####
### 2012 correlation plots ####
#corr plot for compound
corr_matrix_2012 <- cor(data_normalized_2012)
ggcorrplot(corr_matrix_2012)

#corr plot for plant ID
corr_matrix_2012.ID <- cor(data_normalized_2012.ID)
ggcorrplot(corr_matrix_2012.ID)

### 2021 correlation plots ####
#corr plot for compound
corr_matrix_2021 <- cor(data_normalized_2021)
ggcorrplot(corr_matrix_2021)

#corr plot for plant ID
corr_matrix_2021.ID <- cor(data_normalized_2021.ID)
ggcorrplot(corr_matrix_2021.ID) #nothing shows up here... may need to subset to only the 4 plants that would be included for 20% cut off

# PCA ####
## 2012 PCA by compound####
data.pca_2012 <- princomp(corr_matrix_2012)
summary(data.pca_2012)

data.pca_2012$loadings[, 1:2]

fviz_eig(data.pca_2012, addlabels = TRUE) #scree plot is used to visualize the importance of each principal component and can be used to determine the number of principal components to retain.

#With the biplot, it is possible to visualize the similarities and dissimilarities between the samples, and further shows the impact of each attribute on each of the principal components.

# Graph of the compounds
fviz_pca_var(data.pca_2012, col.var = "black")

#Contribution of each compound
fviz_cos2(data.pca_2012, choice = "var", axes = 1:2)

#Biplot combined with cos2
fviz_pca_var(data.pca_2012, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

## 2021 PCA by compound####
data.pca_2021 <- princomp(corr_matrix_2021)
summary(data.pca_2021)

data.pca_2021$loadings[, 1:2]

fviz_eig(data.pca_2021, addlabels = TRUE) #scree plot is used to visualize the importance of each principal component and can be used to determine the number of principal components to retain.

#With the biplot, it is possible to visualize the similarities and dissimilarities between the samples, and further shows the impact of each attribute on each of the principal components.

# Graph of the compounds
fviz_pca_var(data.pca_2021, col.var = "black")

#Contribution of each compound
fviz_cos2(data.pca_2021, choice = "var", axes = 1:2)

#Biplot combined with cos2
fviz_pca_var(data.pca_2021, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

## 2012 PCA by Plant ID ####
data.pca_2012_ID <- princomp(corr_matrix_2012.ID)
summary(data.pca_2012_ID)

data.pca_2012_ID$loadings[, 1:2] 

fviz_eig(data.pca_2012_ID, addlabels = TRUE) #scree plot 

# Graph of the variables
fviz_pca_var(data.pca_2012_ID, col.var = "black") #pca

#Contribution of each compound
fviz_cos2(data.pca_2012_ID, choice = "var", axes = 1:2)

#Biplot combined with cos2
fviz_pca_var(data.pca_2012_ID, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

## 2021 PCA by Plant ID isn't working ####
data.pca_2021_ID <- princomp(corr_matrix_2021.ID)
summary(data.pca_2012_ID)

data.pca_2012_ID$loadings[, 1:2] 

fviz_eig(data.pca_2012_ID, addlabels = TRUE) #scree plot 

# Graph of the variables
fviz_pca_var(data.pca_2012_ID, col.var = "black") #pca

## NMDS plots####
## 2012 NMDS by Compound ####
mdITS_OCG_2012$Subspecies <- as.factor(mdITS_OCG_2012$Subspecies)
mdITS_OCG_2012$Subsp_ploidy <- as.factor(mdITS_OCG_2012$Subsp_ploidy)
mdITS_OCG_2012$Year <- as.factor(mdITS_OCG_2012$Year)
mdITS_OCG_2012$Ploidy <- as.factor(mdITS_OCG_2012$Ploidy)

#turn abundance data frame into a matrix
m_OCG_AUC_2012_subset = as.matrix(OCG_AUC_2012_subset)

set.seed(65)
#OCG_AUC_2012_subset.nmds <- metaMDS(t(m_OCG_AUC_2012_subset), trymax=500) #solution reached.
#save(OCG_AUC_2012_subset.nmds, file = "nmds/OCG_AUC_2012_subset.nmds.rda")
load("nmds/OCG_AUC_2012_subset.nmds.rda")

ordiplot(OCG_AUC_2012_subset.nmds, type = "t",display = "sites",cex = .6) 

## 2012 NMDS by plant ID ####
OCG_AUC_2012_subset.t <- t(OCG_AUC_2012_subset)
m_OCG_AUC_2012_subset.t = as.matrix(OCG_AUC_2012_subset.t)

set.seed(6)
#Plant_ID_OCG_AUC_2012.nmds <- metaMDS(t(m_OCG_AUC_2012_subset.t), trymax=500) #solution reached.
#save(Plant_ID_OCG_AUC_2012.nmds, file = "nmds/Plant_ID_OCG_AUC.nmds.rda")
load("nmds/Plant_ID_OCG_AUC.nmds.rda")

ordiplot(Plant_ID_OCG_AUC.nmds, type = "t",display = "sites",cex = .7)

### Ploidy NMDS####
mdITS_OCG_2012$Ploidy<- droplevels(mdITS_OCG_2012$Ploidy)
levels(mdITS_OCG_2012$Ploidy) #Always smart to check levels before plotting. If there are too many or too few levels, the nmds plots wont line up with the ordispider function.

plot(Plant_ID_OCG_AUC.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Ploidy", 
     col= c("red","blue")[mdITS_OCG_2012$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(Plant_ID_OCG_AUC.nmds,groups = mdITS_OCG_2012$Ploidy, show.groups = "2n", col = "red")
ordispider(Plant_ID_OCG_AUC.nmds,groups = mdITS_OCG_2012$Ploidy, show.groups = "4n", col = "blue")

#### PERMANOVA for 2012 ploidy ####
OCG_AUC_2012_ploidy <- adonis2(OCG_AUC_2012_subset ~ mdITS_OCG_2012$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_ploidy #ploidy significant 0.034

### Subspecies NMDS ####
mdITS_OCG_2012$Subspecies<- droplevels(mdITS_OCG_2012$Subspecies)
levels(mdITS_OCG_2012$Subspecies)

plot(Plant_ID_OCG_AUC.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="Plant chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS_OCG_2012$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(Plant_ID_OCG_AUC.nmds,groups = mdITS_OCG_2012$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(Plant_ID_OCG_AUC.nmds,groups = mdITS_OCG_2012$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(Plant_ID_OCG_AUC.nmds,groups = mdITS_OCG_2012$Subspecies, show.groups = "W", col = "goldenrod")

#### PERMANOVA & pairwaise adonis for subspecies ####
OCG_AUC_2012_subsp <- adonis2(OCG_AUC_2012_subset ~ mdITS_OCG_2012$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_subsp #subspecies significant 0.001

#pairwiseadonis
OCG_AUC_2012_subsp.pw <- pairwise.adonis(OCG_AUC_2012_subset, mdITS_OCG_2012$Subspecies)
OCG_AUC_2012_subsp.pw#T vs V= 0.003, T vs W= 0.120, and W vs V= 0.003.

### Subspecies ploidy NMDS ####
mdITS_OCG_2012$Subsp_ploidy<- droplevels(mdITS_OCG_2012$Subsp_ploidy)
levels(mdITS_OCG_2012$Subsp_ploidy)

plot(Plant_ID_OCG_AUC.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS_OCG_2012$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(Plant_ID_OCG_AUC.nmds,groups = mdITS_OCG_2012$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(Plant_ID_OCG_AUC.nmds,groups = mdITS_OCG_2012$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(Plant_ID_OCG_AUC.nmds,groups = mdITS_OCG_2012$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(Plant_ID_OCG_AUC.nmds,groups = mdITS_OCG_2012$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(Plant_ID_OCG_AUC.nmds,groups = mdITS_OCG_2012$Subsp_ploidy, show.groups = "W_4n", col = "purple")

#### PERMANOVAS and pairwise adonis for subspecies ploidy ####
OCG_AUC_2012_subspploidy <- adonis2(OCG_AUC_2012_subset ~ mdITS_OCG_2012$Subsp_ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_subspploidy #subspecies ploidy significant 0.001

#pairwiseadonis
OCG_AUC_2012_subsp_ploidy.pw <- pairwise.adonis(OCG_AUC_2012_subset, mdITS_OCG_2012$Subsp_ploidy)
OCG_AUC_2012_subsp_ploidy.pw #T_4n vs T_2n=0.02, T_4n vs V_2n = 0.05, T_2n vs V_2n = 0.01, T_2n vs V_4n = 0.01, T_2n vs W_4n = 0.01, V_4n vs W_4n = 0.06, and V_2n vs W_4n= 0.01.

## 2021 NMDS by Compound####
mdITS_OCG_2021$Subspecies <- as.factor(mdITS_OCG_2021$Subspecies)
mdITS_OCG_2021$Subsp_ploidy <- as.factor(mdITS_OCG_2021$Subsp_ploidy)
mdITS_OCG_2021$Year <- as.factor(mdITS_OCG_2021$Year)
mdITS_OCG_2021$Ploidy <- as.factor(mdITS_OCG_2021$Ploidy)

#turn abundance data frame into a matrix
m_OCG_AUC_2021_subset_filtered = as.matrix(OCG_AUC_2021_subset_filtered)

#set.seed(82)
#OCG_AUC_2021_subset_filtered.nmds <- metaMDS(t(m_OCG_AUC_2021_subset_filtered), trymax=500) #solution reached.
#save(OCG_AUC_2021_subset_filtered.nmds, file = "nmds/OCG_AUC_2021_subset_filtered.nmds.rda")
load("nmds/OCG_AUC_2021_subset_filtered.nmds.rda")

#Ordiplot by compound
ordiplot(OCG_AUC_2021_subset_filtered.nmds, type = "t",display = "sites",cex = .6) 

ordiplot(OCG_AUC_2021_subset_filtered.nmds, type = "t",display = "species",cex = .6) 

## 2021 NMDS by Plant ID isn't working####
AUC.OCG <- adonis2(OCG_AUC_2021_subset_filtered ~ mdITS_OCG_2021$Subspecies,by="margin") 
#turn abundance data frame into a matrix
OCG_AUC_2021_subset_filtered.t <- t(OCG_AUC_2021_subset_filtered)

m_OCG_AUC_2021_subset_filtered.t = as.matrix(OCG_AUC_2021_subset_filtered.t)

colSums(is.na(m_OCG_AUC_2021_subset_filtered.t))
set.seed(4)
Plant_ID_OCG_AUC_2021.nmds <- metaMDS(t(m_OCG_AUC_2021_subset_filtered.t), trymax=500) #solution reached.
# save(Plant_ID_OCG_AUC.nmds, file = "nmds/Plant_ID_OCG_AUC.nmds.rda")
# load("nmds/Plant_ID_OCG_AUC.nmds.rda")
# 
# #Ordiplot for plant ID
# ordiplot(Plant_ID_OCG_AUC.nmds, type = "t",display = "sites",cex = .7)


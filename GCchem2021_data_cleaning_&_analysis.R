# Orchard Common Garden Chemistry 2021 GC data cleaning####
## Set working directory and load necessary packages for reading in and tidying the data.####
setwd("~/Documents/Orchard_Common _Garden/commongarden")
if (!require("readr")) {install.packages("readr"); require("readr")}
if (!require("dplyr")) {install.packages("dplyr"); require("dplyr")}
if (!require("tidyr")) {install.packages("tidyr"); require("tidyr")}
if (!require("readxl")) {install.packages("readxl"); require("readxl")} #shouldn't need this hopefully

## Read data in#### 
#read in metadata#
mdITS_OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T)

#read in chemistry data#
OCG_AUC_2021 <- read.csv("data_csv/OCG_2021_GC_data.csv", head=T, check.names = F,stringsAsFactors = T) #This data has AUC, retention time, and peak area percent

#I need to subset the md to only have observations from 2021 to avoid duplicates
mdITS_OCG_2021 <- subset(mdITS_OCG, mdITS_OCG$Year=="2021") #43 observations and 22 variables

sum(duplicated(OCG_AUC_2021)) #no duplicates in the 2021 dataset

View(OCG_AUC_2021)

OCG_AUC_2021_subset <- OCG_AUC_2021[,c(1,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,147,150,153,156,159,162,165,168,171,174,177,180,183,186,189,192,195,198,201,204,207,210,213,216,219,222)] #removing everything except area under the curve for each compound and garden plant ID number

names(OCG_AUC_2021_subset)[names(OCG_AUC_2021_subset) == 'Plant_ID'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

OCG_AUC_2021_subset$`Garden Plant ID` <- as.integer(OCG_AUC_2021_subset$`Garden Plant ID`) #as integer so I can combine them

OCG_AUC_2021_subset <- left_join(mdITS_OCG_2021, OCG_AUC_2021_subset, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 43 of 96 variables

#Going to remove everything except chem data and plant ID

OCG_AUC_2021_subset <- OCG_AUC_2021_subset[,-c(1,3:22)] #43 of 75 variables

##Save this csv so I can re-read in the data with row 1 being plant ID
write.csv(OCG_AUC_2021_subset, file = "data_csv/OCG_AUC_2021.csv",row.names = FALSE)
write.csv(mdITS_OCG_2021, file = "data_csv/md_OCG_2021.csv",row.names = FALSE)

# OCG 2021 chemistry data analyses and visualization####
## Set working directory and load necessary packages for reading in and tidying the data.####
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

##Read in data####
#read in metadata
mdITS_OCG_2021 <- read.csv("data_csv/md_OCG_2021.csv",head=T, check.names = F,stringsAsFactors = T) #96 of 22 variables

#read in chemistry data
OCG_AUC_2021_subset <- read.csv("data_csv/OCG_AUC_2021.csv", row.names = 1) #43 of 74 variables

##Cleaning chem data####

#So now the data needs to be cleaned to only contain compounds that occur in more than 20% of the samples. AKA columns need to be removed that dont have at least 20% of the rows containing a number=NA

#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_AUC_2021_subset))

#Now define the threshold of 20% - update there are no compounds that occur in 20% of samples. 90% had 2 compounds left and 95% had 19 compounds
threshold <- 0.95

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_AUC_2021_subset_filtered <- OCG_AUC_2021_subset[, columns_to_keep]

View(OCG_AUC_2021_subset_filtered)

###Scaling####
#Since PCAs can only be run on numeric data, always check for NAs in the data

# Replace NA values with zeroes
OCG_AUC_2021_subset_filtered[is.na(OCG_AUC_2021_subset_filtered)] <- 0

#only works with numeric data
colSums(is.na(OCG_AUC_2021_subset_filtered)) #checking for null values since pca wont run with NAs. Should have all zeroes for each column which indictes there are no NA values

data_normalized <- scale(OCG_AUC_2021_subset_filtered) #so this will center the data.... able to code this instead of using excel to do this
head(data_normalized) #looks good

## Compound correlation plots####
corr_matrix <- cor(data_normalized)
ggcorrplot(corr_matrix)

## PCAs####
### PCA by compound####
data.pca <- princomp(corr_matrix)
summary(data.pca)

data.pca$loadings[, 1:2]

fviz_eig(data.pca, addlabels = TRUE) #scree plot is used to visualize the importance of each principal component and can be used to determine the number of principal components to retain.

#With the biplot, it is possible to visualize the similarities and dissimilarities between the samples, and further shows the impact of each attribute on each of the principal components.

# Graph of the compounds
fviz_pca_var(data.pca, col.var = "black")

#Contribution of each compound
fviz_cos2(data.pca, choice = "var", axes = 1:2)

#Biplot combined with cos2
fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
## Compound correlation plots####

#First the AUC data needs to be transposed in order to plot the PCA

#transpose the data first to put plant ID as columns and compound as row
OCG_cen_subset.t <- t(OCG_AUC_2021_subset_filtered) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_cen_subset.t)) #checking for null values since pca wont run with NAs

data_normalized.ID <- scale(OCG_cen_subset.t) 
head(data_normalized.ID) #looks good
corr_matrix.ID <- cor(data_normalized.ID)
ggcorrplot(corr_matrix.ID)

### PCA by Plant ID####

data.pca_ID <- princomp(corr_matrix.ID)
summary(data.pca_ID)

data.pca_ID$loadings[, 1:2] #need to talk about this

fviz_eig(data.pca_ID, addlabels = TRUE) #scree plot 

# Graph of the variables
fviz_pca_var(data.pca_ID, col.var = "black") #pca
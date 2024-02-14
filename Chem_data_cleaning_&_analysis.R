# Orchard Common Garden GC 2012 and 2021 Chemistry analysis
# Set working directory and load necessary packages for reading in and tidying the data.####
setwd("~/Documents/Orchard_Common _Garden/commongarden")
if (!require("readr")) {install.packages("readr"); require("readr")}
if (!require("dplyr")) {install.packages("dplyr"); require("dplyr")}
if (!require("tidyr")) {install.packages("tidyr"); require("tidyr")}
if (!require("tidyverse")) {install.packages("tidyverse"); require("tidyverse")}
if (!require("factoextra")) {install.packages("factoextra"); require("factoextra")}
if (!require("vegan")) {install.packages("vegan"); require("vegan")}
if (!require("ggcorrplot")) {install.packages("ggcorrplot"); require("ggcorrplot")}
if (!require("ggpubr")) {install.packages("ggpubr"); require("ggpubr")}
if (!require("devtools")) {install.packages("devtools"); require("devtools")}
if (!require("pairwiseAdonis")) {devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"); require("pairwiseAdonis")}
if (!require("ggplot2")) {install.packages("ggplot2"); require("ggplot2")}
if (!require("xcms"))  {install.packages("xcms"); require("xcms")}
# install.packages("devtools") 
# devtools::install_github("mottensmann/GCalignR", build_vignettes = TRUE) 
library("GCalignR") 
library("xcms")

# Read data in#### 
##Metadata read in#####
mdITS_OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #139 of 22 variables

#subset the md to only have observations from 2012 to avoid duplicates
mdITS_OCG_2012 <- subset(mdITS_OCG, mdITS_OCG$Year=="2012") #96 observations and 22 variables

#subset the md to only have observations from 2021 to avoid duplicates
mdITS_OCG_2021 <- subset(mdITS_OCG, mdITS_OCG$Year=="2021") #43 observations and 22 variables

## 2012 GC raw data read in####
#read in 2012 GC chemistry data
#OCG_AUC_2012 <- read.csv("data_csv/OCG_AUC_2012.csv", head=T, row.names = 1,check.names = F,stringsAsFactors = T) #This data has AUC. 167 pf 53 variables

## 2012 GC Cleaning ####
#names(OCG_AUC_2012)[names(OCG_AUC_2012) == 'garden_Plant_ID'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

#OCG_AUC_2012$`Garden Plant ID` <- as.integer(OCG_AUC_2012$`Garden Plant ID`) #as integer so I can combine them

#OCG_AUC_2012_subset <- left_join(mdITS_OCG_2012, OCG_AUC_2012, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 96 of 73 variables

#OCG_AUC_2012_subset <- OCG_AUC_2012_subset[,-c(1,3:24)] #removing everything except area under the curve for each compound and garden plant ID number

##Save this csv so I can re-read in the data with row 1 being plant ID
#write.csv(OCG_AUC_2012_subset, file = "data_csv/OCG_AUC_subset_2012.csv",row.names = FALSE)
#write.csv(mdITS_OCG_2012, file = "data_csv/md_OCG_2012.csv",row.names = FALSE)

#read in metadata for just 2012 with row 1 being plant ID
#mdITS_OCG_2012 <- read.csv("data_csv/md_OCG_2012.csv",head=T, check.names = F,stringsAsFactors = T) #96 of 22 variables

## 2012 cleaned GC data read in####
#read in subsetted 2012 GC data with row 1 being plant ID
OCG_AUC_2012 <- read.csv("data_csv/OCG_AUC_subset_2012.csv", row.names = 1) #96 of 49 variables

## 2021 GC raw data read in####
#read in 2021 GC chemistry data#
#OCG_AUC_2021 <- read.csv("data_csv/OCG_2021_GCData.csv", head=T, skip = 1) #143 of 223 variables

## 2021 GC Cleaning ####
#sum(duplicated(OCG_AUC_2021)) #no duplicates in the 2021 dataset

#head(OCG_AUC_2021)

#subset to only columns that contain "Peak Area" and "Plant ID". This removes "RT". 
#OCG_AUC_2021 <- OCG_AUC_2021 [, grepl("Peak.Area|Plant_ID", colnames(OCG_AUC_2021))] #143 obs of 149 variables

#subset to remove columns that contain Peak Area Percent and just keep Peak Area.
#OCG_AUC_2021 <- OCG_AUC_2021 [, !grepl("Peak.Area.Percent", colnames(OCG_AUC_2021))] #143 obs of 75 variables

#subset to have just the columns that contain "Peak Area" 1:74
#peak_area_cols <- grep("Peak.Area", colnames(OCG_AUC_2021)) 

#the new column names I want to generate will replace the repeating "Peak.Area" names to be "C001" through "C0074" increasing sequentially. 
#new_col_names <- paste0("C", sprintf("%03d", seq_along(peak_area_cols)))

#Rename the columns containing "Peak Area" to compound number
#colnames(OCG_AUC_2021)[peak_area_cols] <- new_col_names #143 obs and 75 variables

#names(OCG_AUC_2021)[names(OCG_AUC_2021) == 'Plant_ID'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

#mdITS_OCG_2021$`Garden Plant ID` <- as.character(mdITS_OCG_2021$`Garden Plant ID`) #as integer so I can combine them

#OCG_AUC_2021 <- left_join(mdITS_OCG_2021, OCG_AUC_2021, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 53 of 96 variables. Some are duplicates.

#Going to remove everything except chem data and plant ID
#OCG_AUC_2021 <- OCG_AUC_2021[,-c(1,3:22)] #53 of 75 variables

#Since there are few plants that seem to have been sampled twice, I am going to remove them for now. I emailed Deb to try and figure this out. Could try and find the avg between then two (below) but not sure that is an appropriate plan of action.

#plant_ID_removed <- c("211", "466", "199", "441", "387", "53", "99", "192", "158", "450") 
#OCG_AUC_2021$`Garden Plant ID`<- as.numeric(OCG_AUC_2021$`Garden Plant ID`)
#str(OCG_AUC_2021)
#OCG_AUC_2021 <- OCG_AUC_2021[!OCG_AUC_2021$`Garden Plant ID` %in% plant_ID_removed, ] #33 obs of 75 variables

# # Convert columns to numeric since it doesnt want to find a mean unless it is numeric
# OCG_AUC_2021 <- as.data.frame(lapply(OCG_AUC_2021, as.numeric))
# str(OCG_AUC_2021)
# 
# # Check for duplicates in the column 'Garden Plant ID'
# duplicates <- OCG_AUC_2021[duplicated(OCG_AUC_2021$Garden.Plant.ID) | duplicated(OCG_AUC_2021$Garden.Plant.ID, fromLast = TRUE), ] #it says there are 10 duplicates.
# 
# # Function to calculate the average of numeric columns
# average_row <- function(OCG_AUC_2021) {
#   OCG_AUC_2021 %>%
#     group_by(Garden.Plant.ID) %>%
#     summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>%
#     ungroup()
# }
# 
# # Create a new row with average values
# new_row <- average_row(OCG_AUC_2021)
# 
# # Print the result
# print(new_row)

##Save this csv so I can re-read in the data with row 1 being plant ID
#write.csv(OCG_AUC_2021, file = "data_csv/OCG_AUC_2021.csv",row.names = FALSE)
#write.csv(mdITS_OCG_2021, file = "data_csv/md_OCG_2021.csv",row.names = FALSE)

#read in metadata
#mdITS_OCG_2021 <- read.csv("data_csv/md_OCG_2021.csv",head=T, check.names = F,stringsAsFactors = T) #43 of 22 variables

## 2021 cleaned GC data read in####
#read in 2021 chemistry data
OCG_AUC_2021 <- read.csv("data_csv/OCG_AUC_2021.csv", row.names = 1) #33 of 74 variables

## 2012/2021 LCMS raw data read in ####
OCG_LCMS_1uL <- read.csv("data_csv/1uL_Injection_Results_LCMS.csv", head=T, check.names = F,stringsAsFactors = T, skip = 0)

sum(duplicated(OCG_LCMS_1uL)) #no duplicates in the dataset

OCG_LCMS_3uL <- read.csv("data_csv/3uL_injection_results_LCMS.csv", head=T, check.names = F,stringsAsFactors = T)

sum(duplicated(OCG_LCMS_3uL)) #no duplicates in the dataset
## LCMS 1ul cleaning####


#Cleaning for PCA####
## 2012 data cleaning ####
colSums(is.na(OCG_AUC_2012)) #checking for null values since pca wont run with NAs. All zeroes for each column which indicates there are no NA values

## 2021 data cleaning ####
#So now the data needs to be cleaned to only contain compounds that occur in more than 20% of the samples (plants). AKA columns need to be removed that don't have at least 20% of the rows containing a number=NA

#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_AUC_2021))
print(na_proportion)

#Now define the threshold of 20% - there are 8 compounds that remain after this
threshold <- 0.2

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_AUC_2021_subset <- OCG_AUC_2021[, columns_to_keep] #33 of 8 variables

# Replace NA values with zeroes
OCG_AUC_2021_subset[is.na(OCG_AUC_2021_subset)] <- 0

colSums(is.na(OCG_AUC_2021_subset)) #all zeroes for each column which indicates there are no NA values

#Scaling####
## 2012 scaling ####
#Compound
data_normalized_2012 <- scale(OCG_AUC_2012) #so this will center the data.... able to code this instead of using excel to do this
head(data_normalized_2012) #looks good

#Plant ID
#transpose the data first to put plant ID as columns and compound as row
OCG_AUC_2012_subset.t <- t(OCG_AUC_2012) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_AUC_2012_subset.t)) #checking for null values

data_normalized_2012.ID <- scale(OCG_AUC_2012_subset.t) 
head(data_normalized_2012.ID) #looks good

## 2021 scaling ####
#Compound correlation plot
data_normalized_2021 <- scale(OCG_AUC_2021_subset) #so this will center the data
head(data_normalized_2021) #mostly negative which is due to the lack of values

#Plant correlation plot
#transpose the data first to put plant ID as columns and compound as row 
OCG_AUC_2021_subset.t <- t(OCG_AUC_2021_subset) # transpose rows and columns

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
#ggcorrplot(corr_matrix_2012.ID)

### 2021 correlation plots ####
#corr plot for compound
corr_matrix_2021 <- cor(data_normalized_2021)
ggcorrplot(corr_matrix_2021)

#corr plot for plant ID
corr_matrix_2021.ID <- cor(data_normalized_2021.ID)
#ggcorrplot(corr_matrix_2021.ID) #looks awful but so did the first one from 2012
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

## 2021 PCA by Plant ID####
#The error message "princomp can only be used with more units than variables" indicates that there are fewer observations (units) than variables in your dataset. Principal Component Analysis (PCA) requires more observations than variables to perform the analysis. This is why it doesnt work.

corr_matrix_2021.ID[is.na(corr_matrix_2021.ID)] <- 0
anyNA(corr_matrix_2021.ID)
data.pca_2021_ID <- princomp(corr_matrix_2021.ID)
summary(data.pca_2021_ID)
dim(corr_matrix_2021.ID) #33 x 33

data.pca_2021_ID$loadings[, 1:2] 

fviz_eig(data.pca_2021_ID, addlabels = TRUE) #scree plot 

# Graph of the variables
fviz_pca_var(data.pca_2021_ID, col.var = "black") #pca

#Contribution of each compound
fviz_cos2(data.pca_2021_ID, choice = "var", axes = 1:2)

#Biplot combined with cos2
fviz_pca_var(data.pca_2021_ID, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

## NMDS plots####
## 2012 NMDS by Compound ####
mdITS_OCG_2012$Subspecies <- as.factor(mdITS_OCG_2012$Subspecies)
mdITS_OCG_2012$Subsp_ploidy <- as.factor(mdITS_OCG_2012$Subsp_ploidy)
mdITS_OCG_2012$Year <- as.factor(mdITS_OCG_2012$Year)
mdITS_OCG_2012$Ploidy <- as.factor(mdITS_OCG_2012$Ploidy)

#turn abundance data frame into a matrix
m_OCG_AUC_2012_subset = as.matrix(OCG_AUC_2012)

set.seed(65)
#OCG_AUC_2012_subset.nmds <- metaMDS(t(m_OCG_AUC_2012_subset), trymax=500) #solution reached.
#save(OCG_AUC_2012_subset.nmds, file = "nmds/OCG_AUC_2012_subset.nmds.rda")
load("nmds/OCG_AUC_2012_subset.nmds.rda")

ordiplot(OCG_AUC_2012_subset.nmds, type = "t",display = "sites",cex = .6) 

## 2012 NMDS by plant ID ####
OCG_AUC_2012_subset.t <- t(OCG_AUC_2012)
m_OCG_AUC_2012_subset.t = as.matrix(OCG_AUC_2012_subset.t)

set.seed(80)
#Plant_ID_OCG_AUC_2012.nmds <- metaMDS(t(m_OCG_AUC_2012_subset.t), trymax=500) #solution reached.
#save(Plant_ID_OCG_AUC_2012.nmds, file = "nmds/Plant_ID_OCG_AUC_2012.nmds")
#load("nmds/Plant_ID_OCG_AUC_2012.nmds.rda")

Plant_ID_OCG_AUC_2012.nmds

ordiplot(Plant_ID_OCG_AUC_2012.nmds, type = "t",display = "sites",cex = .7)

### Ploidy NMDS####
#mdITS_OCG_2012$Ploidy<- droplevels(mdITS_OCG_2012$Ploidy)
levels(mdITS_OCG_2012$Ploidy) #Always smart to check levels before plotting. If there are too many or too few levels, the nmds plots wont line up with the ordispider function.

plot(Plant_ID_OCG_AUC_2012.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Chem comp of 2012 plant by ploidy", 
     col= c("red","blue")[mdITS_OCG_2012$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(Plant_ID_OCG_AUC_2012.nmds,groups = mdITS_OCG_2012$Ploidy, show.groups = "2n", col = "red")
ordispider(Plant_ID_OCG_AUC_2012.nmds,groups = mdITS_OCG_2012$Ploidy, show.groups = "4n", col = "blue")

#### PERMANOVA for 2012 ploidy ####
OCG_AUC_2012_ploidy <- adonis2(OCG_AUC_2012 ~ mdITS_OCG_2012$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_ploidy #ploidy significant 0.033

### Subspecies NMDS ####
#mdITS_OCG_2012$Subspecies<- droplevels(mdITS_OCG_2012$Subspecies)
levels(mdITS_OCG_2012$Subspecies)

plot(Plant_ID_OCG_AUC_2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="Plant chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS_OCG_2012$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(Plant_ID_OCG_AUC_2012.nmds,groups = mdITS_OCG_2012$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(Plant_ID_OCG_AUC_2012.nmds,groups = mdITS_OCG_2012$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(Plant_ID_OCG_AUC_2012.nmds,groups = mdITS_OCG_2012$Subspecies, show.groups = "W", col = "goldenrod")

#### PERMANOVA & pairwaise adonis for subspecies ####
OCG_AUC_2012_subsp <- adonis2(OCG_AUC_2012 ~ mdITS_OCG_2012$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_subsp #subspecies significant 0.001

#pairwiseadonis
OCG_AUC_2012_subsp.pw <- pairwise.adonis(OCG_AUC_2012, mdITS_OCG_2012$Subspecies)
OCG_AUC_2012_subsp.pw#T vs V= 0.003, T vs W= 0.111, and W vs V= 0.003.

### Subspecies ploidy NMDS ####
#mdITS_OCG_2012$Subsp_ploidy<- droplevels(mdITS_OCG_2012$Subsp_ploidy)
levels(mdITS_OCG_2012$Subsp_ploidy)

plot(Plant_ID_OCG_AUC_2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="Sagebrush chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS_OCG_2012$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(Plant_ID_OCG_AUC_2012.nmds,groups = mdITS_OCG_2012$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(Plant_ID_OCG_AUC_2012.nmds,groups = mdITS_OCG_2012$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(Plant_ID_OCG_AUC_2012.nmds,groups = mdITS_OCG_2012$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(Plant_ID_OCG_AUC_2012.nmds,groups = mdITS_OCG_2012$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(Plant_ID_OCG_AUC_2012.nmds,groups = mdITS_OCG_2012$Subsp_ploidy, show.groups = "W_4n", col = "purple")

#### PERMANOVAS and pairwise adonis for subspecies ploidy ####
OCG_AUC_2012_subspploidy <- adonis2(OCG_AUC_2012 ~ mdITS_OCG_2012$Subsp_ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_subspploidy #subspecies ploidy significant 0.001

#pairwiseadonis
OCG_AUC_2012_subsp_ploidy.pw <- pairwise.adonis(OCG_AUC_2012, mdITS_OCG_2012$Subsp_ploidy)
OCG_AUC_2012_subsp_ploidy.pw #T_4n vs T_2n=0.02, T_4n vs V_2n = 0.05, T_2n vs V_2n = 0.01, T_2n vs V_4n = 0.01, T_2n vs W_4n = 0.01, V_4n vs W_4n = 0.06, and V_2n vs W_4n= 0.01.

## 2021 NMDS by Compound####
mdITS_OCG_2021$Subspecies <- as.factor(mdITS_OCG_2021$Subspecies)
mdITS_OCG_2021$Subsp_ploidy <- as.factor(mdITS_OCG_2021$Subsp_ploidy)
mdITS_OCG_2021$Year <- as.factor(mdITS_OCG_2021$Year)
mdITS_OCG_2021$Ploidy <- as.factor(mdITS_OCG_2021$Ploidy)

#turn abundance data frame into a matrix
m_OCG_AUC_2021_subset= as.matrix(OCG_AUC_2021_subset)

set.seed(82)
#OCG_AUC_2021_subset_filtered.nmds <- metaMDS(t(m_OCG_AUC_2021_subset), trymax=500) #solution reached.
#save(OCG_AUC_2021_subset_filtered.nmds, file = "nmds/OCG_AUC_2021_subset_filtered.nmds.rda")
load("nmds/OCG_AUC_2021_subset_filtered.nmds.rda")

#Ordiplot by compound
ordiplot(OCG_AUC_2021_subset_filtered.nmds, type = "t",display = "sites",cex = .6) 

ordiplot(OCG_AUC_2021_subset_filtered.nmds, type = "t",display = "species",cex = .6) 

## 2021 NMDS by Plant ID isn't working####
#turn abundance data frame into a matrix
OCG_AUC_2021_subset_filtered.t <- t(OCG_AUC_2021) #changed this fro the subset to just see if it would run at all. it doesnt

OCG_AUC_2021_subset_filtered.t = as.data.frame(OCG_AUC_2021_subset_filtered.t)
OCG_AUC_2021_subset_filtered.t <- na.omit(OCG_AUC_2021_subset_filtered.t)
OCG_AUC_2021_subset_filtered.t <- OCG_AUC_2021_subset_filtered.t[complete.cases(OCG_AUC_2021_subset_filtered.t), ]

m_OCG_AUC_2021_subset_filtered.t = as.matrix(OCG_AUC_2021_subset_filtered.t)

#colSums(is.na(m_OCG_AUC_2021_subset_filtered.t)) 
#m_OCG_AUC_2021_subset_filtered.t[is.na(m_OCG_AUC_2021_subset_filtered.t)] <- 0

set.seed(4)
Plant_ID_OCG_AUC_2021.nmds <- metaMDS(t(m_OCG_AUC_2021_subset_filtered.t), trymax=500) #
# save(Plant_ID_OCG_AUC.nmds, file = "nmds/Plant_ID_OCG_AUC.nmds.rda")
# load("nmds/Plant_ID_OCG_AUC.nmds.rda")
# 
# #Ordiplot for plant ID
# ordiplot(Plant_ID_OCG_AUC.nmds, type = "t",display = "sites",cex = .7)


# Calculate Manhattan distance: Computes the sum of the absolute differences between coordinates of corresponding points. Suitable for continuous variables. Not affected by empty rows. Less sensitive to outliers compared to Euclidean distance.
manhattan_distance <- dist(m_OCG_AUC_2021_subset_filtered.t, method = "manhattan")
manhattan_distance
manhattan_distance <- sum(abs(m_OCG_AUC_2021_subset_filtered.t))

barplot(manhattan_distance, names.arg = "Manhattan Distance", main = "Manhattan Distance", ylab = "Distance")



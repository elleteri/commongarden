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
#if (!require("xcms"))  {install.packages("xcms"); require("xcms")}
# install.packages("devtools") 
# devtools::install_github("mottensmann/GCalignR", build_vignettes = TRUE) 
#library("GCalignR") 

# Read data in#### 
##Metadata read in#####
mdITS_OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) 
#139 of 22 variables

#subset the md to only have observations from 2012 to avoid duplicates
mdITS_OCG_2012 <- subset(mdITS_OCG, mdITS_OCG$Year=="2012") 
#96 observations and 22 variables

#subset the md to only have observations from 2021 to avoid duplicates
mdITS_OCG_2021 <- subset(mdITS_OCG, mdITS_OCG$Year=="2021") 
#43 observations and 22 variables

## 2012 GC raw data read in####
#read in 2012 GC chemistry data
OCG_AUC_2012 <- read.csv("data_csv/OCG_AUC_2012.csv", head=T, row.names = 1,check.names = F,stringsAsFactors = T) #This data has AUC. 167 of 52 variables

## 2012 GC Cleaning ####
names(OCG_AUC_2012)[names(OCG_AUC_2012) == 'garden_Plant_ID'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

OCG_AUC_2012$`Garden Plant ID` <- as.integer(OCG_AUC_2012$`Garden Plant ID`) #as integer so I can combine them

OCG_AUC_2012 <- merge(mdITS_OCG_2012, OCG_AUC_2012, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 96 of 73 variables

OCG_AUC_2012 <- OCG_AUC_2012[,-c(2:24)] #removing everything except area under the curve for each compound and garden plant ID number
#96 obs of 50 variables

##Save this csv so I can re-read in the data with row 1 being plant ID
write.csv(OCG_AUC_2012, file = "data_csv/OCG_AUC_2012.csv",row.names = FALSE)

## 2012 cleaned GC data read in####
#read in cleaned 2012 GC data with row 1 being plant ID
OCG_AUC_2012 <- read.csv("data_csv/OCG_AUC_2012.csv", row.names = 1) 
#96 of 49 variables

## 2021 GC raw data read in####
#read in 2021 GC chemistry data#
OCG_AUC_2021 <- read.csv("data_csv/OCG_2021_GCData.csv", head=T, skip = 1) #100 of 225 variables

## 2021 GC Cleaning ####
sum(duplicated(OCG_AUC_2021$Plant_ID)) #3 duplicates (100 and then cocktails and blanks)

#head(OCG_AUC_2021)

#sample 100 was ionized suring the first run but then we ran out of gas and the second run is the 100 we actually want to keep so we can remove the first row
OCG_AUC_2021 <- OCG_AUC_2021[-1, ]  # Remove the first row #99 of 225 variables

# Remove rows where "Blank" or "Cocktail" or 429 or 333 appear in Plant ID column 
OCG_AUC_2021 <- OCG_AUC_2021[!(apply(OCG_AUC_2021, 1, function(row) any(grepl("^Blank_|^Cocktail_|^333|^429|^BLANK|Cocktail ", row)))), ] #87 obs of 225 variables

#subset to only columns that contain "Peak Area" and "Plant ID". This removes "RT". 
OCG_AUC_2021 <- OCG_AUC_2021 [, grepl("Peak.Area|Plant_ID", colnames(OCG_AUC_2021))] #87 obs of 149 variables

#subset to remove columns that contain Peak Area Percent and just keep Peak Area.
OCG_AUC_2021 <- OCG_AUC_2021 [, !grepl("Peak.Area.Percent", colnames(OCG_AUC_2021))] #87 obs of 75 variables

#subset to have just the columns that contain "Peak Area" 1:74
peak_area_cols <- grep("Peak.Area", colnames(OCG_AUC_2021)) 

#the new column names I want to generate will replace the repeating "Peak.Area" names to be "C001" through "C0074" increasing sequentially. 
new_col_names <- paste0("C", sprintf("%03d", seq_along(peak_area_cols)))

#Rename the columns containing "Peak Area" to compound number
colnames(OCG_AUC_2021)[peak_area_cols] <- new_col_names #87 obs and 75 variables

names(OCG_AUC_2021)[names(OCG_AUC_2021) == 'Plant_ID'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

mdITS_OCG_2021$`Garden Plant ID` <- as.character(mdITS_OCG_2021$`Garden Plant ID`) #as character so I can combine them

OCG_AUC_2021 <- merge(mdITS_OCG_2021, OCG_AUC_2021, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 39 obs of 96 variables

#Going to remove everything except chem data and plant ID
OCG_AUC_2021 <- OCG_AUC_2021[,-c(2:22)] #39 obs of 75 variables

##Save this csv so I can re-read in the data with row 1 being plant ID
write.csv(OCG_AUC_2021, file = "data_csv/OCG_AUC_2021.csv",row.names = FALSE)

## 2021 cleaned GC data read in####
#read in 2021 chemistry data
OCG_AUC_2021 <- read.csv("data_csv/OCG_AUC_2021.csv", row.names = 1) 
#39 obs of 74 variables
# Replace NA values with zeroes
OCG_AUC_2021[is.na(OCG_AUC_2021)] <- 0

colSums(is.na(OCG_AUC_2021)) #all zeroes for each column which indicates there are no NA values

## 2012/2021 LCMS raw data read in ####
OCG_LCMS_1uL <- read.csv("data_csv/1uL_Injection_Results_LCMS.csv", head=T, check.names = F,stringsAsFactors = T, skip = 1) #120 of 930 variables

sum(duplicated(OCG_LCMS_1uL$Plant_ID)) #one duplicate: plant #273

OCG_LCMS_3uL <- read.csv("data_csv/3uL_injection_results_LCMS.csv", head=T, check.names = F,stringsAsFactors = T, skip = 1) #120 of 929 variables

sum(duplicated(OCG_LCMS_3uL$Plant_ID)) #no duplicates in the dataset

## LCMS 1ul cleaning####
head(OCG_LCMS_1uL)

#remove row 8 and 61 since that is the 273 sample
OCG_LCMS_1uL <- OCG_LCMS_1uL[-8, ]  # Remove the first row 8 #119 of 930
OCG_LCMS_1uL <- OCG_LCMS_1uL[-61, ]  # Remove the first row 61 #118 of 930

#subset to only columns that contain "Peak Area" and "Plant ID". This removes "RT". 
OCG_LCMS_1uL <- OCG_LCMS_1uL [, grepl("Peak.Area|Plant_ID", colnames(OCG_LCMS_1uL))] #118 obs of 619 variables

#subset to remove columns that contain Peak Area Percent and just keep Peak Area.
OCG_LCMS_1uL <- OCG_LCMS_1uL [, !grepl("Peak.Area.Percent", colnames(OCG_LCMS_1uL))] #118 obs of 310 variables

#now i'm going to have to rename all of the peak areas by their compound name? unless there is a better way to retain those when reading in the data. The compounds are in order 10, 100-109, 11, 110-119, 12, 120-129, 13.... etc. but not all the way throughout so I think I will just create a key until I have found an alternate way to look at the compound name

#subset to have just the columns that contain "Peak Area" 1:309
peak_area_cols <- grep("Peak.Area", colnames(OCG_LCMS_1uL)) 

#the new column names I want to generate will replace the repeating "Peak.Area" names to be "C001" through "C0074" increasing sequentially. 
new_col_names <- paste0("C", sprintf("%03d", seq_along(OCG_LCMS_1uL)))

#Rename the columns containing "Peak Area" to compound number
colnames(OCG_LCMS_1uL)[peak_area_cols] <- new_col_names #118 obs of 310 variables 

#column names are now C001-C309

# Splitting the Plant_ID column and creating new columns
OCG_LCMS_1uL <- OCG_LCMS_1uL %>%
  separate(Plant_ID, into = c("Plant_ID_Number", "Year"), sep = "_")

# Convert 'Year' column to integer
OCG_LCMS_1uL$Year <- as.integer(OCG_LCMS_1uL$Year)

#dataframe now includes the plant ID number and year as their own columns to match with the metadata. 

names(OCG_LCMS_1uL)[names(OCG_LCMS_1uL) == 'Plant_ID_Number'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

#subset for the two years
#2012#
#subset the md to only have observations from 2012 to avoid duplicates
OCG_LCMS_1uL_2012 <- subset(OCG_LCMS_1uL, OCG_LCMS_1uL$Year=="2012") #45 observations and 311 variables

mdITS_OCG_2012$`Garden Plant ID` <- as.character(mdITS_OCG_2012$`Garden Plant ID`) #as character so I can combine them

OCG_LCMS_1uL_2012_subset <- merge(mdITS_OCG_2012, OCG_LCMS_1uL_2012, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadata of OCG to the samples we have. that only leaves 9 obs of 332 variables

OCG_LCMS_1uL_2012_subset <- OCG_LCMS_1uL_2012_subset[,-c(2:23)] #9 obs of 310 variables

#save the csv for 2012
write.csv(OCG_LCMS_1uL_2012_subset, file = "data_csv/OCG_LCMS_1uL_2012.csv",row.names = FALSE)

#rownames(OCG_LCMS_1uL_2012_subset) <- OCG_LCMS_1uL_2012_subset[, 1] #this will set row.names equal to one since it didnt change it when I read it is
#OCG_LCMS_1uL_2012_subset <- OCG_LCMS_1uL_2012_subset[, -1]  # Remove the first column after setting row names #9 obs of 309 variables

#2021#
#subset the md to only have observations from 2021 to avoid duplicates
OCG_LCMS_1uL_2021 <- subset(OCG_LCMS_1uL, OCG_LCMS_1uL$Year=="2021") #72 observations and 311 variables

mdITS_OCG_2021$`Garden Plant ID` <- as.character(mdITS_OCG_2021$`Garden Plant ID`) #as integer so I can combine them

OCG_LCMS_1uL_2021_subset <- merge(mdITS_OCG_2021, OCG_LCMS_1uL_2021, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 

#Going to remove everything except chem data and plant ID
OCG_LCMS_1uL_2021_subset <- OCG_LCMS_1uL_2021_subset[,-c(2:23)] #38 obs of 310 variables

##Save this csv so I can re-read in the data with row 1 being plant ID
write.csv(OCG_LCMS_1uL_2021_subset, file = "data_csv/OCG_LCMS_1uL_2021.csv",row.names = FALSE)

#rownames(OCG_LCMS_1uL_2021_subset) <- OCG_LCMS_1uL_2021_subset[, 1] #this will set row.names equal to one since it didnt change it when I read it is
#OCG_LCMS_1uL_2021_subset <- OCG_LCMS_1uL_2021_subset[, -1]  # Remove the first column after setting row names #38 of 309 variables

## Cleaned 1uL LCMS data read in ####
#2012
OCG_LCMS_1uL_2012 <- read.csv("data_csv/OCG_LCMS_1uL_2012.csv", head=T)
rownames(OCG_LCMS_1uL_2012) <- OCG_LCMS_1uL_2012[, 1] #this will set row.names equal to one since it didnt change it when I read it is
OCG_LCMS_1uL_2012 <- OCG_LCMS_1uL_2012[, -1]  # Remove the first column after setting row names 
#9 obs of 309 variables

#2021
OCG_LCMS_1uL_2021 <- read.csv("data_csv/OCG_LCMS_1uL_2021.csv", head=T)
rownames(OCG_LCMS_1uL_2021) <- OCG_LCMS_1uL_2021[, 1] #this will set row.names equal to one since it didnt change it when I read it is
OCG_LCMS_1uL_2021 <- OCG_LCMS_1uL_2021[, -1]  # Remove the first column after setting row names
#38 obs of 309 variables

## LCMS 3ul cleaning####
head(OCG_LCMS_3uL)

#subset to only columns that contain "Peak Area" and "Plant ID". This removes "RT". 
OCG_LCMS_3uL <- OCG_LCMS_3uL [, grepl("Peak.Area|Plant_ID", colnames(OCG_LCMS_3uL))] #120 obs of 617 variables

#subset to remove columns that contain Peak Area Percent and just keep Peak Area.
OCG_LCMS_3uL <- OCG_LCMS_3uL [, !grepl("Peak.Area.Percent", colnames(OCG_LCMS_3uL))] #118 obs of 309 variables

#now i'm going to have to rename all of the peak areas by their compound name? unless there is a better way to retain those when reading in the data. The compounds are in order 10, 100-109, 11, 110-119, 12, 120-129, 13.... etc. but not all the way throughout so I think I will just create a key until I have found an alternate way to look at the compound name

#subset to have just the columns that contain "Peak Area" 1:308
peak_area_cols <- grep("Peak.Area", colnames(OCG_LCMS_3uL)) 

#the new column names I want to generate will replace the repeating "Peak.Area" names to be "C001" through "C0074" increasing sequentially. 
new_col_names <- paste0("C", sprintf("%03d", seq_along(OCG_LCMS_3uL)))

#Rename the columns containing "Peak Area" to compound number
colnames(OCG_LCMS_3uL)[peak_area_cols] <- new_col_names #120 obs of 309 variables 

#column names are now C001-C309

# Splitting the Plant_ID column and creating new columns
OCG_LCMS_3uL <- OCG_LCMS_3uL %>%
  separate(Plant_ID, into = c("Plant_ID_Number", "Year"), sep = "_")

# Convert 'Year' column to integer
OCG_LCMS_3uL$Year <- as.integer(OCG_LCMS_3uL$Year)

#dataframe now includes the plant ID number and year as their own columns to match with the metadata. 

names(OCG_LCMS_3uL)[names(OCG_LCMS_3uL) == 'Plant_ID_Number'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

#subset for the two years
#2012#
#subset the md to only have observations from 2012 to avoid duplicates
OCG_LCMS_3uL_2012 <- subset(OCG_LCMS_3uL, OCG_LCMS_3uL$Year=="2012") #46 observations and 310 variables

OCG_LCMS_3uL_2012_subset <- merge(mdITS_OCG_2012, OCG_LCMS_3uL_2012, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadata of OCG to the samples we have. that only leaves 9 obs of 332 variables

OCG_LCMS_3uL_2012_subset <- OCG_LCMS_3uL_2012_subset[,-c(2:23)] #10 obs of 309 variables

#save the csv for 2012
write.csv(OCG_LCMS_3uL_2012_subset, file = "data_csv/OCG_LCMS_3uL_2012.csv",row.names = FALSE)

#rownames(OCG_LCMS_1uL_2012_subset) <- OCG_LCMS_1uL_2012_subset[, 1] #this will set row.names equal to one since it didnt change it when I read it is
#OCG_LCMS_1uL_2012_subset <- OCG_LCMS_1uL_2012_subset[, -1]  # Remove the first column after setting row names #9 obs of 309 variables

#2021#
#subset the md to only have observations from 2021 to avoid duplicates
OCG_LCMS_3uL_2021 <- subset(OCG_LCMS_3uL, OCG_LCMS_3uL$Year=="2021") #73 observations and 310 variables

OCG_LCMS_3uL_2021_subset <- merge(mdITS_OCG_2021, OCG_LCMS_3uL_2021, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 

#Going to remove everything except chem data and plant ID
OCG_LCMS_3uL_2021_subset <- OCG_LCMS_3uL_2021_subset[,-c(2:23)] #38 obs of 309 variables

##Save this csv so I can re-read in the data with row 1 being plant ID
write.csv(OCG_LCMS_3uL_2021_subset, file = "data_csv/OCG_LCMS_3uL_2021.csv",row.names = FALSE)
## Cleaned 3uL LCMS data read in ####
#2012
OCG_LCMS_3uL_2012 <- read.csv("data_csv/OCG_LCMS_3uL_2012.csv", head=T)
rownames(OCG_LCMS_3uL_2012) <- OCG_LCMS_3uL_2012[, 1] #this will set row.names equal to one since it didnt change it when I read it is
OCG_LCMS_3uL_2012 <- OCG_LCMS_3uL_2012[, -1]  # Remove the first column after setting row names 
#10 obs of 308 variables

#2021
OCG_LCMS_3uL_2021 <- read.csv("data_csv/OCG_LCMS_3uL_2021.csv", head=T)
rownames(OCG_LCMS_3uL_2021) <- OCG_LCMS_3uL_2021[, 1] #this will set row.names equal to one since it didnt change it when I read it is
OCG_LCMS_3uL_2021 <- OCG_LCMS_3uL_2021[, -1]  # Remove the first column after setting row names 
#38 obs of 308 variables

#Cleaning for PCA####
## 2012 GC data cleaning ####
colSums(is.na(OCG_AUC_2012)) #checking for null values since pca wont run with NAs. All zeroes for each column which indicates there are no NA values

## 2021 GC data cleaning ####
#So now the data needs to be cleaned to only contain compounds that occur in more than 20% of the samples (plants). AKA columns need to be removed that don't have at least 20% of the rows containing a number=NA

#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_AUC_2021))
print(na_proportion)

#Now define the threshold of 20% - there are 15 compounds that remain after this
threshold <- 0.2

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_AUC_2021_subset <- OCG_AUC_2021[, columns_to_keep] #39 of 15 variables

# Replace NA values with zeroes
OCG_AUC_2021_subset[is.na(OCG_AUC_2021_subset)] <- 0

colSums(is.na(OCG_AUC_2021_subset)) #all zeroes for each column which indicates there are no NA values

## 2012 1uL LCMS data cleaning ####
colSums(is.na(OCG_LCMS_1uL_2012)) #checking for null values
# Replace NA values with zeroes
OCG_LCMS_1uL_2012[is.na(OCG_LCMS_1uL_2012)] <- 0

## 2021 1uL LCMS data cleaning ####
colSums(is.na(OCG_LCMS_1uL_2021)) #checking for null values
# Replace NA values with zeroes
OCG_LCMS_1uL_2021[is.na(OCG_LCMS_1uL_2021)] <- 0

## 2012 3uL LCMS data cleaning ####
colSums(is.na(OCG_LCMS_3uL_2012)) #checking for null values
# Replace NA values with zeroes
OCG_LCMS_3uL_2012[is.na(OCG_LCMS_3uL_2012)] <- 0

## 2021 3uL LCMS data cleaning ####
colSums(is.na(OCG_LCMS_3uL_2021)) #checking for null values
# Replace NA values with zeroes
OCG_LCMS_3uL_2021[is.na(OCG_LCMS_3uL_2021)] <- 0

#Scaling####
## 2012 GC scaling ####
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

## 2021 GC scaling ####
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

## 2012 1uL LCMS scaling####
data_LCMS1_normalized_2012 <- scale(OCG_LCMS_1uL_2012) #this will center the data
head(data_LCMS1_normalized_2012) 

#Plant correlation plot
#transpose the data first to put plant ID as columns and compound as row 
OCG_LCMS_1uL_2012.t <- t(OCG_LCMS_1uL_2012) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_LCMS_1uL_2012.t)) #checking for null values

data_LCMS1_normalized_2012.ID <- scale(OCG_LCMS_1uL_2012.t) 
head(data_LCMS1_normalized_2012.ID) 

## 2021 1uL LCMS scaling####
data_LCMS1_normalized_2021 <- scale(OCG_LCMS_1uL_2021) #this will center the data
head(data_LCMS1_normalized_2021) 

#Plant correlation plot
#transpose the data first to put plant ID as columns and compound as row 
OCG_LCMS_1uL_2021.t <- t(OCG_LCMS_1uL_2021) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_LCMS_1uL_2021.t)) #checking for null values

data_LCMS1_normalized_2021.ID <- scale(OCG_LCMS_1uL_2021.t) 
head(data_LCMS1_normalized_2012.ID) 

## 2012 3uL LCMS scaling####
data_LCMS3_normalized_2012 <- scale(OCG_LCMS_3uL_2012) #this will center the data
head(data_LCMS3_normalized_2012) 

#Plant correlation plot
#transpose the data first to put plant ID as columns and compound as row 
OCG_LCMS_3uL_2012.t <- t(OCG_LCMS_3uL_2012) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_LCMS_3uL_2012.t)) #checking for null values

data_LCMS3_normalized_2012.ID <- scale(OCG_LCMS_3uL_2012.t) 
head(data_LCMS3_normalized_2012.ID) 

## 2021 3uL LCMS scaling####
data_LCMS3_normalized_2021 <- scale(OCG_LCMS_3uL_2021) #this will center the data
head(data_LCMS3_normalized_2021) 

#Plant correlation plot
#transpose the data first to put plant ID as columns and compound as row 
OCG_LCMS_3uL_2021.t <- t(OCG_LCMS_3uL_2021) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_LCMS_3uL_2021.t)) #checking for null values

data_LCMS3_normalized_2021.ID <- scale(OCG_LCMS_3uL_2021.t) 
head(data_LCMS3_normalized_2021.ID)

# Data visualization####
## Correlation plots####
### 2012 GC correlation plots ####
#corr plot for compound
corr_matrix_2012 <- cor(data_normalized_2012)
ggcorrplot(corr_matrix_2012)

#corr plot for plant ID
corr_matrix_2012.ID <- cor(data_normalized_2012.ID)
#ggcorrplot(corr_matrix_2012.ID)

### 2021 GC correlation plots ####
#corr plot for compound
corr_matrix_2021 <- cor(data_normalized_2021)
ggcorrplot(corr_matrix_2021)

#corr plot for plant ID
corr_matrix_2021.ID <- cor(data_normalized_2021.ID)
#ggcorrplot(corr_matrix_2021.ID)

### 2012 1uL LCMS correlation plots ####
#corr plot for compound
corr_matrix_LCMS1_2012 <- cor(data_LCMS1_normalized_2012)
ggcorrplot(corr_matrix_LCMS1_2012)

#corr plot for plant ID
corr_matrix_LCMS1_2012.ID <- cor(data_LCMS1_normalized_2012.ID)
#ggcorrplot(corr_matrix_LCMS1_2012.ID)

### 2021 1uL LCMS correlation plots ####
#corr plot for compound
corr_matrix_LCMS1_2021 <- cor(data_LCMS1_normalized_2021)
ggcorrplot(corr_matrix_LCMS1_2021)

#corr plot for plant ID
corr_matrix_LCMS1_2021.ID <- cor(data_LCMS1_normalized_2021.ID)
#ggcorrplot(corr_matrix_2021.ID)

### 2012 3uL LCMS correlation plots ####
#corr plot for compound
corr_matrix_LCMS3_2012 <- cor(data_LCMS3_normalized_2012)
ggcorrplot(corr_matrix_LCMS3_2012)

#corr plot for plant ID
corr_matrix_LCMS3_2012.ID <- cor(data_LCMS3_normalized_2012.ID)
#ggcorrplot(corr_matrix_2021.ID)

### 2021 3uL LCMS correlation plots ####
#corr plot for compound
corr_matrix_LCMS3_2021 <- cor(data_LCMS3_normalized_2021)
ggcorrplot(corr_matrix_LCMS3_2021)

#corr plot for plant ID
corr_matrix_LCMS3_2021.ID <- cor(data_LCMS3_normalized_2021.ID)
#ggcorrplot(corr_matrix_2021.ID)

# PCA ####
## 2012 GC PCA by compound####
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

## 2021 GC PCA by compound####
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

## 2012 GC PCA by Plant ID ####
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

## 2021 GC PCA by Plant ID####
#The error message "princomp can only be used with more units than variables" indicates that there are fewer observations (units) than variables in your dataset. Principal Component Analysis (PCA) requires more observations than variables to perform the analysis. This is why it doesnt work.

#corr_matrix_2021.ID[is.na(corr_matrix_2021.ID)] <- 0
#anyNA(corr_matrix_2021.ID)
data.pca_2021_ID <- princomp(corr_matrix_2021.ID)
summary(data.pca_2021_ID)

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

## 2012 1uL LCMS PCA by Plant ID####
data.pca_LCMS1_2012.ID <- princomp(corr_matrix_LCMS1_2012.ID)
summary(data.pca_LCMS1_2012.ID)

data.pca_LCMS1_2012.ID$loadings[, 1:2]

fviz_eig(data.pca_LCMS1_2012.ID, addlabels = TRUE) #scree plot is used to visualize the importance of each principal component and can be used to determine the number of principal components to retain.

#With the biplot, it is possible to visualize the similarities and dissimilarities between the samples, and further shows the impact of each attribute on each of the principal components.

# Graph of the compounds
fviz_pca_var(data.pca_LCMS1_2012.ID, col.var = "black")

#Contribution of each compound
fviz_cos2(data.pca_LCMS1_2012.ID, choice = "var", axes = 1:2)

#Biplot combined with cos2
fviz_pca_var(data.pca_LCMS1_2012.ID, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

## 2021 1uL LCMS PCA by Plant ID####
data.pca_LCMS1_2021.ID <- princomp(corr_matrix_LCMS1_2021.ID)
summary(data.pca_LCMS1_2021.ID)

data.pca_LCMS1_2021.ID$loadings[, 1:2]

fviz_eig(data.pca_LCMS1_2021.ID, addlabels = TRUE) #scree plot is used to visualize the importance of each principal component and can be used to determine the number of principal components to retain.

#With the biplot, it is possible to visualize the similarities and dissimilarities between the samples, and further shows the impact of each attribute on each of the principal components.

# Graph of the compounds
fviz_pca_var(data.pca_LCMS1_2021.ID, col.var = "black")

#Contribution of each compound
fviz_cos2(data.pca_LCMS1_2021.ID, choice = "var", axes = 1:2)

#Biplot combined with cos2
fviz_pca_var(data.pca_LCMS1_2021.ID, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
## 2012 3uL LCMS PCA by Plant ID####
data.pca_LCMS3_2012.ID <- princomp(corr_matrix_LCMS3_2012.ID)
summary(data.pca_LCMS3_2012.ID)

data.pca_LCMS3_2012.ID$loadings[, 1:2]

fviz_eig(data.pca_LCMS3_2012.ID, addlabels = TRUE) #scree plot is used to visualize the importance of each principal component and can be used to determine the number of principal components to retain.

#With the biplot, it is possible to visualize the similarities and dissimilarities between the samples, and further shows the impact of each attribute on each of the principal components.

# Graph of the compounds
fviz_pca_var(data.pca_LCMS3_2012.ID, col.var = "black")

#Contribution of each compound
fviz_cos2(data.pca_LCMS3_2012.ID, choice = "var", axes = 1:2)

#Biplot combined with cos2
fviz_pca_var(data.pca_LCMS3_2012.ID, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

## 2021 3uL LCMS PCA by Plant ID####
data.pca_LCMS3_2021.ID <- princomp(corr_matrix_LCMS3_2021.ID)
summary(data.pca_LCMS3_2021.ID)

data.pca_LCMS3_2021.ID$loadings[, 1:2]

fviz_eig(data.pca_LCMS3_2021.ID, addlabels = TRUE) #scree plot is used to visualize the importance of each principal component and can be used to determine the number of principal components to retain.

#With the biplot, it is possible to visualize the similarities and dissimilarities between the samples, and further shows the impact of each attribute on each of the principal components.

# Graph of the compounds
fviz_pca_var(data.pca_LCMS3_2021.ID, col.var = "black")

#Contribution of each compound
fviz_cos2(data.pca_LCMS3_2021.ID, choice = "var", axes = 1:2)

#Biplot combined with cos2
fviz_pca_var(data.pca_LCMS3_2021.ID, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

## NMDS plots####
## 2012 GC NMDS by Compound ####
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

## 2012 GC NMDS by plant ID ####
OCG_AUC_2012_subset.t <- t(OCG_AUC_2012)
m_OCG_AUC_2012_subset.t = as.matrix(OCG_AUC_2012_subset.t)

set.seed(57)
Plant_ID_OCG_AUC_2012.nmds <- metaMDS(t(m_OCG_AUC_2012_subset.t), trymax=500) #
#save(Plant_ID_OCG_AUC_2012.nmds, file = "nmds/Plant_ID_OCG_AUC_2012.nmds")
load("nmds/Plant_ID_OCG_AUC_2012.nmds.rda")

Plant_ID_OCG_AUC_2012.nmds

ordiplot(Plant_ID_OCG_AUC_2012.nmds, type = "t",display = "sites",cex = .7)

## 2012 1uL LCMS NMDS by plant ID ####
OCG_LCMS_1uL_2012.t <- t(OCG_LCMS_1uL_2012)
m_OCG_LCMS_1uL_2012.t = as.matrix(OCG_LCMS_1uL_2012.t)

set.seed(53)
OCG_LCMS_1uL_2012_ID.nmds <- metaMDS(t(m_OCG_LCMS_1uL_2012.t), trymax=500) #insufficient data for nmds even though solutiin technically reached

ordiplot(OCG_LCMS_1uL_2012_ID.nmds, type = "t",display = "sites",cex = .7)

## 2021 1uL LCMS NMDS by plant ID ####
OCG_LCMS_1uL_2021.t <- t(OCG_LCMS_1uL_2021)
m_OCG_LCMS_1uL_2021.t = as.matrix(OCG_LCMS_1uL_2021.t)

#set.seed(453)
#OCG_LCMS_1uL_2021_ID.nmds <- metaMDS(t(m_OCG_LCMS_1uL_2021.t), trymax=500) #solution reached
#save(OCG_LCMS_1uL_2021_ID.nmds, file = "nmds/OCG_LCMS_1uL_2021_ID.nmds.rda")
load("nmds/OCG_LCMS_1uL_2021_ID.nmds.rda")

ordiplot(OCG_LCMS_1uL_2021_ID.nmds, type = "t",display = "sites",cex = .7)

## 2012 3uL LCMS NMDS by plant ID ####
OCG_LCMS_3uL_2012.t <- t(OCG_LCMS_3uL_2012)
m_OCG_LCMS_3uL_2012.t = as.matrix(OCG_LCMS_3uL_2012.t)

set.seed(86)
OCG_LCMS_3uL_2012_ID.nmds <- metaMDS(t(m_OCG_LCMS_3uL_2012.t), trymax=500) #insufficient data for nmds even though solution technically reached

ordiplot(OCG_LCMS_3uL_2012_ID.nmds, type = "t",display = "sites",cex = .7)

## 2021 3uL LCMS NMDS by plant ID ####
OCG_LCMS_3uL_2021.t <- t(OCG_LCMS_3uL_2021)
m_OCG_LCMS_3uL_2021.t = as.matrix(OCG_LCMS_3uL_2021.t)

# set.seed(69)
# OCG_LCMS_3uL_2021_ID.nmds <- metaMDS(t(m_OCG_LCMS_3uL_2021.t), trymax=500) #solution reached
# save(OCG_LCMS_3uL_2021_ID.nmds, file = "nmds/OCG_LCMS_3uL_2021_ID.nmds.rda")
OCG_LCMS_3uL_2021_ID.nmds <- load("nmds/OCG_LCMS_3uL_2021_ID.nmds.rda")

ordiplot(OCG_LCMS_3uL_2021_ID.nmds, type = "t",display = "sites",cex = .7)

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
OCG_AUC_2012_ploidy #ploidy is not significant 0.32

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
OCG_AUC_2012_subsp #subspecies is not significant 0.587

#pairwiseadonis
OCG_AUC_2012_subsp.pw <- pairwise.adonis(OCG_AUC_2012, mdITS_OCG_2012$Subspecies)
OCG_AUC_2012_subsp.pw# no significance between subspecies

### Subspecies ploidy NMDS ####
#mdITS_OCG_2012$Subsp_ploidy<- droplevels(mdITS_OCG_2012$Subsp_ploidy)
levels(mdITS_OCG_2012$Subsp_ploidy)

plot(Plant_ID_OCG_AUC_2012.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC chemistry by subspecies and ploidy", 
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
OCG_AUC_2012_subspploidy #subspecies ploidy not significant

#pairwiseadonis
OCG_AUC_2012_subsp_ploidy.pw <- pairwise.adonis(OCG_AUC_2012, mdITS_OCG_2012$Subsp_ploidy)
OCG_AUC_2012_subsp_ploidy.pw #

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

## 2021 NMDS by Plant ID####
#turn abundance data frame into a matrix
OCG_AUC_2021_subset_filtered.t <- t(OCG_AUC_2021) 

set.seed(4)
Plant_ID_OCG_AUC_2021.nmds <- metaMDS(t(m_OCG_AUC_2021_subset_filtered.t), trymax=500) # warning
save(Plant_ID_OCG_AUC_2021.nmds, file = "nmds/Plant_ID_OCG_AUC_2021.nmds.rda")
load("nmds/Plant_ID_OCG_AUC_2021.nmds.rda")

#Ordiplot for plant ID
ordiplot(Plant_ID_OCG_AUC_2021.nmds, type = "t",display = "sites",cex = .7)

# Calculate Manhattan distance: Computes the sum of the absolute differences between coordinates of corresponding points. Suitable for continuous variables. Not affected by empty rows. Less sensitive to outliers compared to Euclidean distance.



#### PERMANOVA for 2021 ploidy ####
mdITS_OCG_2021$Subspecies <- as.factor(mdITS_OCG_2021$Subspecies)
mdITS_OCG_2021$Subsp_ploidy <- as.factor(mdITS_OCG_2021$Subsp_ploidy)
mdITS_OCG_2021$Year <- as.factor(mdITS_OCG_2021$Year)
mdITS_OCG_2021$Ploidy <- as.factor(mdITS_OCG_2021$Ploidy)

dim(OCG_AUC_2021_subset)
dim(mdITS_OCG_2021$Ploidy)

OCG_AUC_2021_ploidy <- adonis2(OCG_AUC_2021_subset ~ mdITS_OCG_2021$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2021_ploidy #ploidy is not significant 0.33

#### PERMANOVA & pairwaise adonis for subspecies ####
OCG_AUC_2021_subsp <- adonis2(OCG_AUC_2021 ~ mdITS_OCG_2021$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2021_subsp #subspecies is not significant 0.587

#pairwiseadonis
OCG_AUC_2021_subsp.pw <- pairwise.adonis(OCG_AUC_2021, mdITS_OCG_2021$Subspecies)
OCG_AUC_2021_subsp.pw# no significance between subspecies




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
library(ggfortify)
library(cluster)
#if (!require("xcms"))  {install.packages("xcms"); require("xcms")}
# install.packages("devtools") 
# devtools::install_github("mottensmann/GCalignR", build_vignettes = TRUE) 
#library("GCalignR") 

# Read data in#### 
##Metadata read in#####
mdITS <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #set to correct file path 505 obs of 16 variables.
mdITS <- mdITS[order(row.names(mdITS)),] # order samples alphabetically

### Subsetting to just the plant in the common garden (OCG)#
mdITS.OCG <- subset(mdITS, mdITS$Project=="OCG") #subsetting to just the Orchard common garden plants #246 of 16 variables

### Remove duplicates from asv#
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1', 'ORT.1.5_2012')
mdITS.OCG <- mdITS.OCG[!rownames(mdITS.OCG) %in% rows_to_remove, ]

## Remove negative control
mdITS.OCG <- mdITS.OCG[!(row.names(mdITS.OCG) == "NEG_8-28-21"),]
mdITS.OCG <- mdITS.OCG[!(row.names(mdITS.OCG) == "NEG_10-2-20"),]

## Remove MTW.3.7.R_2012
mdITS.OCG <- mdITS.OCG[!(row.names(mdITS.OCG) == "MTW.3.7.R_2012"),] #we arent sure what the R represents. 236 of 16 var

#make variables factor to plot
mdITS.OCG$Ploidy <- as.factor(mdITS.OCG$Ploidy)
mdITS.OCG$Subspecies <- as.factor(mdITS.OCG$Subspecies)
mdITS.OCG$Subsp_ploidy <- as.factor(mdITS.OCG$Subsp_ploidy)
mdITS.OCG$Year <- as.factor(mdITS.OCG$Year)
mdITS.OCG$Plant <- as.factor(mdITS.OCG$Plant)
mdITS.OCG <- droplevels(mdITS.OCG)
levels(mdITS.OCG$Subspecies)


#mdITS_OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) 
#154 of 15 variables

#subset the md to only have observations from 2012 to avoid duplicates
mdITS_OCG_2012 <- subset(mdITS.OCG, mdITS.OCG$Year=="2012") 
#160
mdITS_OCG_2012$`Garden Plant ID` <- as.factor(mdITS_OCG_2012$`Garden Plant ID`) #as factor so I can combine them
mdITS_OCG_2012 <- droplevels(mdITS_OCG_2012)
str(mdITS_OCG_2012)

#subset the md to only have observations from 2021 to avoid duplicates
mdITS_OCG_2021 <- subset(mdITS.OCG, mdITS.OCG$Year=="2021") 
#76
mdITS_OCG_2021 <- droplevels(mdITS_OCG_2021)
str(mdITS_OCG_2021)

## 2012 GC raw data read in####
#read in 2012 GC chemistry data
OCG_AUC_2012 <- read.csv("data_csv/OCG_2012_GC.csv", head=T,check.names = F,stringsAsFactors = T, skip = 1) #183 obs of 221 var

## 2012 GC Cleaning ####
# Remove the second column since it is empty
OCG_AUC_2012 <- OCG_AUC_2012[, -2]
names(OCG_AUC_2012)[1] <- "Plant_ID"

sum(duplicated(OCG_AUC_2012$Plant_ID)) #4 duplicates (control, cocktails, and blanks)

# Remove rows where "empty", "Empty", "Ct, "CT", "M", "w"  occur
OCG_AUC_2012 <- OCG_AUC_2012[!(apply(OCG_AUC_2012, 1, function(row) any(grepl("^empty|^Empty|^Ct|^CT |^M|w", row)))), ] #167 of 220 variables

sum(duplicated(OCG_AUC_2012$Plant_ID)) # 0 duplicates
#head(OCG_AUC_2021)

#subset to only columns that contain "Peak Area" and "Plant ID". This removes "RT". 
OCG_AUC_2012 <- OCG_AUC_2012 [, grepl("Peak.Area|Plant_ID", colnames(OCG_AUC_2012))] #167 obs of 147 variables

#subset to remove columns that contain Peak Area Percent and just keep Peak Area.
OCG_AUC_2012 <- OCG_AUC_2012 [, !grepl("Peak.Area.Percent", colnames(OCG_AUC_2012))] #167 obs of 74 variables

#subset to have just the columns that contain "Peak Area" 1:73 #shifted by one compared to 2021 since there was no compound 1. should fix. 
peak_area_cols <- grep("Peak.Area", colnames(OCG_AUC_2012)) 

#the new column names I want to generate will replace the repeating "Peak.Area" names to be "C001" through "C0073" increasing sequentially. 
new_col_names <- paste0("C", sprintf("%03d", seq_along(peak_area_cols)))

#Rename the columns containing "Peak Area" to compound number
colnames(OCG_AUC_2012)[peak_area_cols] <- new_col_names #167 of 74 variables

names(OCG_AUC_2012)[names(OCG_AUC_2012) == 'Plant_ID'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

OCG_AUC_2012 <- merge(mdITS_OCG_2012, OCG_AUC_2012, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 157 obs of 89

OCG_AUC_2012 <- OCG_AUC_2012[,-c(1:15)] #removing everything except area under the curve for each compound and garden plant ID number
#157 obs of 74 variables

##Save this csv so I can re-read in the data with row 1 being plant ID description
write.csv(OCG_AUC_2012, file = "data_csv/OCG_AUC_2012_cleaned.csv",row.names = FALSE)

## 2012 cleaned GC data read in####
#read in cleaned 2012 GC data with row 1 being plant ID
OCG_AUC_2012 <- read.csv("data_csv/OCG_AUC_2012_cleaned.csv", row.names = 1) 
#157 obs of 73 variables
OCG_AUC_2012 <- OCG_AUC_2012[order(row.names(OCG_AUC_2012)),] # order samples alphabetically

OCG_AUC_2012 <- subset(OCG_AUC_2012, row.names(OCG_AUC_2012) %in% row.names(mdITS)) #148 of 73 variables

## 2021 GC raw data read in####
#read in 2021 GC chemistry data#
OCG_AUC_2021 <- read.csv("data_csv/OCG_2021_GCData.csv", head=T, skip = 1) #100 of 225 variables

## 2021 GC Cleaning ####
sum(duplicated(OCG_AUC_2021$Plant_ID)) #3 duplicates (100 and then cocktails and blanks)

#head(OCG_AUC_2021)

#sample 100 was ionized during the first run but then we ran out of gas and the second run is the 100 we actually want to keep so we can remove the first row
OCG_AUC_2021 <- OCG_AUC_2021[-1, ]  # Remove the first row #99 of 225 variables

# Remove rows where "Blank" or "Cocktail" or 429 or 333 appear in Plant ID column 
OCG_AUC_2021 <- OCG_AUC_2021[!(apply(OCG_AUC_2021, 1, function(row) any(grepl("^Blank_|^Cocktail_|^333|^429|^BLANK|Cocktail ", row)))), ] #87 obs of 225 variables

#subset to only columns that contain "Peak Area" and "Plant ID". This removes "RT". 
OCG_AUC_2021 <- OCG_AUC_2021 [, grepl("Peak.Area|Plant_ID", colnames(OCG_AUC_2021))] #87 obs of 149 variables

#subset to remove columns that contain Peak Area Percent and just keep Peak Area.
OCG_AUC_2021 <- OCG_AUC_2021 [, !grepl("Peak.Area.Percent", colnames(OCG_AUC_2021))] #87 obs of 75 variables

#remove column with compound 1 to match to the 2012 GC data that starts at compound C002
OCG_AUC_2021 <- OCG_AUC_2021[,-2]  # Remove the first row #87 of 74 variables

#subset to have just the columns that contain "Peak Area" 1:73
peak_area_cols <- grep("Peak.Area", colnames(OCG_AUC_2021)) 

#the new column names I want to generate will replace the repeating "Peak.Area" names to be "C001" through "C0073" increasing sequentially. 
new_col_names <- paste0("C", sprintf("%03d", seq_along(peak_area_cols)))

#Rename the columns containing "Peak Area" to compound number
colnames(OCG_AUC_2021)[peak_area_cols] <- new_col_names #87 obs and 74 variables

names(OCG_AUC_2021)[names(OCG_AUC_2021) == 'Plant_ID'] <- 'Garden Plant ID' #renaming the ID column to be the same in both datasets: Garden Plant ID

mdITS_OCG_2021$`Garden Plant ID` <- as.character(mdITS_OCG_2021$`Garden Plant ID`) #as character so I can combine them

OCG_AUC_2021 <- merge(mdITS_OCG_2021, OCG_AUC_2021, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 70 obs of 89 variables

#Going to remove everything except chem data and plant ID description
OCG_AUC_2021 <- OCG_AUC_2021[,-c(1:15)] #70 obs of 75 variables

##Save this csv so I can re-read in the data with row 1 being plant ID description
write.csv(OCG_AUC_2021, file = "data_csv/OCG_AUC_2021.csv",row.names = FALSE)

## 2021 cleaned GC data read in####
#read in 2021 chemistry data
OCG_AUC_2021 <- read.csv("data_csv/OCG_AUC_2021.csv", row.names = 1) 
#70 obs of 73 variables
OCG_AUC_2021 <- OCG_AUC_2021[order(row.names(OCG_AUC_2021)),] # order samples alphabetically

OCG_AUC_2021 <- subset(OCG_AUC_2021, row.names(OCG_AUC_2021) %in% row.names(mdITS)) #70 of 73 variables

##Cleaned Full GC data read in ####
OCG_GC <- rbind(OCG_AUC_2012, OCG_AUC_2021) #218 obs of 73 variables
OCG_GC <- OCG_GC[order(row.names(OCG_GC)),] # order samples alphabetically
OCG_GC <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(mdITS)) #218 of 73 variables

## 1uL & 3uL LCMS raw data read in ####
OCG_LCMS_1uL <- read.csv("data_csv/1uL_Injection_Results_LCMS.csv", head=T, check.names = F,stringsAsFactors = T, skip = 1) #120 of 930 variables

sum(duplicated(OCG_LCMS_1uL$Plant_ID)) #one duplicate: plant #273

OCG_LCMS_3uL <- read.csv("data_csv/3uL_injection_results_LCMS.csv", head=T, check.names = F,stringsAsFactors = T, skip = 1) #120 of 929 variables

sum(duplicated(OCG_LCMS_3uL$Plant_ID)) #no duplicates in the dataset

## LCMS 1ul cleaning####
#remove row 8 and 61 since that is the 273 sample
OCG_LCMS_1uL <- OCG_LCMS_1uL[-8, ]  # Remove the row 8 #119 of 930
OCG_LCMS_1uL <- OCG_LCMS_1uL[-61, ]  # Remove the row 61 #118 of 930

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

OCG_LCMS_1uL_2012 <- merge(mdITS_OCG_2012, OCG_LCMS_1uL_2012, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadata of OCG to the samples we have. that only leaves 40 obs of 326 variables

OCG_LCMS_1uL_2012 <- OCG_LCMS_1uL_2012[,-c(1:15,17)] #40 obs of 310 variables

#2021#
#subset the md to only have observations from 2021 to avoid duplicates
OCG_LCMS_1uL_2021 <- subset(OCG_LCMS_1uL, OCG_LCMS_1uL$Year=="2021") #72 observations and 311 variables

mdITS_OCG_2021$`Garden Plant ID` <- as.character(mdITS_OCG_2021$`Garden Plant ID`) #as integer so I can combine them

OCG_LCMS_1uL_2021 <- merge(mdITS_OCG_2021, OCG_LCMS_1uL_2021, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 70 of 326 var

#Going to remove everything except chem data and plant ID
OCG_LCMS_1uL_2021 <- OCG_LCMS_1uL_2021[,-c(1:15,17)] #70 obs of 310 variables

OCG_LCMS_1uL <- data.frame(rbind(OCG_LCMS_1uL_2012,OCG_LCMS_1uL_2021)) #110 of 310 var

write.csv(OCG_LCMS_1uL, file = "data_csv/OCG_LCMS_1uL_cleaned.csv",row.names = FALSE)

## Cleaned 1uL LCMS read in ####
OCG_LCMS_1uL <- read.csv("data_csv/OCG_LCMS_1uL_cleaned.csv", row.names = 1) #110 obs of 309 var
OCG_LCMS_1uL <- OCG_LCMS_1uL[order(row.names(OCG_LCMS_1uL)),] # order samples alphabetically

OCG_LCMS_1uL <- subset(OCG_LCMS_1uL, row.names(OCG_LCMS_1uL) %in% row.names(mdITS)) #109 of 310 variables

## LCMS 3ul cleaning####
#subset to only columns that contain "Peak Area" and "Plant ID". This removes "RT". 
OCG_LCMS_3uL <- OCG_LCMS_3uL [, grepl("Peak.Area|Plant_ID", colnames(OCG_LCMS_3uL))] #120 obs of 617 variables

#subset to remove columns that contain Peak Area Percent and just keep Peak Area.
OCG_LCMS_3uL <- OCG_LCMS_3uL [, !grepl("Peak.Area.Percent", colnames(OCG_LCMS_3uL))] #120 obs of 309 variables

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

OCG_LCMS_3uL_2012 <- merge(mdITS_OCG_2012, OCG_LCMS_3uL_2012, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadata of OCG to the samples we have. that only leaves 41 obs of 332 variables

OCG_LCMS_3uL_2012 <- OCG_LCMS_3uL_2012[,-c(1:15,17)] #41 obs of 309 variables

#2021#
#subset the md to only have observations from 2021 to avoid duplicates
OCG_LCMS_3uL_2021 <- subset(OCG_LCMS_3uL, OCG_LCMS_3uL$Year=="2021") #73 observations and 310 variables

OCG_LCMS_3uL_2021 <- merge(mdITS_OCG_2021, OCG_LCMS_3uL_2021, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 71 obs of 325

#Going to remove everything except chem data and plant ID
OCG_LCMS_3uL_2021 <- OCG_LCMS_3uL_2021[,-c(1:15,17)] #71 obs of 309 variables

OCG_LCMS_3uL <- data.frame(rbind(OCG_LCMS_3uL_2012,OCG_LCMS_3uL_2021)) #112 of 309 var

##Save this csv so I can re-read in the data with row 1 being plant ID
write.csv(OCG_LCMS_3uL, file = "data_csv/OCG_LCMS_3uL_cleaned.csv",row.names = FALSE)

## Cleaned 3uL LCMS data read in ####
OCG_LCMS_3uL <- read.csv("data_csv/OCG_LCMS_3uL_cleaned.csv", row.names = 1) #112 obs of 308 var
OCG_LCMS_3uL <- OCG_LCMS_3uL[order(row.names(OCG_LCMS_3uL)),] # order samples alphabetically

OCG_LCMS_3uL <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(mdITS)) #111 of 308 variables

#Cleaning for PCA####
## 2012 GC data cleaning ####
colSums(is.na(OCG_AUC_2012)) #checking for null values since pca wont run with NAs. All zeroes for each column which indicates there are no NA values

#So now the data needs to be cleaned to only contain compounds that occur in more than 20% of the samples (plants). AKA columns need to be removed that don't have at least 20% of the rows containing a number=NA

#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_AUC_2012))
print(na_proportion)

#Now define the threshold of 20% - there are  compounds that remain after this
threshold <- 0.80

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_AUC_2012_subset <- OCG_AUC_2012[, columns_to_keep] #148 obs of 40 variables

# Replace NA values with zeroes
OCG_AUC_2012_subset[is.na(OCG_AUC_2012_subset)] <- 0

colSums(is.na(OCG_AUC_2012_subset)) #all zeroes for each column which indicates there are no NA values

## 2021 GC data cleaning ####
#So now the data needs to be cleaned to only contain compounds that occur in more than 20% of the samples (plants). AKA columns need to be removed that don't have at least 20% of the rows containing a number=NA
colSums(is.na(OCG_AUC_2012))
#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_AUC_2021))
print(na_proportion)

#Now define the threshold of 20% 
threshold <- 0.8

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_AUC_2021_subset <- OCG_AUC_2021[, columns_to_keep] #70 obs of 36 var variables

# Replace NA values with zeroes
OCG_AUC_2021_subset[is.na(OCG_AUC_2021_subset)] <- 0

colSums(is.na(OCG_AUC_2021_subset)) #all zeroes for each column which indicates there are no NA values

## Full GC data cleaning ####
colSums(is.na(OCG_GC)) #checking for null values since pca wont run with NAs. All zeroes for each column which indicates there are no NA values

#So now the data needs to be cleaned to only contain compounds that occur in more than 20% of the samples (plants). AKA columns need to be removed that don't have at least 20% of the rows containing a number=NA

#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_GC))
print(na_proportion)

#Now define the threshold of 20% - there are  compounds that remain after this
threshold <- 0.80

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_subset <- OCG_GC[, columns_to_keep] #148 obs of 47 variables

# Replace NA values with zeroes
OCG_GC_subset[is.na(OCG_GC_subset)] <- 0

colSums(is.na(OCG_GC_subset)) #all zeroes for each column which indicates there are no NA values

## 1uL LCMS data cleaning ####
colSums(is.na(OCG_LCMS_1uL)) #checking for null values

#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_LCMS_1uL))
print(na_proportion)

#Now define the threshold of 20% - there are  compounds that remain after this
threshold <- 0.80

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_LCMS_1uL_subset <- OCG_LCMS_1uL[, columns_to_keep] #109 of 278 variables

# Replace NA values with zeroes
OCG_LCMS_1uL_subset[is.na(OCG_LCMS_1uL_subset)] <- 0

colSums(is.na(OCG_LCMS_1uL_subset)) #all zeroes for each column which indicates there are no NA values

## 3uL LCMS data cleaning ####
colSums(is.na(OCG_LCMS_3uL)) #checking for null values

#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_LCMS_3uL))
print(na_proportion)

#Now define the threshold of 20% - there are  compounds that remain after this
threshold <- 0.80

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_LCMS_3uL_subset <- OCG_LCMS_3uL[, columns_to_keep] #111 of 283 variables

# Replace NA values with zeroes
OCG_LCMS_3uL_subset[is.na(OCG_LCMS_3uL_subset)] <- 0

colSums(is.na(OCG_LCMS_3uL_subset)) #all zeroes for each column which indicates there are no NA values

#Scaling####
## 2012 GC scaling ####
#Compound
data_normalized_2012 <- scale(OCG_AUC_2012_subset) #this will center the data
#1:148, 1:40

#Plant ID
#transpose the data first to put plant ID as columns and compound as row
OCG_AUC_2012_subset.t <- t(OCG_AUC_2012_subset) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_AUC_2012_subset.t)) #checking for null values

data_normalized_2012.ID <- scale(OCG_AUC_2012_subset.t) 
#1:40, 1:148

## 2021 GC scaling ####
#Compound correlation plot
data_normalized_2021 <- scale(OCG_AUC_2021_subset) #this will center the data
#1:70, 1:36

#Plant correlation plot
#transpose the data first to put plant ID as columns and compound as row 
OCG_AUC_2021_subset.t <- t(OCG_AUC_2021_subset) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_AUC_2021_subset.t)) #checking for null values

data_normalized_2021.ID <- scale(OCG_AUC_2021_subset.t) 
#1:36, 1:70

## Full GC scaling ####
#Compound correlation plot
data_normalized_GC <- scale(OCG_GC_subset) #this will center the data
#1:218, 1:47

#Plant correlation plot
#transpose the data first to put plant ID as columns and compound as row 
OCG_GC_subset.t <- t(OCG_GC_subset) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_GC_subset.t)) #checking for null values

data_normalized_GC.ID <- scale(OCG_GC_subset.t) 
#1:47, 1:218

##1uL LCMS scaling####
data_LCMS1_normalized <- scale(OCG_LCMS_1uL_subset) #this will center the data
#1:109, 1:278

#Plant correlation plot
#transpose the data first to put plant ID as columns and compound as row 
OCG_LCMS_1uL.t <- t(OCG_LCMS_1uL_subset) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_LCMS_1uL.t)) #checking for null values

data_LCMS1_normalized.ID <- scale(OCG_LCMS_1uL.t) 
#1:278, 1:109 

## 3uL LCMS scaling####
data_LCMS3_normalized <- scale(OCG_LCMS_3uL_subset) #this will center the data
#1:111, 1:283

#Plant correlation plot
#transpose the data first to put plant ID as columns and compound as row 
OCG_LCMS_3uL.t <- t(OCG_LCMS_3uL_subset) # transpose rows and columns

#only works with numeric data
colSums(is.na(OCG_LCMS_3uL.t)) #checking for null values

data_LCMS3_normalized.ID <- scale(OCG_LCMS_3uL.t) 
#1:283, 1:111 


# Data visualization####

## Correlation####
### 2012 GC correlation ####
#corr plot for compound
corr_matrix_2012 <- cor(data_normalized_2012) #40 compounds

#corr plot for plant ID
corr_matrix_2012.ID <- cor(data_normalized_2012.ID) #148 plants

### 2021 GC correlation ####
#corr plot for compound
corr_matrix_2021 <- cor(data_normalized_2021) #36 compounds

#corr plot for plant ID
corr_matrix_2021.ID <- cor(data_normalized_2021.ID) # 70 plants

#### GC correlation ####
#corr plot for compound
corr_matrix_GC <- cor(data_normalized_GC) #47 compounds

#corr plot for plant ID
corr_matrix_GC.ID <- cor(data_normalized_GC.ID) # 218 plants

### 1uL LCMS correlation ####
#corr plot for compound
corr_matrix_LCMS1 <- cor(data_LCMS1_normalized)
#large matrix 77284 elements

#corr plot for plant ID
corr_matrix_LCMS1.ID <- cor(data_LCMS1_normalized.ID)
#109

### 3uL LCMS correlation ####
#corr plot for compound
corr_matrix_LCMS3 <- cor(data_LCMS3_normalized)
#large matrix 80089 elements

#corr plot for plant ID
corr_matrix_LCMS3.ID <- cor(data_LCMS3_normalized.ID)
#111 plants

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

fviz_eig(data.pca_2021, addlabels = TRUE) #scree plot

# Graph of the compounds
fviz_pca_var(data.pca_2021, col.var = "black")

#Contribution of each compound
fviz_cos2(data.pca_2021, choice = "var", axes = 1:2)

#Biplot combined with cos2
fviz_pca_var(data.pca_2021, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

## 2012 GC PCA by Plant ID ####
data.pca_2012_ID <- prcomp(data_normalized_2012)
summary(data.pca_2012_ID)
fviz_eig(data.pca_2012_ID, addlabels = TRUE) #scree plot
fviz_cos2(data.pca_2012_ID, choice = "ind", axes = 1:2) #Contribution of each plant
autoplot(data.pca_2012_ID)
autoplot(data.pca_2012_ID, label = TRUE)

plot(data.pca_2012_ID$x[, 1], data.pca_2012_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC comp of 2012 plant by ploidy", 
     col= c("red","blue")[mdITS_OCG_2012$Ploidy],
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
     col= c("pink","brown",'darkgreen')[mdITS_OCG_2012$Subspecies],
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
     xlab="PC 1", ylab="PC 2", 
     main="GC of 2012 plant by subspecies ploidy", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[mdITS_OCG_2012$Subsp_ploidy],
     pch=c(19),
     xlim = range(data.pca_2012_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2012_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")

## 2021 GC PCA by Plant ID####
#prcomp() has improved numerical accuracy, so is preferable to use this function.
data.pca_2021_ID <- prcomp(data_normalized_2021)
summary(data.pca_2021_ID)
fviz_eig(data.pca_2021_ID, addlabels = TRUE) #scree plot
fviz_cos2(data.pca_2021_ID, choice = "ind", axes = 1:2) #Contribution of each plant
autoplot(data.pca_2021_ID)
autoplot(data.pca_2021_ID, label = TRUE)

plot(data.pca_2021_ID$x[, 1], data.pca_2021_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC of 2021 plant by ploidy", 
     col= c("red","blue")[mdITS_OCG_2021$Ploidy],
     pch=c(19),
     xlim = range(data.pca_2021_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2021_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_2021_ID$x[, 1], data.pca_2021_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC of 2021 plant by subspecies", 
     col= c("pink","brown",'darkgreen')[mdITS_OCG_2021$Subspecies],
     pch=c(19),
     xlim = range(data.pca_2021_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2021_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("pink","brown","darkgreen"),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_2021_ID$x[, 1], data.pca_2021_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="GC of 2021 plant by subspecies ploidy", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[mdITS_OCG_2021$Subsp_ploidy],
     pch=c(19),
     xlim = range(data.pca_2021_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2021_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")

## Full GC PCA ####
data.pca_GC_ID <- prcomp(data_normalized_GC)
summary(data.pca_GC_ID)
fviz_eig(data.pca_GC_ID, addlabels = TRUE) #scree plot
fviz_cos2(data.pca_GC_ID, choice = "ind", axes = 1:2) #Contribution of each plant
autoplot(data.pca_GC_ID)
autoplot(data.pca_GC_ID, label = TRUE)

plot(data.pca_GC_ID$x[, 1], data.pca_GC_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="Full GC comp plant by ploidy", 
     col= c("red","blue")[mdITS.OCG$Ploidy],
     pch=c(19),
     xlim = range(data.pca_GC_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_GC_ID$x[, 1], data.pca_GC_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="Full GC comp plant by subspecies", 
     col= c("pink","brown",'darkgreen')[mdITS.OCG$Subspecies],
     pch=c(19),
     xlim = range(data.pca_GC_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("pink","brown","darkgreen"),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_GC_ID$x[, 1], data.pca_GC_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="Full GC comp of plants by subspecies ploidy", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[mdITS.OCG$Subsp_ploidy],
     pch=c(19),
     xlim = range(data.pca_GC_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_GC_ID$x[, 1], data.pca_GC_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="Full GC comp by year", 
     col= rainbow(2)[mdITS.OCG$Year],
     pch=19)
legend("topleft", 
       legend=c("2012","2021"),
       col= c("red","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

##1uL LCMS PCA by plant ID ####
data.pca_LCMS1.ID <- prcomp(data_LCMS1_normalized)
summary(data.pca_LCMS1.ID)
fviz_eig(data.pca_LCMS1.ID, addlabels = TRUE) #13.4%, 8.9%
fviz_cos2(data.pca_LCMS1.ID, choice = "var", axes = 1:2) #Contribution of each compound
fviz_cos2(data.pca_LCMS1.ID, choice = "ind", axes = 1:2) #Contribution of each plant
autoplot(data.pca_LCMS1.ID)
autoplot(data.pca_LCMS1.ID, label = TRUE)

plot(data.pca_LCMS1.ID$x[, 1], data.pca_LCMS1.ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS 1uL plant by ploidy", 
     col= c("red","blue")[mdITS.OCG$Ploidy],
     pch=c(19),
     xlim = range(data.pca_LCMS1.ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS1.ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_LCMS1.ID$x[, 1], data.pca_LCMS1.ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS 1uL plant by subspecies", 
     col= c("pink","brown",'darkgreen')[mdITS.OCG$Subspecies],
     pch=c(19),
     xlim = range(data.pca_LCMS1.ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS1.ID$x[, 2], na.rm = TRUE))
legend("topleft",
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("pink","brown","darkgreen"),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_LCMS1.ID$x[, 1], data.pca_LCMS1.ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS 1uL plant by subspecies ploidy", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[mdITS.OCG$Subsp_ploidy],
     pch=c(19),
     xlim = range(data.pca_LCMS1.ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS1.ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_LCMS1.ID$x[, 1], data.pca_LCMS1.ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS 1uL plant by year", 
     col= c("maroon","cyan")[mdITS.OCG$Year],
     pch=c(19),
     xlim = range(data.pca_LCMS1.ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS1.ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

## 3uL LCMS PCA by Plant ID####
data.pca_LCMS3_ID <- prcomp(data_LCMS3_normalized)
summary(data.pca_LCMS3_ID)
fviz_eig(data.pca_LCMS3_ID, addlabels = TRUE) #14.1% and 8.9%
fviz_cos2(data.pca_LCMS3_ID, choice = "var", axes = 1:2) #Contribution of each compound
fviz_cos2(data.pca_LCMS3_ID, choice = "ind", axes = 1:2) #Contribution of each plant
autoplot(data.pca_LCMS3_ID)
autoplot(data.pca_LCMS3_ID, label = TRUE)

plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS 3uL plant by ploidy", 
     col= c("red","blue")[mdITS.OCG$Ploidy],
     pch=c(19),
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS 3uL plant by subspecies", 
     col= c("pink","brown",'darkgreen')[mdITS.OCG$Subspecies],
     pch=c(19),
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topleft",
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("pink","brown","darkgreen"),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS uL plant by subspecies ploidy", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[mdITS.OCG$Subsp_ploidy],
     pch=c(19),
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")

plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS 3uL plant by year", 
     col= c("maroon","cyan")[mdITS.OCG$Year],
     pch=c(19),
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")


# NMDS plots####
## 2012 GC NMDS plots and stats by plant ID ####
str(mdITS_OCG_2012) #check and make sure levels and factor is correct
# Replace NA with 0
OCG_AUC_2012[is.na(OCG_AUC_2012)] <- 0
#turn abundance data frame into a matrix

OCG_AUC_2012.t <- t(OCG_AUC_2012) #transpose for Plant ID
m_OCG_AUC_2012.t = as.matrix(OCG_AUC_2012.t)

set.seed(7)
#OCG_AUC_2012_ID.nmds <- metaMDS(t(m_OCG_AUC_2012.t), trymax=2000) #solution reached
#save(OCG_AUC_2012_ID.nmds, file = "nmds/OCG_AUC_2012_ID.nmds.rda")
load("nmds/OCG_AUC_2012_ID.nmds.rda")

ordiplot(OCG_AUC_2012_ID.nmds, type = "t",display = "sites",cex = .7)

plot(OCG_AUC_2012_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of 2012 plant by ploidy", 
     col= c("red","blue")[mdITS_OCG_2012$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_AUC_2012_ID.nmds,groups = mdITS_OCG_2012$Ploidy, show.groups = "2n", col = "red")
ordispider(OCG_AUC_2012_ID.nmds,groups = mdITS_OCG_2012$Ploidy, show.groups = "4n", col = "blue")

#### PERMANOVA for 2012 ploidy ##
mdITS.OCG.GC.2012 <- subset(mdITS_OCG_2012, row.names(mdITS_OCG_2012) %in% row.names(OCG_AUC_2012)) #subset md to match AUC samples #147
OCG_AUC_2012 <- subset(OCG_AUC_2012, row.names(OCG_AUC_2012) %in% row.names(mdITS.OCG.GC.2012)) #subset md to match AUC samples #147

OCG_AUC_2012_ploidy <- adonis2(OCG_AUC_2012 ~ mdITS.OCG.GC.2012$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_ploidy #ploidy is significant 0.002

#subspecies nmds
plot(OCG_AUC_2012_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="2012 GC chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.OCG.GC.2012$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_AUC_2012_ID.nmds,groups = mdITS.OCG.GC.2012$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(OCG_AUC_2012_ID.nmds,groups = mdITS.OCG.GC.2012$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(OCG_AUC_2012_ID.nmds,groups = mdITS.OCG.GC.2012$Subspecies, show.groups = "W", col = "goldenrod")

#### PERMANOVA & pairwaise adonis for subspecies ###
OCG_AUC_2012_subsp <- adonis2(OCG_AUC_2012 ~ mdITS.OCG.GC.2012$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_subsp #subspecies is signficant= 0.001

#pairwiseadonis
OCG_AUC_2012_subsp.pw <- pairwise.adonis(OCG_AUC_2012, mdITS.OCG.GC.2012$Subspecies)
OCG_AUC_2012_subsp.pw # sig between all subspecies

# Subspecies ploidy NMDS #
plot(OCG_AUC_2012_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="2012 GC chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS.OCG.GC.2012$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_AUC_2012_ID.nmds,groups = mdITS.OCG.GC.2012$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(OCG_AUC_2012_ID.nmds,groups = mdITS.OCG.GC.2012$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(OCG_AUC_2012_ID.nmds,groups = mdITS.OCG.GC.2012$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(OCG_AUC_2012_ID.nmds,groups = mdITS.OCG.GC.2012$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(OCG_AUC_2012_ID.nmds,groups = mdITS.OCG.GC.2012$Subsp_ploidy, show.groups = "W_4n", col = "purple")

# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_AUC_2012_subspploidy <- adonis2(OCG_AUC_2012 ~ mdITS.OCG.GC.2012$Subsp_ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_subspploidy #subspecies ploidy is significant

#pairwiseadonis
OCG_AUC_2012_subsp_ploidy.pw <- pairwise.adonis(OCG_AUC_2012, mdITS.OCG.GC.2012$Subsp_ploidy)
OCG_AUC_2012_subsp_ploidy.pw 

## 2021 GC NMDS plots and stats by plant ID ####
# Replace NA with 0
OCG_AUC_2021[is.na(OCG_AUC_2021)] <- 0
OCG_AUC_2021 <- OCG_AUC_2021[-which(rownames(OCG_AUC_2021) == "CAT.1.1_2021"), ] #remove this outlier

#turn abundance data frame into a matrix
OCG_AUC_2021.t <- t(OCG_AUC_2021) #transpose for Plant ID
m_OCG_AUC_2021.t = as.matrix(OCG_AUC_2021.t)

set.seed(47)
#OCG_AUC_2021_ID.nmds <- metaMDS(t(m_OCG_AUC_2021.t), trymax=500) #solution reached
#save(OCG_AUC_2021_ID.nmds, file = "nmds/OCG_AUC_2021_ID.nmds.rda")
load("nmds/OCG_AUC_2021_ID.nmds.rda")

ordiplot(OCG_AUC_2021_ID.nmds, type = "t",display = "sites",cex = .7)

mdITS.OCG.GC.2021 <- subset(mdITS_OCG_2021, row.names(mdITS_OCG_2021) %in% row.names(OCG_AUC_2021)) #subset md to match AUC samples #69

plot(OCG_AUC_2021_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of 2021 plant by ploidy", 
     col= c("red","blue")[mdITS.OCG.GC.2021$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_AUC_2021_ID.nmds,groups = mdITS.OCG.GC.2021$Ploidy, show.groups = "2n", col = "red")
ordispider(OCG_AUC_2021_ID.nmds,groups = mdITS.OCG.GC.2021$Ploidy, show.groups = "4n", col = "blue")

#### PERMANOVA for 2021 ploidy ##
OCG_AUC_2021_ploidy <- adonis2(OCG_AUC_2021 ~ mdITS.OCG.GC.2021$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2021_ploidy #ploidy is significant 0.002

#subspecies nmds
plot(OCG_AUC_2021_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="2021 GC chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.OCG.GC.2021$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_AUC_2021_ID.nmds,groups = mdITS.OCG.GC.2021$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(OCG_AUC_2021_ID.nmds,groups = mdITS.OCG.GC.2021$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(OCG_AUC_2021_ID.nmds,groups = mdITS.OCG.GC.2021$Subspecies, show.groups = "W", col = "goldenrod")

#### PERMANOVA & pairwaise adonis for subspecies ###
OCG_AUC_2021_subsp <- adonis2(OCG_AUC_2021 ~ mdITS.OCG.GC.2021$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2021_subsp #subspecies is signficant= 0.001

#pairwiseadonis
OCG_AUC_2021_subsp.pw <- pairwise.adonis(OCG_AUC_2021, mdITS.OCG.GC.2021$Subspecies)
OCG_AUC_2021_subsp.pw # sig between T vs V, T vs W

# Subspecies ploidy NMDS #
plot(OCG_AUC_2021_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="2021 GC chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS.OCG.GC.2021$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_AUC_2021_ID.nmds,groups = mdITS.OCG.GC.2021$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(OCG_AUC_2021_ID.nmds,groups = mdITS.OCG.GC.2021$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(OCG_AUC_2021_ID.nmds,groups = mdITS.OCG.GC.2021$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(OCG_AUC_2021_ID.nmds,groups = mdITS.OCG.GC.2021$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(OCG_AUC_2021_ID.nmds,groups = mdITS.OCG.GC.2021$Subsp_ploidy, show.groups = "W_4n", col = "purple")

# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_AUC_2021_subspploidy <- adonis2(OCG_AUC_2021 ~ mdITS.OCG.GC.2021$Subsp_ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2021_subspploidy #subspecies ploidy is significant

#pairwiseadonis
OCG_AUC_2021_subsp_ploidy.pw <- pairwise.adonis(OCG_AUC_2021, mdITS.OCG.GC.2021$Subsp_ploidy)
OCG_AUC_2021_subsp_ploidy.pw 

## Full GC NMDS plots and stats by plant ID ####
# Replace NA with 0
OCG_GC[is.na(OCG_GC)] <- 0
#OCG_GC <- OCG_GC[-which(rownames(OCG_GC) == "CAT.1.1_2021"), ] #remove this outlier?

#turn abundance data frame into a matrix
OCG_GC.t <- t(OCG_GC) #transpose for Plant ID
m_OCG_GC.t = as.matrix(OCG_GC.t)

set.seed(4)
OCG_GC_ID.nmds <- metaMDS(t(m_OCG_GC.t), trymax=500) #solution reached
save(OCG_GC_ID.nmds, file = "nmds/OCG_GC_ID.nmds.rda")
load("nmds/OCG_GC_ID.nmds.rda")

ordiplot(OCG_GC_ID.nmds, type = "t",display = "sites",cex = .7)

mdITS.OCG.GC <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(OCG_GC)) #subset md to match AUC samples #217
OCG_GC <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(mdITS.OCG))

plot(OCG_GC_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of plant by ploidy", 
     col= c("red","blue")[mdITS.OCG.GC$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_GC_ID.nmds,groups = mdITS.OCG.GC$Ploidy, show.groups = "2n", col = "red")
ordispider(OCG_GC_ID.nmds,groups = mdITS.OCG.GC$Ploidy, show.groups = "4n", col = "blue")

#### PERMANOVA for 2021 ploidy ##
OCG_GC_ploidy <- adonis2(OCG_GC ~ mdITS.OCG.GC$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_GC_ploidy #ploidy is significant 0.003

#subspecies nmds
plot(OCG_GC_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="GC chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.OCG.GC$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_GC_ID.nmds,groups = mdITS.OCG.GC$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(OCG_GC_ID.nmds,groups = mdITS.OCG.GC$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(OCG_GC_ID.nmds,groups = mdITS.OCG.GC$Subspecies, show.groups = "W", col = "goldenrod")

#### PERMANOVA & pairwaise adonis for subspecies ###
OCG_GC_subsp <- adonis2(OCG_GC ~ mdITS.OCG.GC$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_GC_subsp #subspecies is signficant= 0.001

#pairwiseadonis
OCG_GC_subsp.pw <- pairwise.adonis(OCG_GC, mdITS.OCG.GC$Subspecies)
OCG_GC_subsp.pw # sig between T vs V, V vs W

# Subspecies ploidy NMDS #
plot(OCG_GC_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS.OCG.GC$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_GC_ID.nmds,groups = mdITS.OCG.GC$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(OCG_GC_ID.nmds,groups = mdITS.OCG.GC$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(OCG_GC_ID.nmds,groups = mdITS.OCG.GC$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(OCG_GC_ID.nmds,groups = mdITS.OCG.GC$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(OCG_GC_ID.nmds,groups = mdITS.OCG.GC$Subsp_ploidy, show.groups = "W_4n", col = "purple")

# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_GC_subspploidy <- adonis2(OCG_GC ~ mdITS.OCG.GC$Subsp_ploidy,by="margin") # Bray-Curtis is the default metric
OCG_GC_subspploidy #subspecies ploidy is significant

#pairwiseadonis
OCG_GC_subsp_ploidy.pw <- pairwise.adonis(OCG_GC, mdITS.OCG.GC$Subsp_ploidy)
OCG_GC_subsp_ploidy.pw 

#by year
plot(OCG_GC_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of plant by year", 
     col= c("maroon","cyan")[mdITS.OCG.GC$Year],
     pch=c(19))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

OCG_GC_yr <- adonis2(OCG_GC ~ mdITS.OCG.GC$Year,by="margin") # Bray-Curtis is the default metric
OCG_GC_yr #yearis signficant= 0.001

##1uL LCMS NMDS by plant ID ####
# Replace NA with 0
OCG_LCMS_1uL[is.na(OCG_LCMS_1uL)] <- 0
OCG_LCMS_1uL.t <- t(OCG_LCMS_1uL)
m_OCG_LCMS_1uL.t = as.matrix(OCG_LCMS_1uL.t)

set.seed(53)
#OCG_LCMS_1uL_ID.nmds <- metaMDS(t(m_OCG_LCMS_1uL.t), trymax=500) #solution reached
#save(OCG_LCMS_1uL_ID.nmds, file = "nmds/OCG_LCMS_1uL_ID.nmds.rda")
load("nmds/OCG_LCMS_1uL_ID.nmds.rda")

ordiplot(OCG_LCMS_1uL_ID.nmds, type = "t",display = "sites",cex = .7)

mdITS.OCG.LCMS.1 <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(OCG_LCMS_1uL)) #subset md to match AUC samples #109

plot(OCG_LCMS_1uL_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 1uL plant by ploidy", 
     col= c("red","blue")[mdITS.OCG.LCMS.1$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Ploidy, show.groups = "2n", col = "red")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Ploidy, show.groups = "4n", col = "blue")

#### PERMANOVA for ploidy ##
OCG_LCMS1_2021_ploidy <- adonis2(OCG_LCMS_1uL ~ mdITS.OCG.LCMS.1$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_LCMS1_2021_ploidy #ploidy is significant 0.001

#subspecies nmds
plot(OCG_LCMS_1uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 1uL chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.OCG.LCMS.1$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Subspecies, show.groups = "W", col = "goldenrod")

#### PERMANOVA & pairwaise adonis for subspecies ###
OCG_LCMS1_subsp <- adonis2(OCG_LCMS_1uL ~ mdITS.OCG.LCMS.1$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_LCMS1_subsp #subspecies is signficant= 0.001

#pairwiseadonis
OCG_LCMS1_subsp.pw <- pairwise.adonis(OCG_LCMS_1uL, mdITS.OCG.LCMS.1$Subspecies)
OCG_LCMS1_subsp.pw # sig between all subspecies

# Subspecies ploidy NMDS #
plot(OCG_LCMS_1uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 1uL chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS.OCG.LCMS.1$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Subsp_ploidy, show.groups = "W_4n", col = "purple")


# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_LCMS1_subspploidy <- adonis2(OCG_LCMS_1uL ~ mdITS.OCG.LCMS.1$Subsp_ploidy,by="margin") # Bray-Curtis is the default metric
OCG_LCMS1_subspploidy #subspecies ploidy is significant

OCG_LCMS1_subsp_ploidy.pw <- pairwise.adonis(OCG_LCMS_1uL, mdITS.OCG.LCMS.1$Subsp_ploidy)
OCG_LCMS1_subsp_ploidy.pw 

plot(OCG_LCMS_1uL_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 1uL plant by year", 
     col= c("maroon","cyan")[mdITS.OCG.LCMS.1$Year],
     pch=c(19))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Year, show.groups = "2012", col = "maroon")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = mdITS.OCG.LCMS.1$Year, show.groups = "2021", col = "cyan")

# PERMANOVAS for year ##
OCG_LCMS1_year <- adonis2(OCG_LCMS_1uL ~ mdITS.OCG.LCMS.1$Year,by="margin") # Bray-Curtis is the default metric
OCG_LCMS1_year #year is significant

## 3uL LCMS NMDS by plant ID ####
OCG_LCMS_3uL[is.na(OCG_LCMS_3uL)] <- 0
OCG_LCMS_3uL.t <- t(OCG_LCMS_3uL)
m_OCG_LCMS_3uL.t = as.matrix(OCG_LCMS_3uL.t)

set.seed(38)
#OCG_LCMS_3uL_ID.nmds <- metaMDS(t(m_OCG_LCMS_3uL.t), trymax=500) #solution reached
#save(OCG_LCMS_3uL_ID.nmds, file = "nmds/OCG_LCMS_3uL_ID.nmds.rda")
load("nmds/OCG_LCMS_3uL_ID.nmds.rda")

ordiplot(OCG_LCMS_3uL_ID.nmds, type = "t",display = "sites",cex = .7)

mdITS.OCG.LCMS.3 <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(OCG_LCMS_3uL)) #subset md to match AUC samples #111

plot(OCG_LCMS_3uL_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 3uL plant by ploidy", 
     col= c("red","blue")[mdITS.OCG.LCMS.3$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Ploidy, show.groups = "2n", col = "red")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Ploidy, show.groups = "4n", col = "blue")

#### PERMANOVA for ploidy ##
OCG_LCMS3_2021_ploidy <- adonis2(OCG_LCMS_3uL ~ mdITS.OCG.LCMS.3$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_LCMS3_2021_ploidy #ploidy is significant 0.001

#subspecies nmds
plot(OCG_LCMS_3uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 3uL chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[mdITS.OCG.LCMS.3$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Subspecies, show.groups = "W", col = "goldenrod")

#### PERMANOVA & pairwaise adonis for subspecies ###
OCG_LCMS3_subsp <- adonis2(OCG_LCMS_3uL ~ mdITS.OCG.LCMS.3$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_LCMS3_subsp #subspecies is signficant= 0.001

#pairwiseadonis
OCG_LCMS3_subsp.pw <- pairwise.adonis(OCG_LCMS_3uL, mdITS.OCG.LCMS.3$Subspecies)
OCG_LCMS3_subsp.pw # sig between all subspecies

# Subspecies ploidy NMDS #
plot(OCG_LCMS_3uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 3uL chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[mdITS.OCG.LCMS.3$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Subsp_ploidy, show.groups = "W_4n", col = "purple")

# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_LCMS3_subspploidy <- adonis2(OCG_LCMS_3uL ~ mdITS.OCG.LCMS.3$Subsp_ploidy,by="margin") # Bray-Curtis is the default metric
OCG_LCMS3_subspploidy #subspecies ploidy is significant

OCG_LCMS3_subsp_ploidy.pw <- pairwise.adonis(OCG_LCMS_3uL, mdITS.OCG.LCMS.3$Subsp_ploidy)
OCG_LCMS3_subsp_ploidy.pw 

plot(OCG_LCMS_3uL_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 3uL plant by year", 
     col= c("maroon","cyan")[mdITS.OCG.LCMS.3$Year],
     pch=c(19))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Year, show.groups = "2012", col = "maroon")
ordispider(OCG_LCMS_3uL_ID.nmds,groups = mdITS.OCG.LCMS.3$Year, show.groups = "2021", col = "cyan")

# PERMANOVAS for year ##
OCG_LCMS3_year <- adonis2(OCG_LCMS_3uL ~ mdITS.OCG.LCMS.3$Year,by="margin") # Bray-Curtis is the default metric
OCG_LCMS3_year #year is significant


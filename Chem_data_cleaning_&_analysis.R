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
if (!require("ggfortify")) {install.packages("ggfortify"); require("ggfortify")}

# Read data in#### 
##Metadata read in
md <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #505 obs of 16 variables.
md <- md[order(row.names(md)),] #alphabetical

md.OCG <- subset(md, md$Project=="OCG") #subsetting to just the Orchard common garden plants 246 of 16 variables

### Remove duplicates
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1', 'ORT.1.5_2012')
md.OCG <- md.OCG[!rownames(md.OCG) %in% rows_to_remove, ]

## Remove negative control
md.OCG <- md.OCG[!(row.names(md.OCG) == "NEG_8-28-21"),]
md.OCG <- md.OCG[!(row.names(md.OCG) == "NEG_10-2-20"),]

## Remove MTW.3.7.R_2012
md.OCG <- md.OCG[!(row.names(md.OCG) == "MTW.3.7.R_2012"),]#236 of 16 var

#make variables factor to plot
md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant")] <- lapply(md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant")], as.factor)

md.OCG.2012 <- subset(md.OCG, md.OCG$Year=="2012") #159
md.OCG.2012$`Garden Plant ID` <- as.factor(md.OCG.2012$`Garden Plant ID`) #as factor so I can combine them
md.OCG.2012 <- droplevels(md.OCG.2012)
str(md.OCG.2012)

md.ITS.OCG.2012

md.OCG.2021 <- subset(md.OCG, md.OCG$Year=="2021") #76
md.OCG.2021 <- droplevels(md.OCG.2021)
str(md.OCG.2021)

## 2012 GC raw data read in#
OCG_AUC_2012 <- read.csv("data_csv/OCG_2012_GC.csv", head=T,check.names = F,stringsAsFactors = T, skip = 1) #183 obs of 221 var

## 2012 GC Cleaning ####
# Remove the second column since it is empty
OCG_AUC_2012 <- OCG_AUC_2012[, -2]
names(OCG_AUC_2012)[1] <- "Plant_ID"

sum(duplicated(OCG_AUC_2012$Plant_ID)) #4 duplicates (control, cocktails, and blanks)

# Remove rows where "empty", "Empty", "Ct, "CT", "M", "w"  occur
OCG_AUC_2012 <- OCG_AUC_2012[!(apply(OCG_AUC_2012, 1, function(row) any(grepl("^empty|^Empty|^Ct|^CT |^M|w", row)))), ] #167 of 220 variables

sum(duplicated(OCG_AUC_2012$Plant_ID)) # 0 duplicates

#subset to only columns that contain "Peak Area" and "Plant ID". This removes "RT". 
OCG_AUC_2012 <- OCG_AUC_2012 [, grepl("Peak.Area|Plant_ID", colnames(OCG_AUC_2012))]

#subset to remove columns that contain Peak Area Percent and just keep Peak Area.
OCG_AUC_2012 <- OCG_AUC_2012 [, !grepl("Peak.Area.Percent", colnames(OCG_AUC_2012))] 

#subset to have just the columns that contain "Peak Area" 1:73 #shifted by one compared to 2021 since there was no compound 1
peak_area_cols <- grep("Peak.Area", colnames(OCG_AUC_2012)) 

#the new column names will replace the repeating "Peak.Area" names to be "C001" through "C0073" increasing sequentially. 
new_col_names <- paste0("C", sprintf("%03d", seq_along(peak_area_cols)))

#Rename the columns containing "Peak Area" to compound number
colnames(OCG_AUC_2012)[peak_area_cols] <- new_col_names 

names(OCG_AUC_2012)[names(OCG_AUC_2012) == 'Plant_ID'] <- 'Garden Plant ID' 

OCG_AUC_2012 <- merge(md.OCG, OCG_AUC_2012, by="Garden Plant ID") #157 obs of 89

OCG_AUC_2012 <- OCG_AUC_2012[,-c(1:15)] #removing everything except area under the curve 157 obs of 74 variables

##Save this csv
write.csv(OCG_AUC_2012, file = "data_csv/OCG_AUC_2012_cleaned.csv",row.names = FALSE)

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

md.ITS.OCG.2021$`Garden Plant ID` <- as.character(md.ITS.OCG.2021$`Garden Plant ID`) #as character so I can combine them

OCG_AUC_2021 <- merge(md.ITS.OCG.2021, OCG_AUC_2021, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 70 obs of 89 variables

#Going to remove everything except chem data and plant ID description
OCG_AUC_2021 <- OCG_AUC_2021[,-c(1:15)] #70 obs of 74 variables

##Save this csv so I can re-read in the data with row 1 being plant ID description
write.csv(OCG_AUC_2021, file = "data_csv/OCG_AUC_2021.csv",row.names = FALSE)

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

md.ITS.OCG.2012$`Garden Plant ID` <- as.character(md.ITS.OCG.2012$`Garden Plant ID`) #as character so I can combine them

OCG_LCMS_1uL_2012 <- merge(md.ITS.OCG.2012, OCG_LCMS_1uL_2012, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadata of OCG to the samples we have. that only leaves 40 obs of 326 variables

OCG_LCMS_1uL_2012 <- OCG_LCMS_1uL_2012[,-c(1:15,17)] #40 obs of 310 variables

#2021#
#subset the md to only have observations from 2021 to avoid duplicates
OCG_LCMS_1uL_2021 <- subset(OCG_LCMS_1uL, OCG_LCMS_1uL$Year=="2021") #72 observations and 311 variables

md.ITS.OCG.2021$`Garden Plant ID` <- as.character(md.ITS.OCG.2021$`Garden Plant ID`) #as integer so I can combine them

OCG_LCMS_1uL_2021 <- merge(md.ITS.OCG.2021, OCG_LCMS_1uL_2021, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 70 of 326 var

#Going to remove everything except chem data and plant ID
OCG_LCMS_1uL_2021 <- OCG_LCMS_1uL_2021[,-c(1:15,17)] #70 obs of 310 variables

OCG_LCMS_1uL <- data.frame(rbind(OCG_LCMS_1uL_2012,OCG_LCMS_1uL_2021)) #110 of 310 var

write.csv(OCG_LCMS_1uL, file = "data_csv/OCG_LCMS_1uL_cleaned.csv",row.names = FALSE)

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

OCG_LCMS_3uL_2012 <- merge(md.ITS.OCG.2012, OCG_LCMS_3uL_2012, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadata of OCG to the samples we have. that only leaves 41 obs of 332 variables

OCG_LCMS_3uL_2012 <- OCG_LCMS_3uL_2012[,-c(1:15,17)] #41 obs of 309 variables

#2021#
#subset the md to only have observations from 2021 to avoid duplicates
OCG_LCMS_3uL_2021 <- subset(OCG_LCMS_3uL, OCG_LCMS_3uL$Year=="2021") #73 observations and 310 variables

OCG_LCMS_3uL_2021 <- merge(md.ITS.OCG.2021, OCG_LCMS_3uL_2021, by="Garden Plant ID") # Joining the dataframes so I can match/subset metadta of OCG to the samples we have. 71 obs of 325

#Going to remove everything except chem data and plant ID
OCG_LCMS_3uL_2021 <- OCG_LCMS_3uL_2021[,-c(1:15,17)] #71 obs of 309 variables

OCG_LCMS_3uL <- data.frame(rbind(OCG_LCMS_3uL_2012,OCG_LCMS_3uL_2021)) #112 of 309 var

##Save this csv so I can re-read in the data with row 1 being plant ID
write.csv(OCG_LCMS_3uL, file = "data_csv/OCG_LCMS_3uL_cleaned.csv",row.names = FALSE)

#Clear Global environment#############
#Cleaned data read in ####
#Metadata read in 
md <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #set to correct file path 505 obs of 16 variables.
md <- md[order(row.names(md)),] # order samples alphabetically

### Subsetting to just the plant in the common garden (OCG)#
md.OCG <- subset(md, md$Project=="OCG") #subsetting to just the Orchard common garden plants #246 of 16 variables

### Remove duplicates from asv#
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1', 'ORT.1.5_2012')
md.OCG <- md.OCG[!rownames(md.OCG) %in% rows_to_remove, ]

## Remove negative control
md.OCG <- md.OCG[!(row.names(md.OCG) == "NEG_8-28-21"),]
md.OCG <- md.OCG[!(row.names(md.OCG) == "NEG_10-2-20"),]

## Remove MTW.3.7.R_2012
md.OCG <- md.OCG[!(row.names(md.OCG) == "MTW.3.7.R_2012"),] #we arent sure what the R represents. 236 of 16 var

#make variables factor to plot
md.OCG$Ploidy <- as.factor(md.OCG$Ploidy)
md.OCG$Subspecies <- as.factor(md.OCG$Subspecies)
md.OCG$Subsp_ploidy <- as.factor(md.OCG$Subsp_ploidy)
md.OCG$Year <- as.factor(md.OCG$Year)
md.OCG$Plant <- as.factor(md.OCG$Plant)

#subset the md to only have observations from 2012 to avoid duplicates
md.ITS.OCG.2012 <- subset(md.OCG, md.OCG$Year=="2012") 
#159
md.ITS.OCG.2012$`Garden Plant ID` <- as.factor(md.ITS.OCG.2012$`Garden Plant ID`) #as factor so I can combine them
md.ITS.OCG.2012 <- droplevels(md.ITS.OCG.2012)
str(md.ITS.OCG.2012)

#subset the md to only have observations from 2021 to avoid duplicates
md.ITS.OCG.2021 <- subset(md.OCG, md.OCG$Year=="2021") 
#76
md.ITS.OCG.2021 <- droplevels(md.ITS.OCG.2021)
str(md.ITS.OCG.2021)


## 2012 cleaned GC data read in
#read in cleaned 2012 GC data with row 1 being plant ID
OCG_AUC_2012 <- read.csv("data_csv/OCG_AUC_2012_cleaned.csv", row.names = 1) 
#157 obs of 73 variables
OCG_AUC_2012 <- OCG_AUC_2012[order(row.names(OCG_AUC_2012)),] # order samples alphabetically

OCG_AUC_2012 <- subset(OCG_AUC_2012, row.names(OCG_AUC_2012) %in% row.names(md)) #148 of 73 variables

## 2021 cleaned GC data read in
#read in 2021 chemistry data
OCG_AUC_2021 <- read.csv("data_csv/OCG_AUC_2021.csv", row.names = 1) 
#70 obs of 73 variables
OCG_AUC_2021 <- OCG_AUC_2021[order(row.names(OCG_AUC_2021)),] # order samples alphabetically

OCG_AUC_2021 <- subset(OCG_AUC_2021, row.names(OCG_AUC_2021) %in% row.names(md)) #70 of 73 variables

##Cleaned Full GC data read in 
OCG_GC <- rbind(OCG_AUC_2012, OCG_AUC_2021) #218 obs of 73 variables
OCG_GC <- OCG_GC[order(row.names(OCG_GC)),] # order samples alphabetically
OCG_GC <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(md)) #218 of 73 variables
md.OCG.GC <- subset(md.OCG.GC, row.names(md.OCG.GC) %in% row.names(OCG_GC)) #218 of 73 variables

## Cleaned 1uL LCMS read in 
OCG_LCMS_1uL <- read.csv("data_csv/OCG_LCMS_1uL_cleaned.csv", row.names = 1) #110 obs of 309 var
OCG_LCMS_1uL <- OCG_LCMS_1uL[order(row.names(OCG_LCMS_1uL)),] # order samples alphabetically

OCG_LCMS_1uL <- subset(OCG_LCMS_1uL, row.names(OCG_LCMS_1uL) %in% row.names(md)) #109 of 310 variables

## Cleaned 3uL LCMS data read in 
OCG_LCMS_3uL <- read.csv("data_csv/OCG_LCMS_3uL_cleaned.csv", row.names = 1) #112 obs of 308 var
OCG_LCMS_3uL <- OCG_LCMS_3uL[order(row.names(OCG_LCMS_3uL)),] # order samples alphabetically

OCG_LCMS_3uL <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md)) #111 of 308 variables


#Alpha diversity ####
OCG.GC.shannon <- diversity(OCG_GC)
OCG.GC.ef <- exp(OCG.GC.shannon)
OCG.GC.ef.r <- round(OCG.GC.ef)

md.OCG.GC <- cbind(md.OCG.GC, effective_species = OCG.GC.ef.r)

glm.OCG.GC <- glm(effective_species ~ Subspecies + Year, family = poisson, data = md.OCG.GC)
summary(glm.OCG.GC)

plot(allEffects(glm.OCG.GC))

plot(md.OCG.GC$Year,OCG.GC.ef)
plot(md.OCG.GC$Subspecies,OCG.GC.ef)

ggplot(md.OCG.GC, aes(Year, effective_species))+
  geom_boxplot(aes(group = Year, fill = Year))+
  theme_classic()

ggplot(data = md.OCG.GC, mapping = aes(x = Subspecies, y = effective_species, fill = Subspecies)) +
  geom_boxplot() +
  theme_classic()

OCG.LCMS.shannon <- diversity(OCG_LCMS_3uL)
OCG.LCMS.ef <- exp(OCG.LCMS.shannon)
OCG.LCMS.ef.r <- round(OCG.LCMS.ef)

md.OCG.LCMS.3 <- cbind(md.OCG.LCMS.3, effective_species = OCG.LCMS.ef.r)
glm.OCG.LCMS <- glm(effective_species ~ Subspecies + Year, family = poisson, data = md.OCG.LCMS.3)
summary(glm.OCG.LCMS)

plot(allEffects(glm.OCG.LCMS))

plot(md.OCG.LCMS.3$Year,OCG.LCMS.ef)
plot(md.OCG.LCMS.3$Subspecies,OCG.LCMS.ef)

ggplot(md.OCG.LCMS.3, aes(Year, effective_species))+
  geom_boxplot(aes(group = Year, fill = Year))+
  theme_classic()

ggplot(data = md.OCG.LCMS.3, mapping = aes(x = Subspecies, y = effective_species, fill = Subspecies)) +
  geom_boxplot() +
  theme_classic()


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
     col= c("red","blue")[md.ITS.OCG.2012$Ploidy],
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
     col= c("pink","brown",'darkgreen')[md.ITS.OCG.2012$Subspecies],
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
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.ITS.OCG.2012$Subsp_ploidy],
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
     col= c("red","blue")[md.ITS.OCG.2021$Ploidy],
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
     col= c("pink","brown",'darkgreen')[md.ITS.OCG.2021$Subspecies],
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
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.ITS.OCG.2021$Subsp_ploidy],
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
     col= c("red","blue")[md.OCG$Ploidy],
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
     col= c("pink","brown",'darkgreen')[md.OCG$Subspecies],
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
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG$Subsp_ploidy],
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
     col= rainbow(2)[md.OCG$Year],
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
     col= c("red","blue")[md.OCG$Ploidy],
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
     col= c("pink","brown",'darkgreen')[md.OCG$Subspecies],
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
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG$Subsp_ploidy],
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
     col= c("maroon","cyan")[md.OCG$Year],
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
     col= c("red","blue")[md.OCG$Ploidy],
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
     col= c("pink","brown",'darkgreen')[md.OCG$Subspecies],
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
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG$Subsp_ploidy],
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
     col= c("maroon","cyan")[md.OCG$Year],
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
str(md.ITS.OCG.2012) #check and make sure levels and factor is correct
# Replace NA with 0
OCG_AUC_2012[is.na(OCG_AUC_2012)] <- 0
#turn abundance data frame into a matrix

#scaling
OCG_AUC_2012.t <- t(OCG_AUC_2012) #transpose for Plant ID
#OCG_AUC_2012.t <- rarefy(OCG_AUC_2012.t)
m_OCG_AUC_2012.t = as.matrix(OCG_AUC_2012.t)

set.seed(7)
#OCG_AUC_2012_ID.nmds <- metaMDS(t(m_OCG_AUC_2012.t), trymax=1000) #solution reached
#save(OCG_AUC_2012_ID.nmds, file = "nmds/OCG_AUC_2012_ID.nmds.rda")
load("nmds/OCG_AUC_2012_ID.nmds.rda")

md.OCG.GC.2012 <- subset(md.ITS.OCG.2012, row.names(md.ITS.OCG.2012) %in% row.names(OCG_AUC_2012)) #subset md to match AUC samples #147
OCG_AUC_2012 <- subset(OCG_AUC_2012, row.names(OCG_AUC_2012) %in% row.names(md.OCG.GC.2012)) #subset md to match AUC samples #147

ordiplot(OCG_AUC_2012_ID.nmds, type = "t",display = "sites",cex = .7)

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

#### PERMANOVA for 2012 ploidy ##
OCG_AUC_2012_ploidy <- adonis2(OCG_AUC_2012 ~ md.OCG.GC.2012$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_ploidy #ploidy is significant 0.002

#subspecies nmds
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

#### PERMANOVA & pairwaise adonis for subspecies ###
OCG_AUC_2012_subsp <- adonis2(OCG_AUC_2012 ~ md.OCG.GC.2012$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_subsp #subspecies is signficant= 0.001

#pairwiseadonis
OCG_AUC_2012_subsp.pw <- pairwise.adonis(OCG_AUC_2012, md.OCG.GC.2012$Subspecies)
OCG_AUC_2012_subsp.pw # sig between all subspecies

# Subspecies ploidy NMDS #
plot(OCG_AUC_2012_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="2012 GC chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[md.OCG.GC.2012$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")

# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_AUC_2012_subspploidy <- adonis2(OCG_AUC_2012 ~ md.OCG.GC.2012$Subsp_ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2012_subspploidy #subspecies ploidy is significant

#pairwiseadonis
OCG_AUC_2012_subsp_ploidy.pw <- pairwise.adonis(OCG_AUC_2012, md.OCG.GC.2012$Subsp_ploidy)
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

md.OCG.GC.2021 <- subset(md.ITS.OCG.2021, row.names(md.ITS.OCG.2021) %in% row.names(OCG_AUC_2021)) #subset md to match AUC samples #69

plot(OCG_AUC_2021_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of 2021 plant by ploidy", 
     col= c("red","blue")[md.OCG.GC.2021$Ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_AUC_2021_ID.nmds,groups = md.OCG.GC.2021$Ploidy, show.groups = "2n", col = "red")
ordispider(OCG_AUC_2021_ID.nmds,groups = md.OCG.GC.2021$Ploidy, show.groups = "4n", col = "blue")

#### PERMANOVA for 2021 ploidy ##
OCG_AUC_2021_ploidy <- adonis2(OCG_AUC_2021 ~ md.OCG.GC.2021$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2021_ploidy #ploidy is significant 0.002

#subspecies nmds
plot(OCG_AUC_2021_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="2021 GC chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.GC.2021$Subspecies],
     pch=c(19))
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_AUC_2021_ID.nmds,groups = md.OCG.GC.2021$Subspecies, show.groups = "T", col = "olivedrab")
ordispider(OCG_AUC_2021_ID.nmds,groups = md.OCG.GC.2021$Subspecies, show.groups = "V", col = "cadetblue")
ordispider(OCG_AUC_2021_ID.nmds,groups = md.OCG.GC.2021$Subspecies, show.groups = "W", col = "goldenrod")

#### PERMANOVA & pairwaise adonis for subspecies ###
OCG_AUC_2021_subsp <- adonis2(OCG_AUC_2021 ~ md.OCG.GC.2021$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2021_subsp #subspecies is signficant= 0.001

#pairwiseadonis
OCG_AUC_2021_subsp.pw <- pairwise.adonis(OCG_AUC_2021, md.OCG.GC.2021$Subspecies)
OCG_AUC_2021_subsp.pw # sig between T vs V, T vs W

# Subspecies ploidy NMDS #
plot(OCG_AUC_2021_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="2021 GC chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[md.OCG.GC.2021$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_AUC_2021_ID.nmds,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(OCG_AUC_2021_ID.nmds,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(OCG_AUC_2021_ID.nmds,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(OCG_AUC_2021_ID.nmds,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(OCG_AUC_2021_ID.nmds,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "W_4n", col = "purple")

# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_AUC_2021_subspploidy <- adonis2(OCG_AUC_2021 ~ md.OCG.GC.2021$Subsp_ploidy,by="margin") # Bray-Curtis is the default metric
OCG_AUC_2021_subspploidy #subspecies ploidy is significant

#pairwiseadonis
OCG_AUC_2021_subsp_ploidy.pw <- pairwise.adonis(OCG_AUC_2021, md.OCG.GC.2021$Subsp_ploidy)
OCG_AUC_2021_subsp_ploidy.pw 

## Full GC NMDS plots and stats by plant ID ####
# Replace NA with 0
OCG_GC[is.na(OCG_GC)] <- 0

#turn abundance data frame into a matrix
OCG_GC.t <- t(OCG_GC) #transpose for Plant ID
m_OCG_GC.t = as.matrix(OCG_GC.t)

set.seed(4)
#OCG_GC_ID.nmds <- metaMDS(t(m_OCG_GC.t), trymax=500) #solution reached
#save(OCG_GC_ID.nmds, file = "nmds/OCG_GC_ID.nmds.rda")
load("nmds/OCG_GC_ID.nmds.rda")

ordiplot(OCG_GC_ID.nmds, type = "t",display = "sites",cex = .7)

md.OCG.GC <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC)) #subset md to match AUC samples #217
OCG_GC <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(md.OCG))

plot(OCG_GC_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of plant by ploidy", 
     col= c("red","blue")[md.OCG.GC$Ploidy],
     pch=c(17,19)[md.OCG.GC$Year])
text(OCG_GC_ID.nmds$points[,1:2], 
     labels=md.OCG.GC$Description, 
     pos=1, 
     cex=0.4)
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n") 

#### PERMANOVA for 2021 ploidy ##
OCG_GC_ploidy <- adonis2(OCG_GC ~ md.OCG.GC$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_GC_ploidy #ploidy is significant 0.003

OCG_GC_ploidy_yr <- adonis2(OCG_GC ~ md.OCG.GC$Ploidy*md.OCG.GC$Year,by="margin") # Bray-Curtis is the default metric
OCG_GC_ploidy_yr #sig

#subspecies nmds
plot(OCG_GC_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="GC chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.GC$Subspecies],
     pch=c(17,19)[md.OCG.GC$Year])
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n") 

#### PERMANOVA & pairwaise adonis for subspecies ###
OCG_GC_subsp <- adonis2(OCG_GC ~ md.OCG.GC$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_GC_subsp #subspecies is signficant= 0.001

OCG_GC_subsp_yr <- adonis2(OCG_GC ~ md.OCG.GC$Subspecies*md.OCG.GC$Year,by="margin") # Bray-Curtis is the default metric
OCG_GC_subsp_yr #sig 

#pairwiseadonis
OCG_GC_subsp.pw <- pairwise.adonis(OCG_GC, md.OCG.GC$Subspecies)
OCG_GC_subsp.pw # sig between T vs V, V vs W

# Subspecies ploidy NMDS #
plot(OCG_GC_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[md.OCG.GC$Subsp_ploidy],
     pch=c(17,19)[md.OCG.GC$Year])
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n") 

# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_GC_subspploidy <- adonis2(OCG_GC ~ md.OCG.GC$Subsp_ploidy,by="margin") # Bray-Curtis is the default metric
OCG_GC_subspploidy #subspecies ploidy is significant

OCG_GC_subspploidy_yr <- adonis2(OCG_GC ~ md.OCG.GC$Subsp_ploidy+md.OCG.GC$Year,by="margin") # Bray-Curtis is the default metric
OCG_GC_subspploidy_yr #sig

#pairwiseadonis
OCG_GC_subsp_ploidy.pw <- pairwise.adonis(OCG_GC, md.OCG.GC$Subsp_ploidy)
OCG_GC_subsp_ploidy.pw 

#by year
plot(OCG_GC_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of plant by year", 
     col= c("maroon","cyan")[md.OCG.GC$Year],
     pch=c(19))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

OCG_GC_yr <- adonis2(OCG_GC ~ md.OCG.GC$Year,by="margin") # Bray-Curtis is the default metric
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

md.OCG.LCMS.1 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_LCMS_1uL)) #subset md to match AUC samples #109

plotcolor <- c("olivedrab","cadetblue","magenta","blue","pink","green","darkgreen","khaki","goldenrod","yellow","cornflowerblue")

plot(OCG_LCMS_1uL_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 1uL plant by ploidy", 
     col = plotcolor [md.OCG.LCMS.1$Location],
     pch=c(17,19)[md.OCG.LCMS.1$Ploidy])
legend("topleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.6,
       bty = "n")
legend("topright", 
       legend=c("2n","4n"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

#### PERMANOVA for ploidy ##
OCG_LCMS1_ploidy <- adonis2(OCG_LCMS_1uL ~ md.OCG.LCMS.1$Ploidy,by="margin") # Bray-Curtis is the default metric
OCG_LCMS1_ploidy #ploidy is significant 0.001

OCG_LCMS1_loc <- adonis2(OCG_LCMS_1uL ~ md.OCG.LCMS.1$Location,by="margin") # Bray-Curtis is the default metric
OCG_LCMS1_loc #location is significant

#subspecies nmds
plot(OCG_LCMS_1uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 1uL chemistry by subspecies",
     col = plotcolor [md.OCG.LCMS.1$Location],
     pch=c(17,19,15)[md.OCG.LCMS.1$Subspecies])
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       pch=c(17,19,15),
       cex=0.8,
       bty = "n")

#### PERMANOVA & pairwaise adonis for subspecies ###
OCG_LCMS1_subsp <- adonis2(OCG_LCMS_1uL ~ md.OCG.LCMS.1$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_LCMS1_subsp #subspecies is signficant= 0.001

#pairwiseadonis
OCG_LCMS1_subsp.pw <- pairwise.adonis(OCG_LCMS_1uL, md.OCG.LCMS.1$Subspecies)
OCG_LCMS1_subsp.pw # sig between all subspecies

# Subspecies ploidy NMDS #
plot(OCG_LCMS_1uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 1uL chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[md.OCG.LCMS.1$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = md.OCG.LCMS.1$Subsp_ploidy, show.groups = "T_2n", col = "red")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = md.OCG.LCMS.1$Subsp_ploidy, show.groups = "T_4n", col = "orange")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = md.OCG.LCMS.1$Subsp_ploidy, show.groups = "V_2n", col = "green")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = md.OCG.LCMS.1$Subsp_ploidy, show.groups = "V_4n", col = "cyan")
ordispider(OCG_LCMS_1uL_ID.nmds,groups = md.OCG.LCMS.1$Subsp_ploidy, show.groups = "W_4n", col = "purple")


# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_LCMS1_subspploidy <- adonis2(OCG_LCMS_1uL ~ md.OCG.LCMS.1$Subsp_ploidy ,by="margin") # Bray-Curtis is the default metric
OCG_LCMS1_subspploidy #subspecies ploidy is significant

OCG_LCMS1_subsp_ploidy.pw <- pairwise.adonis(OCG_LCMS_1uL, md.OCG.LCMS.1$Subsp_ploidy)
OCG_LCMS1_subsp_ploidy.pw 

plot(OCG_LCMS_1uL_ID.nmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 1uL plant by year", 
     col= c("maroon","cyan")[md.OCG.LCMS.1$Year],
     pch=c(19))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

# PERMANOVAS for year ##
OCG_LCMS1_year <- adonis2(OCG_LCMS_1uL ~ md.OCG.LCMS.1$Year,by="margin") # Bray-Curtis is the default metric
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

md.OCG.LCMS.3 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_LCMS_3uL)) #subset md to match AUC samples #111

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

#### PERMANOVA for ploidy ##
OCG_LCMS3_ploidy <- adonis2(OCG_LCMS_3uL ~ md.OCG.LCMS.3$Ploidy ,by="margin") # Bray-Curtis is the default metric
OCG_LCMS3_ploidy #ploidy is significant 0.001

#subspecies nmds
plot(OCG_LCMS_3uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 3uL chemistry by subspecies",
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

#### PERMANOVA & pairwaise adonis for subspecies ###
OCG_LCMS3_subsp <- adonis2(OCG_LCMS_3uL ~ md.OCG.LCMS.3$Subspecies,by="margin") # Bray-Curtis is the default metric
OCG_LCMS3_subsp #subspecies is signficant= 0.001

#pairwiseadonis
OCG_LCMS3_subsp.pw <- pairwise.adonis(OCG_LCMS_3uL, md.OCG.LCMS.3$Subspecies)
OCG_LCMS3_subsp.pw # sig between all subspecies

# Subspecies ploidy NMDS #
plot(OCG_LCMS_3uL_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="LCMS 3uL chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[md.OCG.LCMS.3$Subsp_ploidy],
     pch=c(19))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")

# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_LCMS3_subspploidy <- adonis2(OCG_LCMS_3uL ~ md.OCG.LCMS.3$Subsp_ploidy + md.OCG.LCMS.3$Location ,by="margin") # Bray-Curtis is the default metric
OCG_LCMS3_subspploidy #subspecies ploidy is significant

OCG_LCMS3_subsp_ploidy.pw <- pairwise.adonis(OCG_LCMS_3uL, md.OCG.LCMS.3$Subsp_ploidy)
OCG_LCMS3_subsp_ploidy.pw 

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

# PERMANOVAS for year ##
OCG_LCMS3_year <- adonis2(OCG_LCMS_3uL ~ md.OCG.LCMS.3$Year ,by="margin") # Bray-Curtis is the default metric
OCG_LCMS3_year #year is significant

###2012 LCMS data ####
str(md.ITS.OCG.2012) #check and make sure levels and factor is correct
#subset to just 2012 data
OCG_LCMS_2012 <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.ITS.OCG.2012)) #subset md to match AUC samples #40
md.OCG.LCMS.2012 <- subset(md.ITS.OCG.2012, row.names(md.ITS.OCG.2012) %in% row.names(OCG_LCMS_2012)) #subset md to match AUC samples #40

# Replace NA with 0
OCG_LCMS_2012[is.na(OCG_LCMS_2012)] <- 0
#turn abundance data frame into a matrix

OCG_LCMS_2012.t <- t(OCG_LCMS_2012) #transpose for Plant ID
m_OCG_LCMS_2012.t = as.matrix(OCG_LCMS_2012.t)

set.seed(450)
OCG_LCMS_2012_ID.nmds <- metaMDS(t(m_OCG_LCMS_2012.t), trymax=500) #solution reached
#save(OCG_LCMS_2012_ID.nmds, file = "nmds/OCG_LCMS_2012_ID.nmds.rda")
load("nmds/OCG_LCMS_2012_ID.nmds.rda")

ordiplot(OCG_LCMS_2012_ID.nmds, type = "t",display = "sites",cex = .7)

#by ploidy
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
# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_LCMS_2012_ploidy <- adonis2(OCG_LCMS_2012 ~ md.OCG.LCMS.2012$Ploidy, by="margin") # Bray-Curtis is the default metric
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
OCG_LCMS_2012_subspploidy <- adonis2(OCG_LCMS_2012 ~ md.OCG.LCMS.2012$Subsp_ploidy, by="margin") # Bray-Curtis is the default metric
OCG_LCMS_2012_subspploidy #subspecies ploidy is significant

#by location 
plotcolor <- c("olivedrab","cadetblue","magenta","blue","orange","green","darkgreen","firebrick","lightgoldenrod","mediumaquamarine","cornflowerblue")
plot(OCG_LCMS_2012_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 2012 3uL chemistry by location",
     col = plotcolor [md.OCG.LCMS.2012$Location],
     pch=17)
legend("bottomleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.6,
       bty = "n")
# PERMANOVAS and pairwise adonis for location ##
OCG_LCMS_2012_loc <- adonis2(OCG_LCMS_2012 ~ md.OCG.LCMS.2012$Location, by="margin") # Bray-Curtis is the default metric
OCG_LCMS_2012_loc #location is significant

###2021 LCMS data ####
str(md.ITS.OCG.2021) #check and make sure levels and factor is correct
#subset to just 2012 data
OCG_LCMS_2021 <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.ITS.OCG.2021)) #subset md to match AUC samples #71
md.OCG.LCMS.2021 <- subset(md.ITS.OCG.2021, row.names(md.ITS.OCG.2021) %in% row.names(OCG_LCMS_2021)) #subset md to match AUC samples #71

# Replace NA with 0
OCG_LCMS_2021[is.na(OCG_LCMS_2021)] <- 0
#turn abundance data frame into a matrix

OCG_LCMS_2021.t <- t(OCG_LCMS_2021) #transpose for Plant ID
m_OCG_LCMS_2021.t = as.matrix(OCG_LCMS_2021.t)

set.seed(21)
OCG_LCMS_2021_ID.nmds <- metaMDS(t(m_OCG_LCMS_2021.t), trymax=500) #solution reached
#save(OCG_LCMS_2021_ID.nmds, file = "nmds/OCG_LCMS_2021_ID.nmds.rda")
load("nmds/OCG_LCMS_2021_ID.nmds.rda")

ordiplot(OCG_LCMS_2021_ID.nmds, type = "t",display = "sites",cex = .7)

#by ploidy
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
# PERMANOVAS and pairwise adonis for ploidy ##
OCG_LCMS_2021_ploidy <- adonis2(OCG_LCMS_2021 ~ md.OCG.LCMS.2021$Ploidy, by="margin") # Bray-Curtis is the default metric
OCG_LCMS_2021_ploidy #ploidy is significant

#by subspecies- only tridentata samples in the 2012 LCMS 3uL data
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
# PERMANOVAS and pairwise adonis for subspecies ##
OCG_LCMS_2021_subsp <- adonis2(OCG_LCMS_2021 ~ md.OCG.LCMS.2021$Subspecies, by="margin") # Bray-Curtis is the default metric
OCG_LCMS_2021_subsp #subspecies is significant

# Subspecies ploidy NMDS there is T_2n and T_4n 
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
# PERMANOVAS and pairwise adonis for subspecies ploidy ##
OCG_LCMS_2021_subspploidy <- adonis2(OCG_LCMS_2021 ~ md.OCG.LCMS.2021$Subsp_ploidy, by="margin") # Bray-Curtis is the default metric
OCG_LCMS_2021_subspploidy #subspecies ploidy is significant

#by location 
plotcolor <- c("olivedrab","cadetblue","magenta","blue","orange","green","darkgreen","firebrick","lightgoldenrod","mediumaquamarine","cornflowerblue")
plot(OCG_LCMS_2021_ID.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 2021 3uL chemistry by location",
     col = plotcolor [md.OCG.LCMS.2021$Location],
     pch=19)
legend("bottomleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.8,
       bty = "n")
# PERMANOVAS and pairwise adonis for location ##
OCG_LCMS_2021_loc <- adonis2(OCG_LCMS_2021 ~ md.OCG.LCMS.2021$Location, by="margin") # Bray-Curtis is the default metric
OCG_LCMS_2021_loc #location is significant

# Binary jicard analysis ####
#LCMS data
set.seed(8)
OCG_LCMS_3uL_ID.jdis <- vegdist(t(m_OCG_LCMS_3uL.t), method = "jaccard", binary = TRUE)
OCG_LCMS_3uL_ID.jnmds <- metaMDS(OCG_LCMS_3uL_ID.jdis, trymax=500) #solution reached

ordiplot(OCG_LCMS_3uL_ID.jnmds, type = "t",display = "sites",cex = .7)

plotcolor <- c("olivedrab","cadetblue","magenta","blue","orange","green","darkgreen","firebrick","lightgoldenrod","mediumaquamarine","cornflowerblue")

#by subspecies
plot(OCG_LCMS_3uL_ID.jnmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="LCMS 3uL chemistry by subspecies and location",
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

#ploidy 
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

#by subspecies ploidy 
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

#GC data
set.seed(6)
# OCG_GC <- OCG_GC[-which(colnames(OCG_GC) == "CAT.1.1_2021"), ] #remove this outlier
OCG_GC_ID.jdis <- vegdist(t(m_OCG_GC.t), method = "jaccard", binary = TRUE)
#OCG_GC_ID.jnmds <- metaMDS(OCG_GC_ID.jdis, trymax=500) #solution reached
#save(OCG_GC_ID.jnmds, file = "jnmds/OCG_GC_ID.jnmds.rda")
load("jnmds/OCG_GC_ID.jnmds.rda")

ordiplot(OCG_GC_ID.jnmds, type = "t",display = "sites",cex = .7)

#by ploidy
plot(OCG_GC_ID.jnmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of plant by ploidy", 
     col= c("red","blue")[md.OCG.GC$Ploidy],
     pch=c(17,19)[md.OCG.GC$Year])
legend("topleft", 
       legend=c("2n","4n"),
       col= c("red","blue"),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n") 

#by subspecies
plot(OCG_GC_ID.jnmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     main="GC chemistry by subspecies",
     col= c("olivedrab","cadetblue","goldenrod")[md.OCG.GC$Subspecies],
     pch=c(17,19)[md.OCG.GC$Year])
legend("topleft", 
       legend=c("Tridentata","Vaseyana","Wyomingensis"),
       col= c("olivedrab","cadetblue","goldenrod"),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n") 

#by subspecies ploidy
plot(OCG_GC_ID.jnmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC chemistry by subspecies and ploidy", 
     col= c("red","orange","green","cyan","purple")[md.OCG.GC$Subsp_ploidy],
     pch=c(17,19)[md.OCG.GC$Year])
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("red","orange","green","cyan","purple"),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

#by year
plot(OCG_GC_ID.jnmds$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     main="GC of plant by year", 
     col= c("maroon","cyan")[md.OCG.GC$Year],
     pch=c(19))
legend("topleft", 
       legend=c("2012","2021"),
       col= c("maroon","cyan"),
       pch=19,
       cex=0.8,
       bty = "n")

# PCoA ####
#GC
GC_dist_matrix <- vegdist(OCG_GC, method = "bray")
GC_pcoa <- cmdscale(GC_dist_matrix) #classic multidimensional scaling (cmdscale)

plot(GC_pcoa[,1], GC_pcoa[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "PCoA Visualization", 
     col= c("red","blue")[md.OCG.GC$Ploidy],
     pch = 16)

plot(GC_pcoa[,1], GC_pcoa[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "GC PCoA Visualization by year", 
     col= c("red","cyan")[md.OCG.GC$Year],
     pch = 16)
legend("topleft", 
       legend=c("2012","2021"),
       col= c("red","cyan"),
       pch=16,
       cex=0.8,
       bty = "n")

#LCMS
LCMS_dist_matrix <- vegdist(OCG_LCMS_3uL, method = "bray")
LCMS_pcoa <- cmdscale(LCMS_dist_matrix)

plot(LCMS_pcoa[,1], LCMS_pcoa[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "PCoA Visualization", 
     col= c("red","blue")[md.OCG.LCMS.3$Ploidy],
     pch = 16)

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

plot(LCMS_pcoa[,1], LCMS_pcoa[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "LCMS PCoA Visualization by subspecies", 
     col= c("pink","brown",'darkgreen')[md.OCG.LCMS.3$Subspecies],
     pch = 16)
legend("topleft", 
       legend=c("T","V","W"),
       col= c("pink","brown",'darkgreen'),
       pch=16,
       cex=0.8,
       bty = "n")

plot(LCMS_pcoa[,1], LCMS_pcoa[,2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS PCoA Visualization by subspecies ploidy", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.LCMS.3$Subsp_ploidy],
     pch=c(16))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=16,
       cex=0.8,
       bty = "n")

plot(LCMS_pcoa[,1], LCMS_pcoa[,2],
     xlab="PC 1", ylab="PC 2", 
     main="LCMS PCoA Visualization by Location", 
     col= plotcolor[md.OCG.LCMS.3$Location],
     pch=c(16))
legend("topleft", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=16,
       cex=0.8,
       bty = "n")

#Stable isotope data analysis#### 
#2012 data read in
stable_iso_data_2012<-read.csv("data_csv/Stable_isotope_2012.csv")
str(stable_iso_data_2012)

# Extracting numbers from Plant_ID column
stable_iso_data_2012$Sample.ID <- as.numeric(gsub("[^0-9]", "", stable_iso_data_2012$ID))
stable_iso_data_2012$Sample.ID <- as.character(stable_iso_data_2012$Sample.ID)
str(stable_iso_data_2012)

#2021 data read in
stable_iso_data_2021<-read.csv("data_csv/Stable_isotope_2021.csv")
str(stable_iso_data_2021)

#renaming the columns to combine the data together
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "Ampl..28"] <- "Ampl28"
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "Ampl..44"] <- "Ampl44"
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "d15N"] <- "Delta15N"
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "d13C"] <- "Delta13C"

#add these two dataframes together
stable_iso_df <- full_join(stable_iso_data_2012,stable_iso_data_2021)

#structure the data to subspecies being T W and V... wont need to do this. I need tocheck what subspecies 134677 is 
stable_iso_data_2021$Subspecies<-recode_factor(stable_iso_data_2021$Subspecies, Tridentata="T", Vaseyana="V", Wyomingensis="W")
levels(stable_iso_data_2021$Subspecies)
levels(stable_iso_data_2021$Subspecies)<- list(Tridentata="T", Vaseyana="V", Wyomingensis="W")
levels(stable_iso_data_2021$Subspecies)


## Remove V-134677 (not sure what subspecies since it says "T/W"). look into later
stable_iso_data_2021 <- stable_iso_data_2021[!(row.names(stable_iso_data_2021) == "31"),]
print(stable_iso_data_2021)

dis<-vegdist(stable_iso_data_2012$Delta15N) #Bray-curtis distance between samples: quantifies differences in the overall taxonomic composition between two samples

groups<-factor(stable_iso_data_2012$Subspecies, labels = c("T","V","W")) #using the md for year 

## Calculate multivariate dispersions
mod<-betadisper(dis,groups)
mod

#Perform test
anova(mod)

#Permutation test for F
permutest(mod, permutations = 99,pairwise = TRUE) 

dis<-vegdist(stable_iso_data_2012$Delta13C) #Bray-curtis distance between samples: quantifies differences in the overall taxonomic composition between two samples

groups<-factor(stable_iso_data_2012$Subspecies, labels = c("T","V","W")) #using the md for year 

## Calculate multivariate dispersions
mod<-betadisper(dis,groups)
mod

#Perform test
anova(mod)

TukeyHSD(p1)

pw.comparison.15N<-aov(Delta15N~Subspecies,data=stable_iso_data_2012)
TukeyHSD(pw.comparison.15N)

pw.comparison.13C<-aov(Delta13C~Subspecies,data=stable_iso_data_2012)
TukeyHSD(pw.comparison.13C)

#Permutation test for F
permutest(mod, permutations = 99,pairwise = TRUE) 

fit<-lm(Delta15N~Subspecies,data = stable_iso_data_2012)
summary(fit) #T vs V= 0.006, T vs W= 0.024, and W vs V= 0.164.

#pairwiseadonis
stable_iso.subsp.pw <- pairwise.adonis(stable_iso_data_2012$Delta15N, stable_iso_data_2012$Subspecies)
# stable_iso.subsp.pw #T vs V= < 2e-16, T vs W= 5.25e-05, and W vs V= 0.0391.

library(pairwiseAdonis)

pairwise.adonis<-pairwise.adonis2(stable_iso_data_2012$Delta15N ~ stable_iso_data_2012$Subspecies)
pairwise.adonis

fit1<-lm(Delta13C~Subspecies,data = stable_iso_data_2012)
summary(fit1)#0.039

fit2<-lm(d15N~Subspecies,data = stable_iso_data_2021)

summary(fit2) #vaseyana significnatly different
fit3<-lm(d13C~Subspecies,data = stable_iso_data_2021)
summary(fit3) #p-value=0.037

#differences between year
stable_isotope_full_data<-read.csv("stable_isotope_full_data")

#model
fit4<-lm(Delta15N~Year,data=)
summary(fit4)

library(ggplot2)
p1<-ggplot(stable_iso_data_2012,aes(Subspecies,Delta15N))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+ggtitle("2012")+ theme(legend.position = "none")

p2<-ggplot(stable_iso_data_2021,aes(Subspecies,d15N))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+
  ylab("Delta15N")+ggtitle("2021")+theme(legend.position = "none")

library(gridExtra)
grid.arrange(p1, p2, nrow = 1)

ggplot(stable_iso_data_2012,aes(Delta13C,Delta15N))+
  geom_point(aes(color=Subspecies))+theme_classic()+
  scale_color_manual(values=c("olivedrab","cadetblue","goldenrod"))+
  ggtitle("2012")+theme(legend.position = "none")+geom_text(aes(label=ID))

ggplot(stable_iso_data_2021,aes(d13C,d15N))+
  geom_point(aes(color=Subspecies))+theme_classic()+
  scale_color_manual(values=c("olivedrab","cadetblue","goldenrod"))+ ggtitle("2021")+theme(legend.position = "none")+
  geom_text(aes(label=Lab))

stable_iso_data_2021$Subspecies

fit<-anova_test(stable_iso_data_2012$Delta13C~stable_iso_data_2012Subspecies,
                detailed = TRUE)
fit<-anova(stable_iso_data_2012$Delta13C,stable_iso_data_2012$Subspecies)

library(ggpubr)
ggscatterhist(
  stable_iso_data_2012, x = "Delta13C", y = "Delta15N", group="Subspecies",
  color = "Subspecies", fill= "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.plot = "boxplot",
  margin.params = list(fill = "Subspecies", color = c("olivedrab", "cadetblue", "goldenrod"), size = 0.2),
  ggtheme = theme_bw()
)+
  geom_text(aes(label=))

ggscatterhist(
  stable_iso_data_2021, x = "d13C", y = "d15N", group="Subspecies",
  color = "Subspecies", fill= "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.plot = "boxplot",
  margin.params = list(fill = "Subspecies", color = c("olivedrab", "cadetblue", "goldenrod"), size = 0.2),
  ggtheme = theme_bw()
)+
  geom_text(aes(label=Sample.ID))

#geom_text(aes(label= Sample.ID))

#mean point, x and y with sd. geom_errorbarh

# Grouped Scatter plot with marginal density plots
ggscatterhist(
  stable_iso_data_2012, x = "Delta13C", y = "Delta15N", group = "Subspecies",
  color = "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.params = list(fill = "Subspecies", color = "black", size = 0.2)
)

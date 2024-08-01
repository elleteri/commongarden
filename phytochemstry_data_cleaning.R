#For Chemistry data cleaning #
setwd("/Users/ellehorwath/Documents/Orchard_Common_Garden/commongarden")

#Packages
cran_packages <- c("dplyr", "readr", "tidyr", "tidyverse")

# Function to install and load CRAN packages
install_and_load_cran <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Install and load all CRAN packages
sapply(cran_packages, install_and_load_cran)

# Read raw data in
##METADATA READ IN
md <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #505 obs of 21 variables.
md <- md[order(row.names(md)),] #alphabetical

### Subsetting to just the plant in the common garden (OCG)#
md.OCG <- subset(md, md$Project=="OCG")#246 of 21 variables

### Remove duplicates from md
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1', 'COVW.2.4_2012')
md.OCG <- md.OCG[!rownames(md.OCG) %in% rows_to_remove, ] #199

md.OCG.2012 <- subset(md.OCG, md.OCG$Year=="2012")
sum(duplicated(md.OCG.2012$`Garden Plant ID`)) #7
md.OCG.2021 <- subset(md.OCG, md.OCG$Year=="2021")
sum(duplicated(md.OCG.2021$`Garden Plant ID`))

write.csv(md.OCG, file = "data_csv/metadata_OCG.csv")

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

OCG_GC_2012 <- merge(md.OCG.2012, OCG_GC_2012, by="Garden Plant ID") #165 obs of 95

OCG_GC_2012 <- OCG_GC_2012[,-c(1:15,17:21)] #removing everything except area under the curve 165 obs of 75 variables

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

#subset to "Peak Area" 1:75
peak_area_cols <- grep("Peak.Area", colnames(OCG_GC_2021)) 

#the new column "C001" through "C0075" increasing sequentially. 
new_col_names <- paste0("C", sprintf("%03d", seq_along(peak_area_cols)))

#Rename
colnames(OCG_GC_2021)[peak_area_cols] <- new_col_names #87 obs and 75 variables

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

OCG_GC_2021 <- merge(md.OCG.2021, OCG_GC_2021, by="Garden Plant ID") #70 obs of 95 variables

#Going to remove everything except AUC and plant ID description and RT
#OCG_GC_2021_w_RT <- OCG_GC_2021[,-c(1:15,17:18)] #70 obs of 149 variables

#Going to remove everything except AUC and plant ID description
OCG_GC_2021 <- OCG_GC_2021[,-c(1:15,17:21)] #70 obs of 75 variables

# #Save csv with RT
# write.csv(OCG_GC_2021_w_RT, file = "data_csv/OCG_GC_2021_cleaned_w_RT.csv",row.names = FALSE)

#Save csv
write.csv(OCG_GC_2021, file = "data_csv/OCG_GC_2021.csv",row.names = FALSE)

#FULL GC COMBINED
OCG_GC <- rbind(OCG_GC_2012, OCG_GC_2021) #235 obs of 75 variables

#FULL GC COMBINED W RT
#OCG_GC_w_RT <- rbind(OCG_GC_2012_w_RT, OCG_GC_2021_w_RT) #226 obs of 149 variables

#restructure the RT GC full data. 
#remove plant ID?
#make the row names the plant id
#rownames(OCG_GC_w_RT) <- OCG_GC_w_RT[,1]
#pivot longer
# df_long <- OCG_GC_w_RT %>%
#   pivot_longer(cols = everything(), 
#                names_to = c(".value", "Compound"), 
#                names_pattern = "(C|RT)0*(\\d+)")
# 
# # Rename columns for clarity
# names(df_long) <- c("Compound", "RT", "PeakArea")
# 
# # ##Save csv
# # write.csv(OCG_GC_w_RT, file = "data_csv/OCG_GC_w_RT_full_clean.csv")
# 
# OCG_GC_w_RT_2012 <- read.csv("data_csv/OCG_GC_2012_cleaned_w_RT.csv", head=T, check.names = F,stringsAsFactors = T, row.names = 1)
# OCG_GC_w_RT_2021 <- read.csv("data_csv/OCG_GC_2021_cleaned_w_RT.csv", head=T, check.names = F,stringsAsFactors = T, row.names = 1)
# OCG_GC_w_RT <- rbind(OCG_GC_w_RT_2012, OCG_GC_w_RT_2021) #226 obs of 149 variables
# OCG_GC_w_RT <- read.csv("data_csv/OCG_GC_w_RT_full_clean.csv", head=T, check.names = F,stringsAsFactors = T)

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

OCG_LCMS_3uL_2012 <- OCG_LCMS_3uL_2012[,-c(1:15,17:22)] #42 obs of 309 variables

#2021 LCMS
OCG_LCMS_3uL_2021 <- subset(OCG_LCMS_3uL, OCG_LCMS_3uL$Year=="2021") #73 observations and 310 variables

OCG_LCMS_3uL_2021 <- merge(md.OCG.2021, OCG_LCMS_3uL_2021, by="Garden Plant ID") #71 obs of 325

OCG_LCMS_3uL_2021 <- OCG_LCMS_3uL_2021[,-c(1:15,17:22)] #71 obs of 309 variables

#FULL LCMS
OCG_LCMS_3uL <- data.frame(rbind(OCG_LCMS_3uL_2012,OCG_LCMS_3uL_2021)) #113 of 309 var

##Save csv
write.csv(OCG_LCMS_3uL, file = "data_csv/OCG_LCMS_3uL_cleaned.csv",row.names = FALSE)

#Clear Global Environment
rm(list = ls())

# Orchard Common Garden Chemistry 2021 LCMS data cleaning
# Set working directory and load necessary packages for reading in and tidying the data.####
setwd("~/Documents/Orchard_Common _Garden/commongarden")
if (!require("readr")) {install.packages("readr"); require("readr")}
if (!require("dplyr")) {install.packages("dplyr"); require("dplyr")}
if (!require("tidyr")) {install.packages("tidyr"); require("tidyr")}
if (!require("xcms"))  {install.packages("xcms"); require("xcms")}


#read in metadata####
mdITS_OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T)

OCG_LCMS_1uL <- read.csv("data_csv/1uL_Injection_Results_LCMS.csv", head=T, check.names = F,stringsAsFactors = T)

sum(duplicated(OCG_LCMS_1uL)) #no duplicates in the dataset

OCG_LCMS_3uL <- read.csv("data_csv/3uL_injection_results_LCMS.csv", head=T, check.names = F,stringsAsFactors = T)

sum(duplicated(OCG_LCMS_3uL)) #no duplicates in the dataset




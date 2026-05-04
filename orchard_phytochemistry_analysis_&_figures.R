# Orchard, ID Common Garden data analysis and visualization for sagebrush phytochemistry 
# Install and load necessary packages ####
library("ANCOMBC")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("vegan")
library("pairwiseAdonis")
library("readr")
library("effects")
library("viridis")
library("fields")
library("BiocManager")
library("SummarizedExperiment")
library("pheatmap")
library("sf")
library("lme4")

# Cleaned data read in####
#METADATA
md_OCG <- read.csv("data_csv/metadata_OCG.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #246 obs of 22 variables.
md_OCG <- md_OCG[order(row.names(md_OCG)),]

### Remove duplicates, negatve controls, double samples, and MTW.3.7.R_2012
md_OCG <- md_OCG[!rownames(md_OCG) %in% c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1', 'ORT.1.5_2012', 'NEG_8-28-21', 'NEG_10-2-20', 'MTW.3.7.R_2012'), ] #234 of 21 var

rownames(md_OCG) <- gsub("v[12]$", "", rownames(md_OCG), ignore.case = TRUE) # remove the v1 or v2 at the end of the row names

# make variables factor to plot and droplevels
md_OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Ecoregion","Description", "Plant Group")] <- lapply(md_OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Ecoregion","Description", "Plant Group")], as.factor)
md_OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Ecoregion","Description","Plant Group")] <- lapply(md_OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Ecoregion","Description", "Plant Group")], droplevels)
str(md_OCG)

#subset the md to only have observations from 2012 to avoid duplicates
md_OCG_2012 <- subset(md_OCG, md_OCG$Year=="2012") #158
str(md_OCG_2012)

#subset the md to only have observations from 2021 to avoid duplicates
md_OCG_2021 <- subset(md_OCG, md_OCG$Year=="2021") #76
str(md_OCG_2021)

##FULL CLEAN GC
OCG_GC <- read.csv("data_csv/OCG_GC_full_clean.csv", row.names = 1)#227 obs of 74 variables
OCG_GC <- OCG_GC[order(row.names(OCG_GC)),]
OCG_GC <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(md_OCG)) #217 of 74 variables
md_OCG_GC <- subset(md_OCG, row.names(md_OCG) %in% row.names(OCG_GC)) #217 of 21 variables
OCG_GC[is.na(OCG_GC)] <- 0

## 2012 CLEAN GC
OCG_GC_2012 <- read.csv("data_csv/OCG_GC_2012_cleaned.csv", row.names = 1, check.names = FALSE) #157 obs of 74 variables
OCG_GC_2012 <- OCG_GC_2012[order(row.names(OCG_GC_2012)),]
OCG_GC_2012 <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md_OCG)) #147 of 74 variables
md_OCG_GC_2012 <- subset(md_OCG, row.names(md_OCG) %in% row.names(OCG_GC_2012)) #147
OCG_GC_2012[is.na(OCG_GC_2012)] <- 0

## 2021 CLEAN GC
OCG_GC_2021 <- read.csv("data_csv/OCG_GC_2021.csv", row.names = 1)#70 obs of 74 variables
OCG_GC_2021 <- OCG_GC_2021[order(row.names(OCG_GC_2021)),] 
OCG_GC_2021 <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md_OCG)) #70 of 74 variables
md_OCG_GC_2021 <- subset(md_OCG, row.names(md_OCG) %in% row.names(OCG_GC_2021)) #70
OCG_GC_2021[is.na(OCG_GC_2021)] <- 0

## CLEAN LCMS 
OCG_LCMS_3uL <- read.csv("data_csv/OCG_LCMS_3uL_cleaned.csv", row.names = 1) #112
OCG_LCMS_3uL <- OCG_LCMS_3uL[order(row.names(OCG_LCMS_3uL)),]
OCG_LCMS_3uL <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md_OCG)) #111
md_OCG_LCMS <- subset(md_OCG, row.names(md_OCG) %in% row.names(OCG_LCMS_3uL)) #111 
OCG_LCMS_3uL[is.na(OCG_LCMS_3uL)] <- 0

###########################################
## GC alpha diversity (GLMM) ####
OCG_GC_shannon <- diversity(OCG_GC)
OCG_GC_ef <- exp(OCG_GC_shannon)
OCG_GC_ef_r <- round(OCG_GC_ef)

md_OCG_GC <- cbind(md_OCG_GC, GC_compounds = OCG_GC_ef_r)

glmm_GC <- glmer(
  formula = GC_compounds ~ Year + Ecoregion + Subsp_ploidy + (1 | Plant),
  data = md_OCG_GC,
  family = Gamma(link = "log"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)  # increase from default (1e4)
  )
)
summary(glmm_GC)
plot(allEffects(glmm_GC))

# Supplementary Figure 1D: Boxplots of compound richness for GC across subspecies and years
ggplot(data = md_OCG_GC, aes(Subsp_ploidy, GC_compounds, fill =Subsp_ploidy)) +
  geom_boxplot() +
  scale_fill_manual(values = c("pink","brown","darkgreen",'tan','lightblue'))+
  labs(y = "Number of Compounds", x = "Subspecies + Ploidy")+
  theme_classic()+
  theme(legend.position = "none")

# Supplementary Figure 1B: Boxplots of compound richness for GC across years
ggplot(md_OCG_GC, aes(Year, GC_compounds))+
  geom_boxplot(aes(group = Year, fill = Year))+
  scale_fill_manual(values  = c("maroon", "cyan"))+
  labs(y = "Number of Compounds")+
  theme_classic()+
  theme(legend.position = "none")

# Tukey HSD pairwise comparison in compound richness across subspecies
gc_compounds_aov <- aov(md_OCG_GC$GC_compounds ~ md_OCG_GC$Subsp_ploidy)
TukeyHSD(gc_compounds_aov)

lcms_compounds_aov <- aov(md_OCG_LCMS$Compounds ~ md_OCG_LCMS$Subsp_ploidy)
TukeyHSD(lcms_compounds_aov)

## LCMS alpha diversity (GLMM) ####
OCG_LCMS_shannon <- diversity(OCG_LCMS_3uL)
OCG_LCMS_ef <- exp(OCG_LCMS_shannon)
OCG_LCMS_ef_r <- round(OCG_LCMS_ef)

md_OCG_LCMS <- cbind(md_OCG_LCMS, Compounds = OCG_LCMS_ef_r)

glmm_LCMS <- glmer(
  formula = Compounds ~ Year + Ecoregion + Subsp_ploidy + (1 | Plant),
  data = md_OCG_LCMS,
  family = Gamma(link = "log"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)  # increase from default (1e4)
  )
)
summary(glmm_LCMS)

# Supplementary Figure 1A: Boxplots of compound richness for LCMS across years
ggplot(md_OCG_LCMS, aes(Year, Compounds))+
  geom_boxplot(aes(group = Year, fill = Year))+
  scale_fill_manual(values = c("maroon","cyan"))+
  labs(y = "Number of Compounds")+
  theme_classic()+
  theme(legend.position = "none")

# Supplementary Figure 1C: Boxplots of compound richness for LCMS across subspecies and years
ggplot(data = md_OCG_LCMS, aes(Subsp_ploidy, Compounds, fill = Subsp_ploidy)) +
  geom_boxplot() +
  scale_fill_manual(values = c("pink","brown","darkgreen",'tan','lightblue'))+
  labs(y = "Number of Compounds", x = "Subspecies + Ploidy")+
  theme_classic()+
  theme(legend.position = "none")

###########################################
## 2012 GC presence threshold defined ####
# seperate each subspecies and ploidy group to apply threshold
# Define the threshold of 20% 
threshold <- 0.80

# T_4n
md_OCG_GC_t4n <- subset(md_OCG_GC, md_OCG_GC$Subsp_ploidy=="T_4n") # 35
OCG_GC_2012_t4n <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md_OCG_GC_t4n)) # 20 of 74 variables

OCG_GC_2012_t4n[OCG_GC_2012_t4n == 0] <- NA
colSums(is.na(OCG_GC_2012_t4n)) #checking for null values since they are used for filtering

#Calculate proportion of NA values that are in each column
na_proportion_gc12_t4n <- colMeans(is.na(OCG_GC_2012_t4n))
print(na_proportion_gc12_t4n)

#identify which columns to keep
columns_to_keep_gc12_t4n <- na_proportion_gc12_t4n <= threshold 

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2012_t4n_subset <- OCG_GC_2012_t4n[, columns_to_keep_gc12_t4n] # 20 of 45 variables

# Replace NA values with zeroes
OCG_GC_2012_t4n_subset[is.na(OCG_GC_2012_t4n_subset)] <- 0
colSums(is.na(OCG_GC_2012_t4n_subset))

# T_2n
md_OCG_GC_t2n <- subset(md_OCG_GC, md_OCG_GC$Subsp_ploidy=="T_2n") # 74
OCG_GC_2012_t2n <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md_OCG_GC_t2n)) # 50 of 74 variables

OCG_GC_2012_t2n[OCG_GC_2012_t2n == 0] <- NA
colSums(is.na(OCG_GC_2012_t2n)) 

#Calculate proportion of NA values that are in each column
na_proportion_gc12_t2n <- colMeans(is.na(OCG_GC_2012_t2n))
print(na_proportion_gc12_t2n)

#identify which columns to keep
columns_to_keep_gc12_t2n <- na_proportion_gc12_t2n <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2012_t2n_subset <- OCG_GC_2012_t2n[, columns_to_keep_gc12_t2n] # 50 of 46 variables

# Replace NA values with zeroes
OCG_GC_2012_t2n_subset[is.na(OCG_GC_2012_t2n_subset)] <- 0
colSums(is.na(OCG_GC_2012_t2n_subset))

# V_4n
md_OCG_GC_v4n <- subset(md_OCG_GC, md_OCG_GC$Subsp_ploidy=="V_4n") # 22
OCG_GC_2012_v4n <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md_OCG_GC_v4n)) # 18 of 74 variables

OCG_GC_2012_v4n[OCG_GC_2012_v4n == 0] <- NA
colSums(is.na(OCG_GC_2012_v4n)) 

#Calculate proportion of NA values that are in each column
na_proportion_gc12_v4n <- colMeans(is.na(OCG_GC_2012_v4n))
print(na_proportion_gc12_v4n)

#identify which columns to keep
columns_to_keep_gc12_v4n <- na_proportion_gc12_v4n <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2012_v4n_subset <- OCG_GC_2012_v4n[, columns_to_keep_gc12_v4n] # 18 of 46 variables

# Replace NA values with zeroes
OCG_GC_2012_v4n_subset[is.na(OCG_GC_2012_v4n_subset)] <- 0
colSums(is.na(OCG_GC_2012_v4n_subset))

# V_2n
md_OCG_GC_v2n <- subset(md_OCG_GC, md_OCG_GC$Subsp_ploidy=="V_2n") # 28
OCG_GC_2012_v2n <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md_OCG_GC_v2n)) # 27 of 74 variables

OCG_GC_2012_v2n[OCG_GC_2012_v2n == 0] <- NA
colSums(is.na(OCG_GC_2012_v2n)) 

#Calculate proportion of NA values that are in each column
na_proportion_gc12_v2n <- colMeans(is.na(OCG_GC_2012_v2n))
print(na_proportion_gc12_v2n)

#identify which columns to keep
columns_to_keep_gc12_v2n <- na_proportion_gc12_v2n <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2012_v2n_subset <- OCG_GC_2012_v2n[, columns_to_keep_gc12_v2n] # 27 of 44 variables

# Replace NA values with zeroes
OCG_GC_2012_v2n_subset[is.na(OCG_GC_2012_v2n_subset)] <- 0
colSums(is.na(OCG_GC_2012_v2n_subset))

# W_4n
md_OCG_GC_w4n <- subset(md_OCG_GC, md_OCG_GC$Subsp_ploidy=="W_4n") # 65
OCG_GC_2012_w4n <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md_OCG_GC_w4n)) # 39 of 74 variables

OCG_GC_2012_w4n[OCG_GC_2012_w4n == 0] <- NA
colSums(is.na(OCG_GC_2012_w4n)) 

#Calculate proportion of NA values that are in each column
na_proportion_gc12_w4n <- colMeans(is.na(OCG_GC_2012_w4n))
print(na_proportion_gc12_w4n)

#identify which columns to keep
columns_to_keep_gc12_w4n <- na_proportion_gc12_w4n <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2012_w4n_subset <- OCG_GC_2012_w4n[, columns_to_keep_gc12_w4n] # 39 of 45

# Replace NA values with zeroes
OCG_GC_2012_w4n_subset[is.na(OCG_GC_2012_w4n_subset)] <- 0
colSums(is.na(OCG_GC_2012_w4n_subset))

# Combine all 5 subspecies groups into one GC 2012 dataset
OCG_GC_2012_subset <- bind_rows(
  OCG_GC_2012_t4n_subset,
  OCG_GC_2012_t2n_subset,
  OCG_GC_2012_v4n_subset,
  OCG_GC_2012_v2n_subset,
  OCG_GC_2012_w4n_subset
) # 154 of 62

# write GC 2012 data as csv
# write.csv(OCG_GC_2012_subset, "data_csv/OCG_GC_2012_thresholded.csv")

md_OCG_GC_2012 <- subset(md_OCG, row.names(md_OCG) %in% row.names(OCG_GC_2012_subset)) #112

## 2021 GC presence threshold defined ####
# T_4n
OCG_GC_2021_t4n <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md_OCG_GC_t4n)) # 15 of 74 variables

OCG_GC_2021_t4n[OCG_GC_2021_t4n == 0] <- NA
colSums(is.na(OCG_GC_2021_t4n)) #checking for null values since they are used for filtering

#Calculate proportion of NA values that are in each column
na_proportion_gc21_t4n <- colMeans(is.na(OCG_GC_2021_t4n))
print(na_proportion_gc21_t4n)

#identify which columns to keep
columns_to_keep_gc21_t4n <- na_proportion_gc21_t4n <= threshold 

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2021_t4n_subset <- OCG_GC_2021_t4n[, columns_to_keep_gc21_t4n] # 15 of 35 variables

# Replace NA values with zeroes
OCG_GC_2021_t4n_subset[is.na(OCG_GC_2021_t4n_subset)] <- 0
colSums(is.na(OCG_GC_2021_t4n_subset))

# T_2n
OCG_GC_2021_t2n <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md_OCG_GC_t2n)) # 24 of 74 variables

OCG_GC_2021_t2n[OCG_GC_2021_t2n == 0] <- NA
colSums(is.na(OCG_GC_2021_t2n)) 

#Calculate proportion of NA values that are in each column
na_proportion_gc21_t2n <- colMeans(is.na(OCG_GC_2021_t2n))
print(na_proportion_gc21_t2n)

#identify which columns to keep
columns_to_keep_gc21_t2n <- na_proportion_gc21_t2n <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2021_t2n_subset <- OCG_GC_2021_t2n[, columns_to_keep_gc21_t2n] # 24 of 40 variables

# Replace NA values with zeroes
OCG_GC_2021_t2n_subset[is.na(OCG_GC_2021_t2n_subset)] <- 0
colSums(is.na(OCG_GC_2021_t2n_subset))

# V_4n
OCG_GC_2021_v4n <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md_OCG_GC_v4n)) # 4 of 74 variables

OCG_GC_2021_v4n[OCG_GC_2021_v4n == 0] <- NA
colSums(is.na(OCG_GC_2021_v4n)) 

#Calculate proportion of NA values that are in each column
na_proportion_gc21_v4n <- colMeans(is.na(OCG_GC_2021_v4n))
print(na_proportion_gc21_v4n)

#identify which columns to keep
columns_to_keep_gc21_v4n <- na_proportion_gc21_v4n <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2021_v4n_subset <- OCG_GC_2021_v4n[, columns_to_keep_gc21_v4n] # 4 of 41 variables

# Replace NA values with zeroes
OCG_GC_2021_v4n_subset[is.na(OCG_GC_2021_v4n_subset)] <- 0
colSums(is.na(OCG_GC_2021_v4n_subset))

# V_2n
OCG_GC_2021_v2n <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md_OCG_GC_v2n)) # 1 of 74 variables

OCG_GC_2021_v2n[OCG_GC_2021_v2n == 0] <- NA
colSums(is.na(OCG_GC_2021_v2n)) 

#Calculate proportion of NA values that are in each column
na_proportion_gc21_v2n <- colMeans(is.na(OCG_GC_2021_v2n))

print(na_proportion_gc21_v2n)

#identify which columns to keep
columns_to_keep_gc21_v2n <- na_proportion_gc21_v2n <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2021_v2n_subset <- OCG_GC_2021_v2n[, columns_to_keep_gc21_v2n] # 1 of 23 variables

# Replace NA values with zeroes
OCG_GC_2021_v2n_subset[is.na(OCG_GC_2021_v2n_subset)] <- 0
colSums(is.na(OCG_GC_2021_v2n_subset))

# W_4n
OCG_GC_2021_w4n <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md_OCG_GC_w4n)) # 26 of 74 variables

OCG_GC_2021_w4n[OCG_GC_2021_w4n == 0] <- NA
colSums(is.na(OCG_GC_2021_w4n)) 

#Calculate proportion of NA values that are in each column
na_proportion_gc21_w4n <- colMeans(is.na(OCG_GC_2021_w4n))
print(na_proportion_gc21_w4n)

#identify which columns to keep
columns_to_keep_gc21_w4n <- na_proportion_gc21_w4n <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2021_w4n_subset <- OCG_GC_2021_w4n[, columns_to_keep_gc21_w4n] # 36 of 42

# Replace NA values with zeroes
OCG_GC_2021_w4n_subset[is.na(OCG_GC_2021_w4n_subset)] <- 0
colSums(is.na(OCG_GC_2021_w4n_subset))

# Combine all 5 subspecies groups into one GC 2021 dataset
OCG_GC_2021_subset <- bind_rows(
  OCG_GC_2021_t4n_subset,
  OCG_GC_2021_t2n_subset,
  OCG_GC_2021_v4n_subset,
  OCG_GC_2021_v2n_subset,
  OCG_GC_2021_w4n_subset
) # 70 of 47

# write GC 2021 data as csv
# write.csv(OCG_GC_2021_subset, "data_csv/OCG_GC_2021_thresholded.csv")

md_OCG_GC_2021 <- subset(md_OCG, row.names(md_OCG) %in% row.names(OCG_GC_2021_subset)) #70

## Full GC threshold defined #### 
OCG_GC_subset <- bind_rows(
  OCG_GC_2012_subset,
  OCG_GC_2021_subset
) # 224 of 66

# write csv with OCG GC data
# write.csv(OCG_GC_subset, "data_csv/OCG_GC_thresholded.csv")

md_OCG_GC <- subset(md_OCG, row.names(md_OCG) %in% row.names(OCG_GC_subset)) #224

## LC-MS presence threshold defined ####
# T_4n
md_OCG_LCMS_t4n <- subset(md_OCG_LCMS, md_OCG_LCMS$Subsp_ploidy=="T_4n") # 27
OCG_LCMS_t4n <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md_OCG_LCMS_t4n)) # 27 of 308

OCG_LCMS_t4n[OCG_LCMS_t4n == 0] <- NA
colSums(is.na(OCG_LCMS_t4n)) #checking for null values since they are used for filtering

#Calculate proportion of NA values that are in each column
na_proportion_lcms_t4n <- colMeans(is.na(OCG_LCMS_t4n))
print(na_proportion_lcms_t4n)

#identify which columns to keep
columns_to_keep_lcms_t4n <- na_proportion_lcms_t4n <= threshold 

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_LCMS_t4n_subset <- OCG_LCMS_t4n[, columns_to_keep_lcms_t4n] # 27 of 294 variables

# Replace NA values with zeroes
OCG_LCMS_t4n_subset[is.na(OCG_LCMS_t4n_subset)] <- 0
colSums(is.na(OCG_LCMS_t4n_subset))

# T_2n
md_OCG_LCMS_t2n <- subset(md_OCG_LCMS, md_OCG_LCMS$Subsp_ploidy=="T_2n") # 54
OCG_LCMS_t2n <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md_OCG_LCMS_t2n)) # 54 of 308

OCG_LCMS_t2n[OCG_LCMS_t2n == 0] <- NA
colSums(is.na(OCG_LCMS_t2n)) 

#Calculate proportion of NA values that are in each column
na_proportion_lcms_t2n <- colMeans(is.na(OCG_LCMS_t2n))
print(na_proportion_lcms_t2n)

#identify which columns to keep
columns_to_keep_lcms_t2n <- na_proportion_lcms_t2n <= threshold 

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_LCMS_t2n_subset <- OCG_LCMS_t2n[, columns_to_keep_lcms_t2n] # 54 of 278

# Replace NA values with zeroes
OCG_LCMS_t2n_subset[is.na(OCG_LCMS_t2n_subset)] <- 0
colSums(is.na(OCG_LCMS_t2n_subset))

# V_4n
md_OCG_LCMS_v4n <- subset(md_OCG_LCMS, md_OCG_LCMS$Subsp_ploidy=="V_4n") # 4
OCG_LCMS_v4n <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md_OCG_LCMS_v4n)) # 54 of 308

OCG_LCMS_v4n[OCG_LCMS_v4n == 0] <- NA
colSums(is.na(OCG_LCMS_v4n)) 

#Calculate proportion of NA values that are in each column
na_proportion_lcms_v4n <- colMeans(is.na(OCG_LCMS_v4n))
print(na_proportion_lcms_v4n)

#identify which columns to keep
columns_to_keep_lcms_v4n <- na_proportion_lcms_v4n <= threshold 

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_LCMS_v4n_subset <- OCG_LCMS_v4n[, columns_to_keep_lcms_v4n] # 4 of 290

# Replace NA values with zeroes
OCG_LCMS_v4n_subset[is.na(OCG_LCMS_v4n_subset)] <- 0
colSums(is.na(OCG_LCMS_v4n_subset))

# V_2n
md_OCG_LCMS_v2n <- subset(md_OCG_LCMS, md_OCG_LCMS$Subsp_ploidy=="V_2n") # 1
OCG_LCMS_v2n <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md_OCG_LCMS_v2n)) # 1 of 308

OCG_LCMS_v2n[OCG_LCMS_v2n == 0] <- NA
colSums(is.na(OCG_LCMS_v2n)) 

#Calculate proportion of NA values that are in each column
na_proportion_lcms_v2n <- colMeans(is.na(OCG_LCMS_v2n))
print(na_proportion_lcms_v2n)

#identify which columns to keep
columns_to_keep_lcms_v2n <- na_proportion_lcms_v2n <= threshold 

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_LCMS_v2n_subset <- OCG_LCMS_v2n[, columns_to_keep_lcms_v2n] # 1 of 234

# Replace NA values with zeroes
OCG_LCMS_v2n_subset[is.na(OCG_LCMS_v2n_subset)] <- 0
colSums(is.na(OCG_LCMS_v2n_subset))

# W_4n
md_OCG_LCMS_w4n <- subset(md_OCG_LCMS, md_OCG_LCMS$Subsp_ploidy=="W_4n") # 26
OCG_LCMS_w4n <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md_OCG_LCMS_w4n)) # 26 of 308

OCG_LCMS_w4n[OCG_LCMS_w4n == 0] <- NA
colSums(is.na(OCG_LCMS_w4n)) 

#Calculate proportion of NA values that are in each column
na_proportion_lcms_w4n <- colMeans(is.na(OCG_LCMS_w4n))
print(na_proportion_lcms_w4n)

#identify which columns to keep
columns_to_keep_lcms_w4n <- na_proportion_lcms_w4n <= threshold 

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_LCMS_w4n_subset <- OCG_LCMS_w4n[, columns_to_keep_lcms_w4n] # 26 of 299

# Replace NA values with zeroes
OCG_LCMS_w4n_subset[is.na(OCG_LCMS_w4n_subset)] <- 0
colSums(is.na(OCG_LCMS_w4n_subset))

# Combine all 5 subspecies groups into one LCMS dataset
OCG_LCMS_subset <- bind_rows(
  OCG_LCMS_t4n_subset,
  OCG_LCMS_t2n_subset,
  OCG_LCMS_v4n_subset,
  OCG_LCMS_v2n_subset,
  OCG_LCMS_w4n_subset
) # 112 of 308

# write LCMS data as csv
# write.csv(OCG_LCMS_subset, "data_csv/OCG_LCMS_thresholded.csv")

md_OCG_LCMS <- subset(md_OCG, row.names(md_OCG) %in% row.names(OCG_LCMS_subset)) #112

# drop levels
md_OCG_LCMS$Ecoregion <- droplevels(md_OCG_LCMS$Ecoregion)

########################################
## 2012 GC scaling ####
#Compound
data_normalized_2012 <- scale(OCG_GC_2012_subset) 
data_normalized_2012[is.na(data_normalized_2012)] <- 0
data_normalized_2012 <- data_normalized_2012[order(row.names(data_normalized_2012)),]

## 2021 GC scaling ####
#Compound 
data_normalized_2021 <- scale(OCG_GC_2021_subset) #1:70, 1:42
data_normalized_2021[is.na(data_normalized_2021)] <- 0
data_normalized_2021 <- data_normalized_2021[order(row.names(data_normalized_2021)),]

## Full GC scaling ####
#Compound
data_normalized_GC <- scale(OCG_GC_subset) #1:217, 1:54
data_normalized_GC[is.na(data_normalized_GC)] <- 0
data_normalized_GC <- data_normalized_GC[order(row.names(data_normalized_GC)),]

## LCMS scaling ####
#Compound
data_LCMS_normalized <- scale(OCG_LCMS_subset) #1:111, 1:302
data_LCMS_normalized[is.na(data_LCMS_normalized)] <- 0
data_LCMS_normalized <- data_LCMS_normalized[order(row.names(data_LCMS_normalized)),]

# Betadispersion test####
gc_betadisp <- betadisper(vegdist(data_normalized_GC, method = "euclidean"), md_OCG_GC$Ecoregion)
anova(gc_betadisp) 
TukeyHSD(gc_betadisp) 

lcms_betadisp <- betadisper(vegdist(data_LCMS_normalized, method = "euclidean"), md_OCG_LCMS$Ecoregion)
anova(lcms_betadisp) 
TukeyHSD(lcms_betadisp)

###########################################
## Full GC PCA ####
pca_GC <- prcomp(data_normalized_GC)
summary(pca_GC)
rownames(data_normalized_GC) == rownames(md_OCG_GC)

### PCA plot of full GC subspecies ploidy and year ####
# Figure 2C
plot(pca_GC$x[, 1], pca_GC$x[, 2],
     xlab="PC 1 (11.96%)", ylab="PC 2 (9.71%)", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md_OCG_GC$Subsp_ploidy],
     pch=c(17,19)[md_OCG_GC$Year],
     xlim = range(pca_GC$x[, 1], na.rm = TRUE),
     ylim = range(pca_GC$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
legend("bottomleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")
ordiellipse(pca_GC,groups = md_OCG_GC$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordiellipse(pca_GC,groups = md_OCG_GC$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordiellipse(pca_GC,groups = md_OCG_GC$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordiellipse(pca_GC,groups = md_OCG_GC$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordiellipse(pca_GC,groups = md_OCG_GC$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

### PCA plot of full GC ecoregion and year ####
gc_eco_palette <- c("#CF597E", "#E27170", "#E68969", "#E89A69", "#E79069", "#E9AD6D", "#EAC87F", "#EADA97", "#DDDE96", "#C3D78C", "#EAE29C", "#BED68A", "#8BC982", "#64C084", "#52BA88", "#1CA890", "#089392")

# Figure 2D 
plot(pca_GC$x[, 1], pca_GC$x[, 2],
     xlab="PC 1 (11.96%)", ylab="PC 2 (9.71%)", 
     col= gc_eco_palette[md_OCG_GC$Ecoregion],
     pch=c(17,19)[md_OCG_GC$Year],
     xlim = range(pca_GC$x[, 1], na.rm = TRUE),
     ylim = range(pca_GC$x[, 2], na.rm = TRUE))
legend("bottomleft", 
       legend=c("AZ/NM Plateau","AZ/NM Mountains", "Blue Mountains", "California Coastal Sage", "Central Basin & Range","Colorado Plateaus", "Columbia Plateau", "Idaho Batholith", "Middle Rockies", "Mojave Basin & Range", "Northern Basin & Range", "Northwestern Great Plains", "Snake River Plain", "Southern & Baja California", "Southern Rockies", "Wasatch & Uinta Mountains", "Wyoming Basin"),
       col= gc_eco_palette,
       pch=19,
       cex=0.6,
       bty = "n")
legend("topleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")
#### PERMANOVA for GC subspecies ploidy, ecoregion and year ####
pca_perm_GC_subsppl <- adonis2(data_normalized_GC ~ md_OCG_GC$Subsp_ploidy + md_OCG_GC$Year + md_OCG_GC$Ecoregion + md_OCG_GC$Site %in% md_OCG_GC$Ecoregion, data = GC_PCA_df, method = "euclidean", by = "term", permutations = 999)
pca_perm_GC_subsppl 

###########################################
## LCMS PCA ####
pca_LCMS <- prcomp(data_LCMS_normalized)
summary(pca_LCMS) # 12.99 9.44
rownames(data_LCMS_normalized) == rownames(md_OCG_LCMS)

### PCA plots of LCMS subspecies ploidy and year ####
# Figure 2A
plot(pca_LCMS$x[, 1], pca_LCMS$x[, 2],
     xlab="PC 1 (12.99%)", ylab="PC 2 (9.44%)", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md_OCG_LCMS$Subsp_ploidy],
     pch=c(17,19)[md_OCG_LCMS$Year],
     xlim = range(pca_LCMS$x[, 1], na.rm = TRUE),
     ylim = range(pca_LCMS$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown",'darkgreen','tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
legend("topright", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")
ordispider(pca_LCMS,groups = md_OCG_LCMS$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(pca_LCMS,groups = md_OCG_LCMS$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(pca_LCMS,groups = md_OCG_LCMS$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(pca_LCMS,groups = md_OCG_LCMS$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(pca_LCMS,groups = md_OCG_LCMS$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

### PCA plots of LCMS ecoregion and year ####
lcms_eco_palette <- c("#CF597E", "#E27170", "#E79069", "#E9AD6D", "#EAC87F", "#EAE29C", "#BED68A", "#8BC982", "#52BA88", "#1CA890", "#089392")

# Figure 2B
plot(pca_LCMS$x[, 1], pca_LCMS$x[, 2],
     xlab="PC 1 (12.99%)", ylab="PC 2 (9.44%)", 
     col= lcms_eco_palette[md_OCG_LCMS$Ecoregion],
     pch=c(17,19)[md_OCG_LCMS$Year],
     xlim = range(pca_LCMS$x[, 1], na.rm = TRUE),
     ylim = range(pca_LCMS$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("AZ/NM Plateau","AZ/NM Mountains","Central Basin & Range","Colorado Plateaus", "Columbia Plateaus", "Northern Basin & Range", "Northwestern Great Plains", "Snake River Plain", "Southern Rockies", "Wasatch & Uinta Mountains", "Wyoming Basin"),
       col= lcms_eco_palette,
       pch=19,
       cex=0.6,
       bty = "n")
legend("topright", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

#### PERMANOVA for LCMS subspecies, ecoregion, and year ####
pca_perm_lcms_subs <- adonis2(data_LCMS_normalized ~ md_OCG_LCMS$Subsp_ploidy + md_OCG_LCMS$Year + md_OCG_LCMS$Ecoregion + md_OCG_LCMS$Site %in% md_OCG_LCMS$Ecoregion, data = LCMS_PCA_df, by = "terms", method = "euclidean", permutations = 999)
pca_perm_lcms_subs

lcms_subsp.pw.r <- pairwise.adonis(data_LCMS_normalized, md_OCG_LCMS$Ecoregion, sim.method = "euclidean", perm = 9999)
lcms_subsp.pw.r 

###########################################
# PCA using the climate data ####
# read in climate data
climate <- read.csv("data_csv/ARTRspline_climate.csv", row.names = 1)

climate_subset <- subset(climate, rownames(climate) %in% md_OCG$`Plant Group`) # plant group = population

# subset the climate data to just the climate related variables 
climate_vars <- climate_subset[, which(colnames(climate_subset) == "elev"):which(colnames(climate_subset) == "mapmtcm")]

climate_vars <- climate_vars %>%
  rownames_to_column(var = "Population")

md_clim <- left_join(md_OCG, climate_vars, by = c("Plant Group" = "Population"))
rownames(md_clim) <- rownames(md_OCG)  # preserve original rownames

climate_vars <- climate_vars %>%
  column_to_rownames(var = "Population")

climate_pca <- prcomp(climate_vars, scale = TRUE, center = TRUE)
summary(climate_pca)

# Supplementary Figure 2 - climate PCA
biplot(climate_pca, xlab = "PC1 (62.7%)", ylab = "PC2 (13.8%)", cex = 0.65)

loadings <- climate_pca$rotation 

loadings_pc1_pc2 <- loadings[,1:2]

arrow_lengths <- sqrt(loadings_pc1_pc2[, "PC1"]^2 + loadings_pc1_pc2[, "PC2"]^2)

###########################################
# Spatial data ####
site_coords <- md_OCG[, c("Longitude", "Latitude")]
site_coords <- site_coords[order(row.names(site_coords)),]

## Partial Mantel for spatial and climate data for 2012 GC - T_2n ####
# GC
gc12_t2n_geo_samples <- intersect(rownames(OCG_GC_2012_t2n_subset), rownames(site_coords)) # 50 samples

# Subset and reorder both
gc12_t2n_geo <- OCG_GC_2012_t2n_subset[gc12_t2n_geo_samples, ]

# subset spatial df to just t_2n 2012 samples
spatial_df_GC12_t2n <- subset(site_coords, row.names(site_coords) %in% row.names(OCG_GC_2012_t2n_subset))
# spatial_gc12_t2n <- spatial_df[gc12_t2n_geo_samples, ]
spatial_df_GC12_t2n_dist <- vegdist(spatial_df_GC12_t2n, method = "euclidean")

gc12_t2n_dist <- vegdist(OCG_GC_2012_t2n_subset, method = "euclidean")

# subset climate metadata to just t_2n 2012 samples
md_clim_GC12_t2n <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC_2012_t2n_subset))
# clim_gc12_t2n <- subset(climate, rownames(climate) %in% md_clim_GC12_t2n$`Plant Group`) # plant group = population

clim_gc12_t2n <- md_clim_GC12_t2n[, which(colnames(md_clim_GC12_t2n) == "elev"):which(colnames(md_clim_GC12_t2n) == "mapmtcm")]
# climate_vars_t_2n <- climate_vars_t_2n %>%
#   rownames_to_column(var = "Population")

row.names(clim_gc12_t2n) == row.names(OCG_GC_2012_t2n_subset) # 50 samples
row.names(clim_gc12_t2n) == row.names(spatial_df_GC12_t2n)

clim_gc12_t2n_dist <- vegdist(clim_gc12_t2n, method = "euclidean")

# Partial Mantel test for 2012 GC T_2n samples controlling for climate
partial_mantel_gc12_t2n_spatial <- mantel.partial(gc12_t2n_dist, spatial_df_GC12_t2n_dist, clim_gc12_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc12_t2n_spatial) 

# Partial Mantel test for 2012 gc T_2n samples controlling for spatial
partial_mantel_gc12_t2n_clim <- mantel.partial(gc12_t2n_dist, clim_gc12_t2n_dist, spatial_df_GC12_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc12_t2n_clim) 

###########################################
## Partial Mantel for spatial and climate data for 2021 GC - T_2n ####
# 2021 GC
gc21_t2n_geo_samples <- intersect(rownames(OCG_GC_2021_t2n_subset), rownames(site_coords)) # 24 samples

# Subset and reorder both
gc21_t2n_geo <- OCG_GC_2021_t2n_subset[gc21_t2n_geo_samples, ]

spatial_gc21_t2n <- subset(site_coords, row.names(site_coords) %in% row.names(OCG_GC_2021_t2n_subset))
# spatial_gc21_t2n <- site_coords[gc21_t2n_geo_samples, ]
spatial_gc21_t2n_dist <- vegdist(spatial_gc21_t2n, method = "euclidean")

gc21_t2n_dist <- vegdist(OCG_GC_2021_t2n_subset, method = "euclidean")

# subset metadata to just t2n 2021 samples
md_clim_GC21_t2n <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC_2021_t2n_subset))

clim_gc21_t2n <- md_clim_GC21_t2n[, which(colnames(md_clim_GC21_t2n) == "elev"):which(colnames(md_clim_GC21_t2n) == "mapmtcm")]

row.names(clim_gc21_t2n) == row.names(OCG_GC_2021_t2n_subset) # 50 samples
row.names(clim_gc21_t2n) == row.names(spatial_gc21_t2n)

clim_gc21_t2n_dist <- vegdist(clim_gc21_t2n, method = "euclidean")

# Partial Mantel test for 2021 GC t_2n samples controlling for climate
partial_mantel_gc21_t2n_spatial <- mantel.partial(gc21_t2n_dist, spatial_gc21_t2n_dist, clim_gc21_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc21_t2n_spatial) 

# Partial Mantel test for 2021 GC t2n samples controlling for spatial
partial_mantel_gc21_t2n_clim <- mantel.partial(gc21_t2n_dist, clim_gc21_t2n_dist, spatial_gc21_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc21_t2n_clim) 

###########################################
## Partial Mantel for spatial and climate data for 2012 GC - W_4n ####
# 2012 GC
gc12_w4n_geo_samples <- intersect(rownames(OCG_GC_2012_w4n_subset), rownames(site_coords)) # 39 samples

# Subset and reorder both
# gc12_w4n_geo <- OCG_GC_2012_w4n_subset[gc12_w4n_geo_samples, ]
# spatial_gc12_w4n <- spatial_df[gc12_w4n_geo_samples, ]

spatial_df_GC12_w4n <- subset(site_coords, row.names(site_coords) %in% row.names(OCG_GC_2012_w4n_subset))
spatial_gc12_w4n_dist <- vegdist(spatial_df_GC12_w4n, method = "euclidean")

gc12_w4n_dist <- vegdist(OCG_GC_2012_w4n_subset, method = "euclidean")

# subset metadata to just w4n 2012 samples
md_clim_GC12_w4n <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC_2012_w4n_subset))

clim_gc12_w4n <- md_clim_GC12_w4n[, which(colnames(md_clim_GC12_w4n) == "elev"):which(colnames(md_clim_GC12_w4n) == "mapmtcm")]

row.names(clim_gc12_w4n) == row.names(OCG_GC_2012_w4n_subset) # 50 samples
row.names(clim_gc12_w4n) == row.names(spatial_df_GC12_w4n)

clim_gc12_w4n_dist <- vegdist(clim_gc12_w4n, method = "euclidean")

# Partial Mantel test for 2012 GC w_4n samples controlling for climate
partial_mantel_gc12_w4n_spatial <- mantel.partial(gc12_w4n_dist, spatial_gc12_w4n_dist, clim_gc12_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc12_w4n_spatial) 

# Partial Mantel test for 2012 GC w4n samples controlling for spatial
partial_mantel_gc12_w4n_clim <- mantel.partial(gc12_w4n_dist, clim_gc12_w4n_dist, spatial_gc12_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc12_w4n_clim) 

###########################################
## Partial Mantel for spatial and climate data for 2021 GC - W_4n ####
# 2021 GC
gc21_w4n_geo_samples <- intersect(rownames(OCG_GC_2021_w4n_subset), rownames(site_coords)) # 26 samples

# Subset and reorder both
gc21_w4n_geo <- OCG_GC_2021_w4n_subset[gc21_w4n_geo_samples, ]
# spatial_gc21_w4n <- spatial_df[gc21_w4n_geo_samples, ]

spatial_df_GC21_w4n <- subset(site_coords, row.names(site_coords) %in% row.names(OCG_GC_2021_w4n_subset))
spatial_gc21_w4n_dist <- vegdist(spatial_df_GC21_w4n, method = "euclidean")

gc21_w4n_dist <- vegdist(gc21_w4n_geo, method = "euclidean")

# subset metadata to just w4n 2021 samples
md_clim_GC21_w4n <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC_2021_w4n_subset))

clim_gc21_w4n <- md_clim_GC21_w4n[, which(colnames(md_clim_GC21_w4n) == "elev"):which(colnames(md_clim_GC21_w4n) == "mapmtcm")]

row.names(clim_gc21_w4n) == row.names(gc21_w4n_geo) # 50 samples
row.names(clim_gc21_w4n) == row.names(spatial_df_GC21_w4n)

clim_gc21_w4n_dist <- vegdist(clim_gc21_w4n, method = "euclidean")

# Partial Mantel test for 2021 gc w_4n samples controlling for climate
partial_mantel_gc21_w4n_spatial <- mantel.partial(gc21_w4n_dist, spatial_gc21_w4n_dist, clim_gc21_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc21_w4n_spatial)

# Partial Mantel test for 2021 gc w4n samples controlling for spatial
partial_mantel_gc21_w4n_clim <- mantel.partial(gc21_w4n_dist, clim_gc21_w4n_dist, spatial_gc21_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc21_w4n_clim) 

###########################################
## Partial Mantel for spatial and climate data for 2012 LC-MS - T_2n ####
# LCMS
# subset lcms to 2012 T_2n samples
md_clim_lcms_2012 <- subset(md_clim, md_clim$Year=="2012") #158
# md_clim_lcms_2012_t_2n <- subset(md_clim_lcms_2012, md_clim_lcms_2012$Subsp_ploidy=="T_2n") #30
md_clim_2012_lcms_t_2n <- subset(md_clim_lcms_2012, row.names(md_clim_lcms_2012) %in% row.names(OCG_LCMS_t2n_subset))  #30
OCG_LCMS_2012_t2n_subset <- subset(OCG_LCMS_t2n_subset, row.names(OCG_LCMS_t2n_subset) %in% row.names(md_clim_2012_lcms_t_2n)) # 30
lcms12_t2n_geo_samples <- intersect(rownames(OCG_LCMS_2012_t2n_subset), rownames(site_coords)) # 30 samples

# Subset and reorder both
lcms12_t2n_geo <- OCG_LCMS_2012_t2n_subset[lcms12_t2n_geo_samples, ]
# spatial_lcms12_t2n <- spatial_df[lcms12_t2n_geo_samples, ]

# subset spatial df to just t_2n 2012 samples
spatial_df_lcms12_t2n <- subset(site_coords, row.names(site_coords) %in% row.names(OCG_LCMS_2012_t2n_subset))
spatial_df_lcms12_t2n_dist <- vegdist(spatial_df_lcms12_t2n, method = "euclidean")

lcms12_t2n_dist <- vegdist(lcms12_t2n_geo, method = "euclidean")

# subset climate metadata to just t_2n 2012 samples
md_clim_lcms12_t2n <- subset(md_clim, row.names(md_clim) %in% row.names(lcms12_t2n_geo))

clim_lcms12_t2n <- md_clim_lcms12_t2n[, which(colnames(md_clim_lcms12_t2n) == "elev"):which(colnames(md_clim_lcms12_t2n) == "mapmtcm")]

row.names(clim_lcms12_t2n) == row.names(spatial_df_lcms12_t2n) 
row.names(clim_lcms12_t2n) == row.names(lcms12_t2n_geo)

clim_lcms12_t2n_dist <- vegdist(clim_lcms12_t2n, method = "euclidean")

# Partial Mantel test for 2012 lcms T_2n samples controlling for climate
partial_mantel_lcms12_t2n_spatial <- mantel.partial(lcms12_t2n_dist, spatial_df_lcms12_t2n_dist, clim_lcms12_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_lcms12_t2n_spatial) #significant

# Partial Mantel test for 2012 lcms T_2n samples controlling for spatial
partial_mantel_lcms12_t2n_clim <- mantel.partial(lcms12_t2n_dist, clim_lcms12_t2n_dist, spatial_df_lcms12_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_lcms12_t2n_clim) # significant

###########################################
## Partial Mantel test for spatial and climate data for 2021 LCMS- T_2n ####
# subset lcms to 2021 t_2n samples
md_clim_lcms_2021 <- subset(md_clim, md_clim$Year=="2021") #158
md_clim_lcms_2021_t_2n <- subset(md_clim_lcms_2021, md_clim_lcms_2021$Subsp_ploidy=="T_2n") #24
OCG_LCMS_2021_t2n_subset <- subset(OCG_LCMS_t2n_subset, row.names(OCG_LCMS_t2n_subset) %in% row.names(md_clim_lcms_2021_t_2n)) # 24
lcms21_t2n_geo_samples <- intersect(rownames(OCG_LCMS_2021_t2n_subset), rownames(site_coords)) # 24 samples

# Subset and reorder both
lcms21_t2n_geo <- OCG_LCMS_2021_t2n_subset[lcms21_t2n_geo_samples, ]
# spatial_lcms21_t2n <- spatial_df[lcms21_t2n_geo_samples, ]

# subset spatial df to just t_2n 2021 samples
spatial_df_lcms21_t2n <- subset(site_coords, row.names(site_coords) %in% row.names(OCG_LCMS_2021_t2n_subset))
spatial_df_lcms21_t2n_dist <- vegdist(spatial_df_lcms21_t2n, method = "euclidean")

lcms21_t2n_dist <- vegdist(lcms21_t2n_geo, method = "euclidean")

# subset climate metadata to just t_2n 2021 samples
md_clim_lcms21_t2n <- subset(md_clim, row.names(md_clim) %in% row.names(lcms21_t2n_geo))

clim_lcms21_t2n <- md_clim_lcms21_t2n[, which(colnames(md_clim_lcms21_t2n) == "elev"):which(colnames(md_clim_lcms21_t2n) == "mapmtcm")]

row.names(clim_lcms21_t2n) == row.names(spatial_df_lcms21_t2n) 
row.names(clim_lcms21_t2n) == row.names(lcms21_t2n_geo)

clim_lcms21_t2n_dist <- vegdist(clim_lcms21_t2n, method = "euclidean")

# Partial Mantel test for 2021 lcms T_2n samples controlling for climate
partial_mantel_lcms21_t2n_spatial <- mantel.partial(lcms21_t2n_dist, spatial_df_lcms21_t2n_dist, clim_lcms21_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_lcms21_t2n_spatial) # sig

# Partial Mantel test for 2021 lcms T_2n samples controlling for spatial
partial_mantel_lcms21_t2n_clim <- mantel.partial(lcms21_t2n_dist, clim_lcms21_t2n_dist, spatial_df_lcms21_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_lcms21_t2n_clim) # sig

###########################################
## Partial Mantel test for spatial and climate data for 2021 LCMS- W_4n ####
# subset lcms to 2021 w_4n samples
md_clim_lcms_2021_w_4n <- subset(md_clim_lcms_2021, md_clim_lcms_2021$Subsp_ploidy=="W_4n") #26
OCG_LCMS_2021_w4n_subset <- subset(OCG_LCMS_w4n_subset, row.names(OCG_LCMS_w4n_subset) %in% row.names(md_clim_lcms_2021_w_4n)) # 26
lcms21_w4n_geo_samples <- intersect(rownames(OCG_LCMS_2021_w4n_subset), rownames(site_coords)) # 26 samples

# Subset and reorder both
lcms21_w4n_geo <- OCG_LCMS_2021_w4n_subset[lcms21_w4n_geo_samples, ]
# spatial_lcms21_w4n <- spatial_df[lcms21_w4n_geo_samples, ]

# subset spatial df to just w_4n 2021 samples
spatial_df_lcms21_w4n <- subset(site_coords, row.names(site_coords) %in% row.names(lcms21_w4n_geo))
spatial_df_lcms21_w4n_dist <- vegdist(spatial_df_lcms21_w4n, method = "euclidean")

lcms21_w4n_dist <- vegdist(lcms21_w4n_geo, method = "euclidean")

# subset climate metadata to just w_4n 2021 samples
md_clim_lcms21_w4n <- subset(md_clim, row.names(md_clim) %in% row.names(lcms21_w4n_geo))

clim_lcms21_w4n <- md_clim_lcms21_w4n[, which(colnames(md_clim_lcms21_w4n) == "elev"):which(colnames(md_clim_lcms21_w4n) == "mapmtcm")]

row.names(clim_lcms21_w4n) == row.names(spatial_df_lcms21_w4n) 
row.names(clim_lcms21_w4n) == row.names(lcms21_w4n_geo)

clim_lcms21_w4n_dist <- vegdist(clim_lcms21_w4n, method = "euclidean")

# Partial Mantel test for 2021 lcms w_4n samples controlling for climate
partial_mantel_lcms21_w4n_spatial <- mantel.partial(lcms21_w4n_dist, spatial_df_lcms21_w4n_dist, clim_lcms21_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_lcms21_w4n_spatial) # not sig

# Partial Mantel test for 2021 w4n samples controlling for spatial
partial_mantel_lcms21_w4n_clim <- mantel.partial(lcms21_w4n_dist, clim_lcms21_w4n_dist, spatial_df_lcms21_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_lcms21_w4n_clim) # sig

###########################################
## PERMANOVA for climate variables for LC-MS 2012 T_2n ####
row.names(OCG_LCMS_2012_t2n_subset) == row.names(md_clim_lcms12_t2n)

cor(md_clim_lcms12_t2n[, c("pratio","mtwm","mtcm")], use = "pairwise.complete.obs")

permanova_lcms_2012_t2n <- adonis2(OCG_LCMS_2012_t2n_subset ~ pratio + mtcm + mtwm, data = md_clim_lcms12_t2n, permutations = 999, by = "margin")
print(permanova_lcms_2012_t2n) 

# PCA on T_2n 2012 LCMS samples for ordisurf plots ####
OCG_LCMS_2012_t2n_subset_scaled <- scale(OCG_LCMS_2012_t2n_subset) 
OCG_LCMS_2012_t2n_subset_scaled[is.na(OCG_LCMS_2012_t2n_subset_scaled)] <- 0
pca_lcms_2012_t2n <- prcomp(OCG_LCMS_2012_t2n_subset_scaled)
summary(pca_lcms_2012_t2n) # 19.70 9.30
rownames(OCG_LCMS_2012_t2n_subset_scaled) == rownames(md_clim_lcms12_t2n) # TRUE

#drop levels for ecoregion
md_clim_lcms12_t2n$Ecoregion <- droplevels(md_clim_lcms12_t2n$Ecoregion)

# palette for ecoregions
lcms_eco_t_2n_2012_palette <- c("#E79069","#E9AD6D",'#EAC87F','#EAE29C','#8BC982', '#1CA890')


# pratio
ordi_lcms_12_t_2n_prat <- ordisurf(pca_lcms_2012_t2n ~ md_clim_lcms12_t2n$pratio,
                              col = "black", main = "pratio", add = TRUE)
summary(ordi_lcms_12_t_2n_prat) # sig

# Figure 3 - panel C - pratio
plot(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
     type = "n",  # don't draw points yet
     xlab = "PC 1 (19.7%)", ylab = "PC 2 (9.29%)",
     xlim = range(pca_lcms_2012_t2n$x[, 1], na.rm = TRUE),
     ylim = range(pca_lcms_2012_t2n$x[, 2], na.rm = TRUE))
ordisurf(pca_lcms_2012_t2n$x[, 1:2], md_clim_lcms12_t2n$pratio,
         col = "black", add = TRUE)
points(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
       col = lcms_eco_t_2n_2012_palette[md_clim_lcms12_t2n$Ecoregion],
       pch = 19)
legend("topright",
       legend=c("Central Basin & Range","Colorado Plateaus","Columbia Plateau","Northern Basin & Range","Snake River Plain", "Wasatch & Uinta Mountains"),
       col= lcms_eco_t_2n_2012_palette,
       pch=19,
       cex=0.6,
       bty = "n")

# mtwm
ordi_lcms_12_t_2n_mtwm <- ordisurf(pca_lcms_2012_t2n, md_clim_lcms12_t2n$mtwm, 
                                   col = "black", main = "mtwm", add = TRUE)
summary(ordi_lcms_12_t_2n_mtwm)

# Figure 3 - panel B - mtwm 
plot(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
     type = "n",  # don't draw points yet
     xlab = "PC 1 (19.7%)", ylab = "PC 2 (9.29%)",
     xlim = range(pca_lcms_2012_t2n$x[, 1], na.rm = TRUE),
     ylim = range(pca_lcms_2012_t2n$x[, 2], na.rm = TRUE))
ordisurf(pca_lcms_2012_t2n$x[, 1:2], md_clim_lcms12_t2n$mtwm, 
         col = "black", add = TRUE)
points(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
       col = lcms_eco_t_2n_2012_palette[md_clim_lcms12_t2n$Ecoregion],
       pch = 19)
legend("topright",
       legend=c("Central Basin & Range","Colorado Plateaus","Columbia Plateau","Northern Basin & Range","Snake River Plain", "Wasatch & Uinta Mountains"),
       col= lcms_eco_t_2n_2012_palette,
       pch=19,
       cex=0.6,
       bty = "n")
# mtcm
ordi_lcms_12_t_2n_mtcm <- ordisurf(pca_lcms_2012_t2n, md_clim_lcms12_t2n$mtcm, 
                                   col = "black", main = "mtcm", add = TRUE)
summary(ordi_lcms_12_t_2n_mtcm)

# Figure 3 - panel A - mtcm
plot(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
     type = "n",  # don't draw points yet
     xlab = "PC 1 (19.7%)", ylab = "PC 2 (9.29%)",
     xlim = range(pca_lcms_2012_t2n$x[, 1], na.rm = TRUE),
     ylim = range(pca_lcms_2012_t2n$x[, 2], na.rm = TRUE))
ordisurf(pca_lcms_2012_t2n$x[, 1:2], md_clim_lcms12_t2n$mtcm, 
         col = "black", add = TRUE)
points(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
       col = lcms_eco_t_2n_2012_palette[md_clim_lcms12_t2n$Ecoregion],
       pch = 19)
legend("topright",
       legend=c("Central Basin & Range","Colorado Plateaus","Columbia Plateau","Northern Basin & Range","Snake River Plain", "Wasatch & Uinta Mountains"),
       col= lcms_eco_t_2n_2012_palette,
       pch=19,
       cex=0.6,
       bty = "n")


###########################################
## PERMANOVA for climate variables for LC-MS 2021 T_2n ####
rownames(OCG_LCMS_2021_t2n_subset) == rownames(md_clim_lcms21_t2n) # TRUE

# rda_lcms_2021_t2n <- rda(OCG_LCMS_2021_t2n_subset ~ map + mtwm + mtcm, data = md_clim_lcms21_t2n)
# anova(rda_lcms_2021_t2n, by = "term")
# vif.cca(rda_lcms_2021_t2n) 

permanova_lcms_2021_t2n <- adonis2(OCG_LCMS_2021_t2n_subset ~ pratio + mtcm + mtwm, data = md_clim_lcms21_t2n, permutations = 999, by = "margin")
print(permanova_lcms_2021_t2n)

cor(md_clim_lcms21_t2n[, c("map","mtwm","mtcm")], use = "pairwise.complete.obs")

# PCA on T_2n 2021 LCMS samples for ordisurf plot ####
# subset the lcms data to match the metadata
OCG_LCMS_2021_t2n_subset <- subset(OCG_LCMS_t2n_subset, row.names(OCG_LCMS_t2n_subset) %in% row.names(md_clim_lcms21_t2n)) # 24

OCG_LCMS_2021_t2n_subset_scaled <- scale(OCG_LCMS_2021_t2n_subset) #24
OCG_LCMS_2021_t2n_subset_scaled[is.na(OCG_LCMS_2021_t2n_subset_scaled)] <- 0

# remove columns with zero standard deviation
col_sds <- apply(OCG_LCMS_2021_t2n_subset_scaled, 2, sd, na.rm = TRUE)
OCG_LCMS_2021_t2n_subset_filtered <- OCG_LCMS_2021_t2n_subset_scaled[, col_sds > 0]

pca_lcms_2021_t2n <- prcomp(OCG_LCMS_2021_t2n_subset_filtered)
summary(pca_lcms_2021_t2n) 

rownames(OCG_LCMS_2021_t2n_subset_filtered) == rownames(md_clim_lcms21_t2n) # TRUE

#drop levels for ecoregion
md_clim_lcms21_t2n$Ecoregion <- droplevels(md_clim_lcms21_t2n$Ecoregion)

# palette for ecoregions
lcms_eco_t_2n_2021_palette <- c("#E79069",'#EAC87F','#EAE29C','#8BC982', '#1CA890', '#089392')

# pratio
ordi_lcms_21_t_2n_prat <- ordisurf(pca_lcms_2021_t2n, md_clim_lcms21_t2n$pratio, 
                                   col = "black", add = TRUE)
summary(ordi_lcms_21_t_2n_prat)

# MTcm
ordi_lcms_21_t_2n_mtcm <- ordisurf(pca_lcms_2021_t2n, md_clim_lcms21_t2n$mtcm, 
                                   col = "black", add = TRUE)
summary(ordi_lcms_21_t_2n_mtcm)

# MTWM
ordi_lcms_21_t_2n_mtwm <- ordisurf(pca_lcms_2021_t2n, md_clim_lcms21_t2n$mtwm, 
         col = "black", add = TRUE)
summary(ordi_lcms_21_t_2n_mtwm)

# Figure 3 - panel D - mtwm
plot(pca_lcms_2021_t2n$x[, 1], pca_lcms_2021_t2n$x[, 2],
     type = "n",  # don't draw points yet
     xlab = "PC 1 (20.89%)", ylab = "PC 2 (10.85%)",
     xlim = range(pca_lcms_2021_t2n$x[, 1], na.rm = TRUE),
     ylim = range(pca_lcms_2021_t2n$x[, 2], na.rm = TRUE))
ordisurf(pca_lcms_2021_t2n$x[, 1:2], md_clim_lcms21_t2n$mtwm, 
         col = "black", add = TRUE)
points(pca_lcms_2021_t2n$x[, 1], pca_lcms_2021_t2n$x[, 2],
       col = lcms_eco_t_2n_2021_palette[md_clim_lcms21_t2n$Ecoregion],
       pch = 19)
legend("bottomleft",
       legend = c("Central Basin & Range", "Columbia Plateau", "Northern Basin & Range","Snake River Plain", "Wasatch & Uinta Mountains", "Wyoming Basin"),
       col = lcms_eco_t_2n_2021_palette,
       pch = 19,
       cex = 0.6,
       bty = "n")

###########################################
# ANCOM BC for LCMS tridentata 2n Ecoregion####
table(md_OCG.LCMS.tri_2n$Ecoregion)
md_OCG.LCMS.tri.abc <- md_OCG.LCMS.tri_2n[!md_OCG.LCMS.tri_2n$Ecoregion %in% c("Wyoming Basin", "Colorado Plateaus"), ] # ecoregions with too small of sample sizes
md_OCG.LCMS.tri.abc$Ecoregion <- droplevels(md_OCG.LCMS.tri.abc$Ecoregion)
levels(md_OCG.LCMS.tri.abc$Ecoregion) # 5 ecoregions
OCG_LCMS_tri.abc <- subset(OCG_LCMS_t2n, row.names(OCG_LCMS_t2n) %in% row.names(md_OCG.LCMS.tri.abc)) 

rownames(OCG_LCMS_tri.abc) == rownames(md_OCG.LCMS.tri.abc) # TRUE
rounded_matrix.lcmstri <- as.matrix(OCG_LCMS_tri.abc)
rounded_matrix.lcmstri <- round(rounded_matrix.lcmstri)
rounded_matrix.lcmstri<- rounded_matrix.lcmstri
rounded_matrix.lcmstri[is.na(rounded_matrix.lcmstri)] <- 0
rounded_matrix.lcmstri <- as.data.frame(rounded_matrix.lcmstri)


# Create the tse object
assays.lcmstri = S4Vectors::SimpleList(counts = t(rounded_matrix.lcmstri)) 
smd.lcmstri = S4Vectors::DataFrame(md_OCG.LCMS.tri.abc)
tse.tri = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays.lcmstri, colData = smd.lcmstri)

# run model 
output.tri = ancombc2(data = tse.tri, assay_name = "counts", tax_level = NULL,
                      fix_formula = "Ecoregion", rand_formula = NULL,
                      p_adj_method = "fdr", pseudo_sens = TRUE,
                      prv_cut = 0.20, lib_cut = 1000, s0_perc = 0.05,
                      group = "Ecoregion", struc_zero = FALSE, neg_lb = FALSE,
                      alpha = 0.001, n_cl = 2, verbose = TRUE,
                      global = TRUE, pairwise = TRUE, 
                      dunnet = FALSE, trend = FALSE,
                      iter_control = list(tol = 1e-5, max_iter = 20, 
                                          verbose = FALSE),
                      em_control = list(tol = 1e-5, max_iter = 100),
                      lme_control = NULL, 
                      mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100), 
                      trend_control = NULL)

# # save model
# saveRDS(output.tri, "ancombc phytochemistry models/ancombc_lcms_loc_fdr.tri.RDS")
# output.tri <- readRDS("ancombc phytochemistry models/ancombc_lcms_loc_fdr.tri.RDS")

# Extract raw abundance data (counts)
abundance_data_tri <- assay(tse.tri, "counts")
View(abundance_data_tri)

# Extract p-values from ANCOM-BC 2 results
p_vals.tri <- as.data.frame(output.tri$res_global$q_val)

# identify sig compounds = 50
significant_compounds.tri <- output.tri$res_global %>%
  dplyr::filter(q_val < 0.05) %>%
  dplyr::select(taxon) %>%
  dplyr::pull()

abundance_sig.tri <- abundance_data_tri[significant_compounds.tri, , drop = FALSE] # subset abundance data to match significant compounds 

# log transform abundances 
abundance_log.tri <- log1p(abundance_sig.tri)  # log(1 + abundance)

# heat map with metadata call in 
annotation_col.tri <- data.frame(Ecoregion = md_OCG.LCMS.tri.abc$Ecoregion,
                                 row.names = rownames(md_OCG.LCMS.tri.abc))
rownames(annotation_col.tri) == colnames(abundance_log.tri)  # TRUE

Ecoregion_cluster <- hclust(dist(as.numeric(md_OCG.LCMS.tri.abc$Ecoregion)))

pheatmap(abundance_log.tri, 
         cluster_rows = TRUE, 
         cluster_cols = Ecoregion_cluster,
         show_colnames = TRUE,
         annotation_legend =  TRUE,  # Add metadata annotation
         annotation_names_col = FALSE,
         scale = "row",
         color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100),
         main = "Differentially Abundant Compounds - Tri ecoregion")


# extracting values 
res.tri <- output.tri$res

abundance_data <- output.tri$res_global

sig_abundance_data.tri <- abundance_data[output.tri$res_global$p_val < 0.05,]
sig_abundance_data.tri <- abundance_data[output.tri$res_global$diff_abn == "TRUE",]
rownames(sig_abundance_data.tri) <- as.vector(sig_abundance_data.tri[,1])
sig_abundance_data.tri <- sig_abundance_data.tri[,-1] # removes the first row with names that was moved in line before
# sig_abundance_data.tri <- sig_abundance_data.tri[, -c(1, 2, 4)] # removes diff column 1,2 4

# log_abundance_data.tri <- log(sig_abundance_data.tri + 1)
# 
# log_abundance_data.tri <- as.matrix(log_abundance_data.tri)  
# 
# log_abundance_data.tri <- log_abundance_data.tri[,1:length(md_OCG.LCMS.tri.abc$Location)] ## stops working here 
# 
# pheatmap(log_abundance_data.tri, cluster_rows = TRUE, cluster_cols = TRUE,
#          annotation_col = data.frame(Location = tse.tri$Location),
#          scale = "row", color = viridis::viridis(100))

# Rename columns containing "(Intercept)" to "AZ"
colnames(res.tri) <- gsub("\\(Intercept\\)", "EcoregionCentral Basin and Range", colnames(res.tri))
pair.tri <- output.tri$res_pair
dim(pair.tri) 
#256 61

str(res.tri)

#subset to the significant compounds between ecoregions
diff_cols <- grep("^diff_", names(res.tri), value = TRUE)
all_pair_sig.tri <- res.tri[rowSums(res.tri[, diff_cols] == TRUE) > 0, ] # 63 of 7 variables

## Figure for ANCOMBC LOCATION for tridentata 2n ####
lfc_long.tri <- all_pair_sig.tri %>%
  dplyr::select(taxon, starts_with("lfc_")) %>%  # Select compound and LFC columns
  tidyr::pivot_longer(cols = starts_with("lfc_"), 
                      names_to = "comparison", 
                      values_to = "lfc") %>%
  dplyr::mutate(comparison = gsub("lfc_", "", comparison)) %>%  # Clean comparison names
  dplyr::filter(!grepl("_", comparison)) %>%
  dplyr::mutate(comparison = factor(comparison, levels = c ("EcoregionCentral Basin and Range","EcoregionColumbia Plateau", "EcoregionNorthern Basin and Range", "EcoregionSnake River Plain", "EcoregionWasatch and Uinta Mountains")))# Keep rows without underscores

lfc_long.tri <- lfc_long.tri %>% 
  dplyr::mutate(lfc = ifelse(is.na(lfc), 0, lfc))

lfc_long.tri$comparison <- sub("^Ecoregion", "", lfc_long.tri$comparison) #removing location from the start of all location names 
# lfc_long.tri <- lfc_long.tri %>%
#   mutate(comparison = recode(comparison,
#                              "Central Basin and Range" = "Central B & R",
#                              "Columbia Plateau" = "Columbia Plateau",
#                              "Northern Basin and Range" = "Northern B & R",
#                              "Snake River Plain" = "Snake River Plain",
#                              "Wasatch and Uinta Mountains" = "Wasatch & Uinta Mtns"))

# Pivot the data to a matrix format
heatmap_data.tri <- reshape2::dcast(lfc_long.tri, taxon ~ comparison, value.var = "lfc")
rownames(heatmap_data.tri) <- heatmap_data.tri$taxon
heatmap_matrix.tri <- as.matrix(heatmap_data.tri[, -1])  # Remove the `taxon` column

pheatmap(heatmap_matrix.tri,
         color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100),  # Custom color scale
         clustering_distance_rows = "euclidean",  # Distance metric for rows
         clustering_distance_cols = "euclidean",  # Distance metric for columns
         clustering_method = "complete",         # Clustering method
         scale = "none",                        # No scaling of data
         angle = 0)
####

significant_compounds.tri <- lfc_long.tri %>%
  dplyr::select(taxon,comparison) 

significant_abundance_data.tri <- OCG_LCMS_tri.abc[,colnames(OCG_LCMS_tri.abc) %in% significant_compounds.tri$taxon]

significant_abundance_data.tri <- as.data.frame(t(significant_abundance_data.tri))
significant_abundance_data.tri[is.na(significant_abundance_data.tri)] <- 0

abundance_log.tri <- log1p(significant_abundance_data.tri)  # log(1 + abundance)

# Take the mean for columns grouped by state abbreviation
mean_abundance_tri <- abundance_log.tri %>%
  # Transpose and bind row names for easier manipulation
  t() %>%
  as.data.frame() %>%
  mutate(group = annotation_col.tri$Ecoregion) %>%
  group_by(group) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "group")  # Set group names as row names

mean_abundance_tri_t <- t(mean_abundance_tri)

# Figure 4 - HEATMAP for ANCOMBC results for tridentata 2n ecoregion
pheatmap(mean_abundance_tri_t,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "mediumorchid2", "mediumorchid4"))(100), 
         show_rownames = TRUE,
         show_colnames = TRUE,
         angle_col = 0,
         fontsize_row = 8,
         fontsize_col = 7)

###########################################
# Procrustes for LCMS against GC data ####
# subset original OCG uncleaned data with original GC data 
OCG_GCp <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(OCG_LCMS_3uL)) #108 of 74 variables 
OCG_LCMS_3uLp <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(OCG_GC)) #108 of 307 variables 

rownames(OCG_GCp) == rownames(OCG_LCMS_3uLp) #TRUE

# Replace NA with 0
OCG_LCMS_3uLp[is.na(OCG_LCMS_3uLp)] <- 0
OCG_GCp[is.na(OCG_GCp)] <- 0

# Procrustes test 
GC_v_LCMS <- protest(OCG_LCMS_3uLp, OCG_GCp, scale = TRUE) 

# Create the data frame for plotting (Correctly using $X and $Yrot):
procrustes_df <- data.frame(
  x = GC_v_LCMS$X[, 1],       # LCMS data (the one being transformed) - Dimension 1
  y = GC_v_LCMS$X[, 2],       # LCMS data - Dimension 2
  xend = GC_v_LCMS$Yrot[, 1],  # GC data - Dimension 1
  yend = GC_v_LCMS$Yrot[, 2],  # GC data - Dimension 2
  Sample = rownames(GC_v_LCMS$X) # Sample names (use either dataset's rownames)
)

## Supplementary Figure 3 #
ggplot(procrustes_df) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.6, color = "gray") +
  geom_point(aes(x = x, y = y, color = "LC-MS"), size = 2, shape = 19) +  # LCMS points
  geom_point(aes(x = xend, y = yend, color = "GC"), size = 2, shape = 17) + # Transformed GC points
  labs(x = "Procrustes axis 1", y = "Procrustes axis 2", color = "Data Source") +
  scale_color_manual(values = c("LC-MS" = "maroon", "GC" = "lightseagreen")) +
  theme_classic()+
  theme(legend.title = element_blank())

###########################################
#Creating a map of sites ####
# Read in shapefile
eco_shapefile <- st_read("data_csv/us_eco_l3_state_boundaries/us_eco_l3_state_boundaries.shp")

western_states <- c("Arizona", "New Mexico", "California", "Idaho", "Nevada", "Montana", "Utah", "Wyoming", "Colorado", "Oregon", "Washington")

eco_west <- eco_shapefile %>%
  filter(STATE_NAME %in% western_states)

# Turn site coordinates into an sf object
points_sf <- st_as_sf(md_OCG, coords = c("Longitude", "Latitude"), crs = 4326)

# Add the coordinates for the orchard location ID
orchard_location <- data.frame(
  Longitude = -115.998, # Note: Longitude is negative for west
  Latitude = 43.322
)

# turn ocg coor to sf
orchard_sf <- st_as_sf(orchard_location, coords = c("Longitude", "Latitude"), crs = 4326)

# subset to match subsetted shapefile
points_sf <- st_transform(points_sf, st_crs(eco_west))
orchard_sf <- st_transform(orchard_sf, st_crs(eco_west))

# Join points to ecoregions to find which ecoregion each site is in
site_ecos <- st_join(points_sf, eco_west)

# Get unique ecoregion names where sites are located
site_ecoregions <- unique(site_ecos$US_L3NAME)

# Number of highlighted ecoregions
n_highlight <- length(site_ecoregions)

# Generate palette
highlight_colors <- c("#E4796D", "#B2D387", "#E99F69", "grey80", "#CF679E", "#CBB47F", "#CCEBC5", "#EAC17A", "#6CC382", "#E1D8AB", "#089392", "grey80", "#D3DBC2", "#CF597E", "#FDDAEC", "#EAE29C", "#B3CDEA", "#29AD8E", "#A4BFAD")

# Create named color vector including grey for "Other"
final_colors <- c(highlight_colors, "grey80")
names(final_colors) <- c(site_ecoregions, "Other")

# Make highlight_group a factor to ensure order matches color names
eco_west$highlight_group <- factor(
  eco_west$highlight_group,
  levels = c(site_ecoregions, "Other")
)

# Figure 1 - Map of sites
ggplot() +
  geom_sf(data = eco_west, aes(fill = highlight_group), color = "white", size = 0.3) +
  scale_fill_manual(values = final_colors, name = "Ecoregions") +
  geom_sf(data = points_sf, color = "darkslategrey", size = 1.5) +
  geom_sf(data = orchard_sf, shape = 7, color = "white", size = 3) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf()







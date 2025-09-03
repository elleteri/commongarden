# APPROACH - 
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
md.OCG <- read.csv("data_csv/metadata_OCG.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #246 obs of 22 variables.
md.OCG <- md.OCG[order(row.names(md.OCG)),]

### Remove duplicates, negatve controls, and MTW.3.7.R_2012
md.OCG <- md.OCG[!rownames(md.OCG) %in% c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1', 'ORT.1.5_2012', 'NEG_8-28-21', 'NEG_10-2-20', 'MTW.3.7.R_2012'), ] #234 of 21 var

rownames(md.OCG) <- gsub("v[12]$", "", rownames(md.OCG), ignore.case = TRUE)

#make variables factor to plot and droplevels
md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Ecoregion","Description", "Plant Group")] <- lapply(md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Ecoregion","Description", "Plant Group")], as.factor)
md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Ecoregion","Description","Plant Group")] <- lapply(md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Ecoregion","Description", "Plant Group")], droplevels)
str(md.OCG)

#subset the md to only have observations from 2012 to avoid duplicates
md.OCG.2012 <- subset(md.OCG, md.OCG$Year=="2012") #158
str(md.OCG.2012)

#subset the md to only have observations from 2021 to avoid duplicates
md.OCG.2021 <- subset(md.OCG, md.OCG$Year=="2021") #76
str(md.OCG.2021)

##FULL CLEAN GC
OCG_GC <- read.csv("data_csv/OCG_GC_full_clean.csv", row.names = 1)#227 obs of 74 variables
OCG_GC <- OCG_GC[order(row.names(OCG_GC)),]
OCG_GC <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(md.OCG)) #217 of 74 variables
md.OCG.GC <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC)) #217 of 21 variables
OCG_GC[is.na(OCG_GC)] <- 0

## 2012 CLEAN GC
OCG_GC_2012 <- read.csv("data_csv/OCG_GC_2012_cleaned.csv", row.names = 1, check.names = FALSE) #157 obs of 74 variables
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

###########################################
# Alpha diversity ####
## GC alpha diversity ####
OCG.GC.shannon <- diversity(OCG_GC)
OCG.GC.ef <- exp(OCG.GC.shannon)
OCG.GC.ef.r <- round(OCG.GC.ef)

md.OCG.GC <- cbind(md.OCG.GC, GC_compounds = OCG.GC.ef.r)

# glmm.OCG.GC <- glmer(
#   GC_compounds ~ Year + Subsp_ploidy + Ecoregion + (1 | Plant),
#   family = Gamma(link = "log"),
#   data = md.OCG.GC
# )
# summary(glmm.OCG.GC) #year and subspecies sig

glmm.OCG.GC <- glmer(
  formula = GC_compounds ~ Year + Ecoregion + Subsp_ploidy + (1 | Plant),
  data = md.OCG.GC,
  family = Gamma(link = "log"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)  # increase from default (1e4)
  )
)
summary(glmm.OCG.GC)

plot(allEffects(glmm.OCG.GC))

# ocg_gc_aov <- aov(effective_species ~ Ecoregion + Subsp_ploidy + Year, data = md.OCG.GC)
# summary(ocg_gc_aov)
# 
# ocg_gc_k <- kruskal.test(effective_species ~ Ecoregion, data = md.OCG.GC)
# 
# shapiro.test(residuals(ocg_gc_aov))
# qqnorm(ocg_gc_aov$residuals)
# qqline(ocg_gc_aov$residuals)
# 
# post_hoc_gc <- pairwise.wilcox.test(md.OCG.GC$effective_species, md.OCG.GC$Ecoregion,
#                                     p.adjust.method = "BH")
# 
# p_mat <- as.data.frame(as.table(post_hoc_gc$p.value))
# 
# names(p_mat)[1:3] <- c("Group1", "Group2", "p_value")
# 
# p_summary <- p_mat %>%
#   filter(!is.na(p_value)) %>%
#   arrange(p_value)
# 
# print(p_summary)
# 
# plot(allEffects(glm.OCG.GC))

ggplot(data = md.OCG.GC, aes(Ecoregion, GC_compounds, fill = Ecoregion)) +
  geom_boxplot() +
  labs(y = "Number of Compounds", x = "Ecoregions")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))

ggplot(data = md.OCG.GC, aes(Subsp_ploidy, GC_compounds, fill =Subsp_ploidy)) +
  geom_boxplot() +
  scale_fill_manual(values = c("pink","brown","darkgreen",'tan','lightblue'))+
  labs(y = "Number of Compounds", x = "Subspecies + Ploidy")+
  theme_classic()+
  theme(legend.position = "none")

ggplot(md.OCG.GC, aes(Year, GC_compounds))+
  geom_boxplot(aes(group = Year, fill = Year))+
  scale_fill_manual(values  = c("maroon", "cyan"))+
  labs(y = "Number of Compounds")+
  theme_classic()+
  theme(legend.position = "none")

## LCMS alpha diversity ####
OCG.LCMS.shannon <- diversity(OCG_LCMS_3uL)
OCG.LCMS.ef <- exp(OCG.LCMS.shannon)
OCG.LCMS.ef.r <- round(OCG.LCMS.ef)

md.OCG.LCMS.3 <- cbind(md.OCG.LCMS.3, Compounds = OCG.LCMS.ef.r)
# glmm.OCG.LCMS <- glmer(
#   Compounds ~ Year + Subsp_ploidy + (1 | Plant),
#   family = Gamma(link = "log"),
#   data = md.OCG.LCMS.3
# )
glmm.OCG.LCMS <- glmer(
  formula = Compounds ~ Year + Ecoregion + Subsp_ploidy + (1 | Plant),
  data = md.OCG.LCMS.3,
  family = Gamma(link = "log"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)  # increase from default (1e4)
  )
)
summary(glmm.OCG.LCMS)

# ggplot(md.OCG.LCMS.3, aes(Ecoregion, Compounds))+
#   geom_boxplot(aes(group = Ecoregion, fill = Ecoregion))+
#   labs(y = "Number of Compounds")+
#   theme_classic()+
#   theme(legend.position = "none")

plot(allEffects(glmm.OCG.LCMS))

ggplot(md.OCG.LCMS.3, aes(Year, Compounds))+
  geom_boxplot(aes(group = Year, fill = Year))+
  scale_fill_manual(values = c("maroon","cyan"))+
  labs(y = "Number of Compounds")+
  theme_classic()+
  theme(legend.position = "none")

ggplot(md.OCG.LCMS.3, aes(Ecoregion, Compounds))+
  geom_boxplot(aes(group = Ecoregion, fill = Ecoregion))+
  labs(y = "Number of Compounds")+
  theme_classic()+
  theme(legend.position = "none")

ggplot(data = md.OCG.LCMS.3, aes(Subsp_ploidy, Compounds, fill = Subsp_ploidy)) +
  geom_boxplot() +
  scale_fill_manual(values = c("pink","brown","darkgreen",'tan','lightblue'))+
  labs(y = "Number of Compounds", x = "Subspecies + Ploidy")+
  theme_classic()+
  theme(legend.position = "none")

# ggplot(md.OCG.LCMS.3, aes(Ploidy, Compounds, fill = Ploidy))+
#   geom_boxplot()+
#   scale_fill_viridis_d(option = "plasma")+
#   labs(y = "Number of Compounds")+
#   theme_classic() +
#   theme(legend.position = "none")
# 
# ocg_lcms_aov <- aov(Compounds ~ Ecoregion + Subsp_ploidy + Year, data = md.OCG.LCMS.3)
# summary(ocg_lcms_aov)
# 
# ocg_lcms_k <- kruskal.test(Compounds ~ Ecoregion, data = md.OCG.LCMS.3)
# 
# shapiro.test(residuals(ocg_lcms_aov))
# qqnorm(ocg_lcms_aov$residuals)
# qqline(ocg_lcms_aov$residuals)
# 
# post_hoc_lcms <- pairwise.wilcox.test(md.OCG.LCMS.3$Compounds, md.OCG.LCMS.3$Ecoregion,
#                                       p.adjust.method = "BH")
# 
# p_mat <- as.data.frame(as.table(post_hoc_lcms$p.value))
# 
# names(p_mat)[1:3] <- c("Group1", "Group2", "p_value")
# 
# p_summary <- p_mat %>%
#   filter(!is.na(p_value)) %>%
#   arrange(p_value)
# 
# print(p_summary)

###########################################
# Cleaning for PCAs####

## 2012 GC presence threshold defined ####
# seperate each subspecies and ploidy group to apply threshold
# Define the threshold of 10% 
threshold <- 0.80

# T_4n
md.OCG.GC_t4n <- subset(md.OCG.GC, md.OCG.GC$Subsp_ploidy=="T_4n") # 35
OCG_GC_2012_t4n <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md.OCG.GC_t4n)) # 20 of 74 variables

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
md.OCG.GC_t2n <- subset(md.OCG.GC, md.OCG.GC$Subsp_ploidy=="T_2n") # 74
OCG_GC_2012_t2n <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md.OCG.GC_t2n)) # 50 of 74 variables

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
md.OCG.GC_v4n <- subset(md.OCG.GC, md.OCG.GC$Subsp_ploidy=="V_4n") # 22
OCG_GC_2012_v4n <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md.OCG.GC_v4n)) # 18 of 74 variables

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
md.OCG.GC_v2n <- subset(md.OCG.GC, md.OCG.GC$Subsp_ploidy=="V_2n") # 28
OCG_GC_2012_v2n <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md.OCG.GC_v2n)) # 27 of 74 variables

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
md.OCG.GC_w4n <- subset(md.OCG.GC, md.OCG.GC$Subsp_ploidy=="W_4n") # 65
OCG_GC_2012_w4n <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(md.OCG.GC_w4n)) # 39 of 74 variables

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

md.OCG.GC.2012 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC_2012_subset)) #112


## 2021 GC presence threshold defined ####
# T_4n
OCG_GC_2021_t4n <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md.OCG.GC_t4n)) # 15 of 74 variables

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
OCG_GC_2021_t2n <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md.OCG.GC_t2n)) # 24 of 74 variables

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
OCG_GC_2021_v4n <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md.OCG.GC_v4n)) # 4 of 74 variables

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
OCG_GC_2021_v2n <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md.OCG.GC_v2n)) # 1 of 74 variables

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
OCG_GC_2021_w4n <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md.OCG.GC_w4n)) # 26 of 74 variables

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

md.OCG.GC.2021 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC_2021_subset)) #70

## Full GC threshold defined #### 
OCG_GC_subset <- bind_rows(
  OCG_GC_2012_subset,
  OCG_GC_2021_subset
) # 224 of 66

md.OCG.GC <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_GC_subset)) #224

## LC-MS presence threshold defined ####
# T_4n
md.OCG.LCMS_t4n <- subset(md.OCG.LCMS.3, md.OCG.LCMS.3$Subsp_ploidy=="T_4n") # 27
OCG_LCMS_t4n <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.OCG.LCMS_t4n)) # 27 of 308

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
md.OCG.LCMS_t2n <- subset(md.OCG.LCMS.3, md.OCG.LCMS.3$Subsp_ploidy=="T_2n") # 54
OCG_LCMS_t2n <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.OCG.LCMS_t2n)) # 54 of 308

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
md.OCG.LCMS_v4n <- subset(md.OCG.LCMS.3, md.OCG.LCMS.3$Subsp_ploidy=="V_4n") # 4
OCG_LCMS_v4n <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.OCG.LCMS_v4n)) # 54 of 308

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
md.OCG.LCMS_v2n <- subset(md.OCG.LCMS.3, md.OCG.LCMS.3$Subsp_ploidy=="V_2n") # 1
OCG_LCMS_v2n <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.OCG.LCMS_v2n)) # 1 of 308

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
md.OCG.LCMS_w4n <- subset(md.OCG.LCMS.3, md.OCG.LCMS.3$Subsp_ploidy=="W_4n") # 26
OCG_LCMS_w4n <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(md.OCG.LCMS_w4n)) # 26 of 308

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

# Combine all 5 subspecies groups into one GC 2012 dataset
OCG_LCMS_subset <- bind_rows(
  OCG_LCMS_t4n_subset,
  OCG_LCMS_t2n_subset,
  OCG_LCMS_v4n_subset,
  OCG_LCMS_v2n_subset,
  OCG_LCMS_w4n_subset
) # 112 of 308

md.OCG.LCMS.3 <- subset(md.OCG, row.names(md.OCG) %in% row.names(OCG_LCMS_subset)) #112

# drop levels
md.OCG.LCMS.3$Ecoregion <- droplevels(md.OCG.LCMS.3$Ecoregion)

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

###########################################
## Full GC PCA ####
pca_GC <- prcomp(data_normalized_GC)
summary(pca_GC)
rownames(data_normalized_GC) == rownames(md.OCG.GC)

## PCA plot of full GC subspecies ploidy and year ####
plot(pca_GC$x[, 1], pca_GC$x[, 2],
     xlab="PC 1 (11.96%)", ylab="PC 2 (9.71%)", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC$Subsp_ploidy],
     pch=c(17,19)[md.OCG.GC$Year],
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
ordiellipse(pca_GC,groups = md.OCG.GC$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordiellipse(pca_GC,groups = md.OCG.GC$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordiellipse(pca_GC,groups = md.OCG.GC$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordiellipse(pca_GC,groups = md.OCG.GC$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordiellipse(pca_GC,groups = md.OCG.GC$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

## PCA plot of full GC ecoregion and year ####
gc_eco_palette <- c("#CF597E", "#E27170", "#E68969", "#E89A69", "#E79069", "#E9AD6D", "#EAC87F", "#EADA97", "#DDDE96", "#C3D78C", "#EAE29C", "#BED68A", "#8BC982", "#64C084", "#52BA88", "#1CA890", "#089392")


plot(pca_GC$x[, 1], pca_GC$x[, 2],
     xlab="PC 1 (11.96%)", ylab="PC 2 (9.71%)", 
     col= gc_eco_palette[md.OCG.GC$Ecoregion],
     pch=c(17,19)[md.OCG.GC$Year],
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
### PERMANOVA for GC subspecies ploidy, ecoregion and year ####
pca_scores_GC <- pca_GC$x[]
GC_PCA_df <- as.data.frame(pca_scores_GC)
pca_perm_GC_subsppl <- adonis2(pca_scores_GC ~ md.OCG.GC$Subsp_ploidy + md.OCG.GC$Year + md.OCG.GC$Ecoregion, data = GC_PCA_df, method = "euclidean", by = "margin")
pca_perm_GC_subsppl # subsppl: p = 0.001, R2 = 0.05532, for year R2 = 0.07923, p = 0.001 for ecoregions R2 = 0.12674 p = 0.001

OCG_GC_subsp.pw.r <- pairwise.adonis(pca_scores_GC, md.OCG.GC$Subsp_ploidy, sim.method = "euclidean")
OCG_GC_subsp.pw.r

###########################################
## LCMS PCA ####
pca_LCMS <- prcomp(data_LCMS_normalized)
summary(pca_LCMS) # 12.99 9.44
rownames(data_LCMS_normalized) == rownames(md.OCG.LCMS.3)

## PCA plots of LCMS subspecies ploidy and year ####
plot(pca_LCMS$x[, 1], pca_LCMS$x[, 2],
     xlab="PC 1 (12.99%)", ylab="PC 2 (9.44%)", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.LCMS.3$Subsp_ploidy],
     pch=c(17,19)[md.OCG.LCMS.3$Year],
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
ordispider(pca_LCMS,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(pca_LCMS,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(pca_LCMS,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(pca_LCMS,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(pca_LCMS,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

## PCA plots of LCMS ecoregion and year ####
lcms_eco_palette <- c("#CF597E", "#E27170", "#E79069", "#E9AD6D", "#EAC87F", "#EAE29C", "#BED68A", "#8BC982", "#52BA88", "#1CA890", "#089392")

plot(pca_LCMS$x[, 1], pca_LCMS$x[, 2],
     xlab="PC 1 (12.99%)", ylab="PC 2 (9.44%)", 
     col= lcms_eco_palette[md.OCG.LCMS.3$Ecoregion],
     pch=c(17,19)[md.OCG.LCMS.3$Year],
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

### PERMANOVA for LCMS subspecies, ecoregion, and year ####
pca_scores_LCMS <- pca_LCMS$x[,]
LCMS_PCA_df <- as.data.frame(pca_scores_LCMS)

pca_perm_lcms_subs <- adonis2(pca_scores_LCMS ~ md.OCG.LCMS.3$Subsp_ploidy + md.OCG.LCMS.3$Year + md.OCG.LCMS.3$Ecoregion, data = LCMS_PCA_df, by = "margin", method = "euclidean")
pca_perm_lcms_subs # subsp_pl p = 0.001, R2 = 0.08667; year p = 0.001, R2 = 0.03379, Ecoregion p = 0.001, R2 = 0.20327

lcms_subsp.pw.r <- pairwise.adonis(pca_scores_LCMS, md.OCG.LCMS.3$Subsp_ploidy, sim.method = "euclidean")
lcms_subsp.pw.r # sig between tridentatas, tridentata 4n and vaseyana and wyomingensis 4n, and tridentata 2n and vaseyana 2n and wyomingensis 2n

###########################################
# PCA using the climate data ####
# read in climate data
climate <- read.csv("data_csv/ARTRspline_climate.csv", row.names = 1)

# #remove non-numeric variables
# climate <- climate[, sapply(climate, is.numeric)]
# 
# #test for correlation between variables 
# str(climate)
# cor_climate_matrix <- cor(climate)
# 
# # Visualize the correlation matrix
# corrplot(cor_climate_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45, order = "hclust")

climate_subset <- subset(climate, rownames(climate) %in% md.OCG$`Plant Group`) # plant group = population

# subset the climate data to just the climate related variables 
climate_vars <- climate_subset[, which(colnames(climate_subset) == "elev"):which(colnames(climate_subset) == "mapmtcm")]
climate_vars <- climate_vars %>%
  rownames_to_column(var = "Population")

md_clim <- left_join(md.OCG, climate_vars, by = c("Plant Group" = "Population"))
rownames(md_clim) <- rownames(md.OCG)  # preserve original rownames

climate_vars <- climate_vars %>%
  column_to_rownames(var = "Population")

climate_scaled <- scale(climate_vars)
climate_pca <- prcomp(climate_scaled, center = TRUE)
summary(climate_pca)
biplot(climate_pca, scale = 0, xlab = "PC1 (58.65%)", ylab = "PC2 (15.53%)")

# Extract scores for populations
climate_scores <- as.data.frame(climate_pca$x[, 1:2])  # first 2 PCs
climate_scores$population <- rownames(climate_pca$x)

# Extract loadings (variable contributions)
climate_loadings <- as.data.frame(climate_pca$rotation[, 1:2])
climate_loadings$variable <- rownames(climate_loadings)

# # Scale loadings so arrows fit well on plot
mult <- 5   # adjust scaling factor
climate_loadings$PC1 <- climate_loadings$PC1 * mult
climate_loadings$PC2 <- climate_loadings$PC2 * mult

# Mark highlighted variables
climate_loadings$highlight <- ifelse(climate_loadings$variable %in% c("long", "mtcm", "mtwm"),
                             "highlight", "normal")

# Percent variance explained for axis labels
climate_pca_var <- summary(climate_pca)$importance[2, 1:2] * 100

# Build plot
ggplot() +
  # Plot populations
  # geom_point(data = climate_scores, aes(x = PC1, y = PC2), color = "black") +
  geom_text_repel(data = climate_scores, aes(x = PC1, y = PC2, label = population), size = 3) +
  theme(legend.position = "none") +
  
  # Arrows for normal variables
  geom_segment(data = subset(climate_loadings, highlight == "normal"),
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgreen", size = 0.4) +
  theme(legend.position = "none") +
  
  # Arrows for highlighted variables (thicker + longer)
  geom_segment(data = subset(climate_loadings, highlight == "highlight"),
               aes(x = 0, y = 0, xend = PC1*1.5, yend = PC2*1.5), # extend length
               arrow = arrow(length = unit(0.3, "cm")),
               color = "green", size = 1.2) +
  theme(legend.position = "none") +
  
  # Labels for variables
  geom_text_repel(data = climate_loadings,
                  aes(x = PC1, y = PC2, label = variable, color = highlight),
                  size = 3, segment.color = NA) +
  theme(legend.position = "none") +
  
  scale_color_manual(values = c("highlight" = "red", "normal" = "darkred")) +
  labs(x = paste0("PC1 (", round(climate_pca_var[1], 2), "%)"),
       y = paste0("PC2 (", round(climate_pca_var[2], 2), "%)")) +
  theme(legend.position = "none") +
  theme_classic()

# extract top 3 PCs
population_scores <- as.data.frame(climate_pca$x[, 1:3])
colnames(population_scores) <- c("climPC1", "climPC2", "climPC3")
population_scores$Population <- rownames(population_scores)

md_clim <- left_join(md_clim, population_scores, by = c("Plant Group" = "Population"))
rownames(md_clim) <- rownames(md.OCG)  # preserve original rownames

# remove any samples with NAs from PC1 column 
md_clim <- md_clim[!is.na(md_clim$climPC1), ] # 229 samples with climate data

# # GC
# # subset to make the GC PCA and the metadata with clim variables match
# pca_scores_GC_clim <- subset(GC_PCA_df, row.names(GC_PCA_df) %in% row.names(md_clim)) #224 of 74 variables
# md.OCG.GC_clim <- subset(md_clim, row.names(md_clim) %in% row.names(GC_PCA_df)) #224 of 27 variables
#
# # PERMANOVA with all the PCs, ecoregion, subspecies ploidy, 
# pca_perm_GC_clim <- adonis2(pca_scores_GC_clim ~ md.OCG.GC_clim$Subsp_ploidy + md.OCG.GC_clim$Ecoregion + md.OCG.GC_clim$Year + md.OCG.GC_clim$climPC1 + md.OCG.GC_clim$climPC2 + md.OCG.GC_clim$climPC3, data = pca_scores_GC_clim, method = "euclidean", by = "margin")
# pca_perm_GC_clim
# 
# # GC 2012
# # subset to make the GC 2012 PCA scores and the metadata with clim variables match
# pca_scores_GC12_clim <- subset(GC_12_PCA_df, row.names(GC_12_PCA_df) %in% row.names(md_clim)) #149 of 51variables
# md.OCG.GC12_clim <- subset(md_clim, row.names(md_clim) %in% row.names(GC_12_PCA_df)) #224 of 27 variables
# 
# # PERMANOVA with all the PCs, ecoregion, subspecies ploidy, 
# pca_perm_GC12_clim <- adonis2(pca_scores_GC12_clim ~ md.OCG.GC12_clim$Subsp_ploidy + md.OCG.GC12_clim$Ecoregion + md.OCG.GC12_clim$PC1 + md.OCG.GC12_clim$PC2 + md.OCG.GC12_clim$PC3, data = pca_scores_GC12_clim, method = "euclidean", by = "margin")
# pca_perm_GC12_clim
# 
# # GC 2021
# # subset to make the GC 2021 PCA scores and the metadata with clim variables match
# pca_scores_GC21_clim <- subset(GC_21_PCA_df, row.names(GC_21_PCA_df) %in% row.names(md_clim)) #70 of 45 variables
# md.OCG.GC21_clim <- subset(md_clim, row.names(md_clim) %in% row.names(GC_21_PCA_df)) #224 of 27 variables
# 
# # PERMANOVA with all the PCs, ecoregion, subspecies ploidy, 
# pca_perm_GC21_clim <- adonis2(pca_scores_GC21_clim ~ md.OCG.GC21_clim$Subsp_ploidy + md.OCG.GC21_clim$Ecoregion + md.OCG.GC21_clim$PC1 + md.OCG.GC21_clim$PC2 + md.OCG.GC21_clim$PC3, data = pca_scores_GC21_clim, method = "euclidean", by = "margin")
# pca_perm_GC21_clim
# 
# # LCMS
# # subset to make the LCMS PCA scores and the metadata with clim variables match
# pca_scores_lcms_clim <- subset(LCMS_PCA_df , row.names(LCMS_PCA_df) %in% row.names(md_clim)) #112 
# md.OCG.lcms_clim <- subset(md_clim, row.names(md_clim) %in% row.names(LCMS_PCA_df)) 
# 
# # PERMANOVA with all the PCs, ecoregion, subspecies ploidy, 
# pca_perm_lcms_clim <- adonis2(pca_scores_lcms_clim ~ md.OCG.lcms_clim$Subsp_ploidy + md.OCG.lcms_clim$Ecoregion + md.OCG.lcms_clim$Year + md.OCG.lcms_clim$PC1 + md.OCG.lcms_clim$PC2 + md.OCG.lcms_clim$PC3, data = pca_scores_lcms_clim, method = "euclidean", by = "margin")
# pca_perm_lcms_clim

# RDA using the climate related PCs with the chemistry data ####
# # GC 
# OCG_GC_subset_clim <- OCG_GC_subset[order(row.names(OCG_GC_subset)),]
# md.OCG.GC_clim <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC_subset_clim)) 
# OCG_GC_subset_clim <- subset(OCG_GC_subset_clim, row.names(OCG_GC_subset_clim) %in% row.names(md.OCG.GC_clim)) 
# rownames(OCG_GC_subset_clim) == rownames(md.OCG.GC_clim)
# 
# gc_dist_clim <- vegdist(OCG_GC_subset_clim, method = "euclidean")
# 
# rda_gc_clim <- vegan::rda(OCG_GC_subset_clim ~ Subsp_ploidy + Ecoregion + Year + PC1 + PC2 + PC3, data = md.OCG.GC_clim)
# anova(rda_gc_clim, by = "term") 
# 
# scores_gc_clim <- scores(rda_gc_clim, display = "sites")
# scores_gc_clim <- as.data.frame(scores_gc_clim)
# 
# ordiplot(rda_gc_clim, type = "t",display = "sites",cex = .6)
# 
# plot(scores_gc_clim[, 1], scores_gc_clim[, 2],
#      xlab="axis 1", ylab="axis 2", 
#      col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC_clim$Subsp_ploidy],
#      pch=c(17,19)[md.OCG.GC_clim$Year])
# legend("topright", 
#        legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
#        col= c("pink","brown",'darkgreen','tan','lightblue'),
#        pch=19,
#        cex=0.8,
#        bty = "n")
# legend("topleft", 
#        legend=c("2012","2021"),
#        col= "black",
#        pch=c(17,19),
#        cex=0.8,
#        bty = "n")
# 
# # GC 2012
# OCG_GC12_subset_clim <- OCG_GC_2012_subset[order(row.names(OCG_GC_2012_subset)),]
# md.OCG.GC12_clim <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC12_subset_clim)) 
# OCG_GC12_subset_clim <- subset(OCG_GC12_subset_clim, row.names(OCG_GC12_subset_clim) %in% row.names(md.OCG.GC12_clim)) 
# rownames(OCG_GC12_subset_clim) == rownames(md.OCG.GC12_clim)
# 
# rda_gc12_clim <- vegan::rda(OCG_GC12_subset_clim ~ Subsp_ploidy + Ecoregion + PC1 + PC2 + PC3, data = md.OCG.GC12_clim)
# anova(rda_gc12_clim, by = "term")
# 
# scores_gc12_clim <- scores(dbrda_gc12_clim, display = "sites")
# scores_gc12_clim <- as.data.frame(scores_gc12_clim)
# 
# ordiplot(dbrda_gc12_clim, type = "t",display = "sites",cex = .6)
# 
# plot(scores_gc12_clim[, 1], scores_gc12_clim[, 2],
#      xlab="axis 1", ylab="axis 2", 
#      col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC12_clim$Subsp_ploidy],
#      pch=c(17))
# legend("topleft", 
#        legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
#        col= c("pink","brown",'darkgreen','tan','lightblue'),
#        pch=19,
#        cex=0.8,
#        bty = "n")
# 
# # GC 2021
# OCG_GC21_subset_clim <- OCG_GC_2021_subset[order(row.names(OCG_GC_2021_subset)),]
# md.OCG.GC21_clim <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC21_subset_clim)) 
# OCG_GC21_subset_clim <- subset(OCG_GC21_subset_clim, row.names(OCG_GC21_subset_clim) %in% row.names(md.OCG.GC21_clim)) 
# rownames(OCG_GC21_subset_clim) == rownames(md.OCG.GC21_clim)
# 
# rda_gc21_clim <- vegan::rda(OCG_GC21_subset_clim ~ Subsp_ploidy + Ecoregion + PC1 + PC2 + PC3, data = md.OCG.GC21_clim)
# anova(rda_gc21_clim, by = "term")
# 
# scores_gc21_clim <- scores(rda_gc21_clim, display = "sites")
# scores_gc21_clim <- as.data.frame(scores_gc21_clim)
# 
# ordiplot(rda_gc21_clim, type = "t",display = "sites",cex = .6)
# 
# plot(scores_gc21_clim[, 1], scores_gc21_clim[, 2],
#      xlab="axis 1", ylab="axis 2", 
#      col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC21_clim$Subsp_ploidy],
#      pch=c(17))
# legend("bottomleft", 
#        legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
#        col= c("pink","brown",'darkgreen','tan','lightblue'),
#        pch=17,
#        cex=0.8,
#        bty = "n")
# 
# # LCMS
# OCG_LCMS_subset_clim <- OCG_LCMS_subset[order(row.names(OCG_LCMS_subset)),]
# md.OCG.LCMS_clim <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_LCMS_subset_clim)) 
# OCG_LCMS_subset_clim <- subset(OCG_LCMS_subset_clim, row.names(OCG_LCMS_subset_clim) %in% row.names(md.OCG.LCMS_clim)) 
# rownames(OCG_LCMS_subset_clim) == rownames(md.OCG.LCMS_clim)
# 
# rda_lcms_clim <- vegan::rda(OCG_LCMS_subset_clim ~ Subsp_ploidy + Ecoregion + Year + PC1 + PC2 + PC3, data = md.OCG.LCMS_clim)
# anova(rda_lcms_clim, by = "term")
# 
# scores_lcms_clim <- scores(dbrda_lcms_clim, display = "sites")
# scores_lcms_clim <- as.data.frame(scores_lcms_clim)
# 
# arrow_scores_lcms_clim <- scores(dbrda_lcms_clim, display = "bp")
# row.names(arrow_scores_lcms_clim)
# 
# # # only plot the arrows from PC 2 and PC 3 
# # arrow_scores_lcms_clim <- arrow_scores_lcms_clim[c("PC2", "PC3"), ]
# 
# ordiplot(dbrda_lcms_clim, type = "t",display = "sites",cex = .6)
# 
# plot(scores_lcms_clim[, 1], scores_lcms_clim[, 2],
#      xlab="axis 1", ylab="axis 2", 
#      col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.LCMS_clim$Subsp_ploidy],
#      pch=c(17,19)[md.OCG.LCMS_clim$Year])
# legend("topleft", 
#        legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
#        col= c("pink","brown",'darkgreen','tan','lightblue'),
#        pch=19,
#        cex=0.8,
#        bty = "n")
# legend("bottomleft", 
#        legend=c("2012","2021"),
#        col= "black",
#        pch=c(17,19),
#        cex=0.8,
#        bty = "n")
# 
# 

###########################################
# Mantel using spatial matrix against GC and LCMS #### 
# Extract coordinates matrix
site_coords <- md.OCG[, c("Longitude", "Latitude")]
site_coords <- site_coords[order(row.names(site_coords)),]

# match row names for mantel test
# LCMS
lcms_geo_samples <- intersect(rownames(OCG_LCMS_subset), rownames(site_coords)) # 41 samples

# Subset and reorder both
lcms_geo <- OCG_LCMS_subset[lcms_geo_samples, ]
site_coords_lcms <- site_coords[lcms_geo_samples, ]

# create geographic distance matrix
library(geodist)
geo_dist_lcms <- geodist(site_coords_lcms, measure = "geodesic")
geo_dist_lcms <- as.dist(geo_dist_lcms)

# create lcms data matrix
lcms_dist <- vegdist(lcms_geo, method = "euclidean")

# mantel test 
mantel_geo_lcms <- mantel(geo_dist_lcms, lcms_dist, method = "spearman", permutations = 9999)
print(mantel_geo_lcms) #significant. Mantel statistic r: 0.3979; Significance: 1e-04 

# GC dist matrix
# Subset and reorder both
gc_geo_samples <- intersect(rownames(OCG_GC), rownames(site_coords)) # 224 samples
gc_geo <- OCG_GC[gc_geo_samples, ]
site_coords_gc <- site_coords[gc_geo_samples, ]

# create gc distance matrix
gc_dist <- vegdist(gc_geo, method = "euclidean")

# create geographic distance matrix for GC samples
geo_dist_gc <- geodist(site_coords_gc, measure = "geodesic")
geo_dist_gc <- as.dist(geo_dist_gc)

# mantel test 
mantel_geo_gc <- mantel(geo_dist_gc, gc_dist, method = "spearman", permutations = 9999)
print(mantel_geo_gc) # not significant Mantel statistic r: -0.01069; Significance: 0.6239 

# 2012 GC 
# Subset and reorder both
gc12_geo_samples <- intersect(rownames(OCG_GC_2012_subset), rownames(site_coords)) 
gc12_geo <- OCG_GC_2012_subset[gc12_geo_samples, ]
site_coords_gc12 <- site_coords[gc12_geo_samples, ]

# create gc distance matrix
gc12_dist <- vegdist(gc12_geo, method = "euclidean")

# create geographic distance matrix for GC samples
geo_dist_gc12 <- geodist(site_coords_gc12, measure = "geodesic")
geo_dist_gc12 <- as.dist(geo_dist_gc12)

# mantel test 
mantel_geo_gc12 <- mantel(geo_dist_gc12, gc12_dist, method = "spearman", permutations = 9999)
print(mantel_geo_gc12) # not significant Mantel statistic r: 0.006646; Significance: 0.3994 

# 2021 GC
# Subset and reorder both
gc21_geo_samples <- intersect(rownames(OCG_GC_2021_subset), rownames(site_coords)) 
gc21_geo <- OCG_GC_2021_subset[gc21_geo_samples, ]
site_coords_gc21 <- site_coords[gc21_geo_samples, ]

# create gc distance matrix
gc21_dist <- vegdist(gc21_geo, method = "euclidean")

# create geographic distance matrix for GC samples
geo_dist_gc21 <- geodist(site_coords_gc21, measure = "geodesic")
geo_dist_gc21 <- as.dist(geo_dist_gc21)

# mantel test 
mantel_geo_gc21 <- mantel(geo_dist_gc21, gc21_dist, method = "spearman", permutations = 9999)
print(mantel_geo_gc21) #significant. Mantel statistic r: -0.02557; Significance: 0.6931

# correlation plot using mantel test results 
# convert to distance matrices 
geo_vec <- as.vector(geo_dist_lcms)
lcms_vec <- as.vector(lcms_dist)

# create a dataframe
lcms_geo_df <- data.frame(geo = geo_vec, lcms = lcms_vec)

ggplot(lcms_geo_df, aes(x = geo, y = lcms)) +
  geom_point(alpha = 0.3, color = "#FF73EE") +
  geom_smooth(method = "lm", se = TRUE, color = "#ba0033") +
  labs(title = "Mantel Test: Geographic vs. LCMS Distance",
       x = "Geographic Distance (based on population)",
       y = "LCMS Euclidean Distance") +
  theme_minimal() +
  theme( plot.title = element_text(hjust = 0.5))
###########################################
## RDA for spatial and climate data for GC - T_2n ####
md.clim_GC_t2n <- subset(md_clim, md.OCG.GC$Subsp_ploidy=="T_2n") # 71
OCG_GC_t2n <- subset(OCG_GC_subset_clim, row.names(OCG_GC_subset_clim) %in% row.names(md.clim_GC_t2n)) # 71
md.clim_GC_t2n <- subset(md.clim_GC_t2n, row.names(md.clim_GC_t2n) %in% row.names(OCG_GC_t2n)) # 71
OCG_GC_t2n[is.na(OCG_GC_t2n)] <- 0
OCG_GC_t2n <- OCG_GC_t2n[order(row.names(OCG_GC_t2n)),]
rownames(OCG_GC_t2n) == rownames(md.clim_GC_t2n)

# create geographic distance matrix for GC T_2n samples
library(geodist)
geo_dist_gc_t2n <- geodist(md.clim_GC_t2n[, c("Longitude", "Latitude")], measure = "geodesic")

rda_gc_clim_spatial<- vegan::rda(OCG_GC_t2n ~ PC1 + PC2 + PC3 + geo_dist_gc_t2n, data = md.clim_GC_t2n)
anova(rda_gc_clim_spatial, by = "term")

# using the three climate variables from climate as predictors
rda_gc_clim_spatial2 <- vegan::rda(OCG_GC_t2n ~ long + mtwm + mtcm + geo_dist_gc_t2n, data = md.clim_GC_t2n)
anova(rda_gc_clim_spatial2, by = "term")
# vegan::vif.cca(rda_gc_clim_spatial)

scores_gc_clim_spatial <- scores(rda_gc_clim_spatial, display = "sites")
scores_gc_clim_spatial <- as.data.frame(scores_gc_clim_spatial)

ordiplot(rda_gc_clim_spatial, type = "t",display = "sites",cex = .6)

plot(scores_gc_clim_spatial[, 1], scores_gc_clim_spatial[, 2],
     xlab="axis 1", ylab="axis 2", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.clim_GC_t2n$Subsp_ploidy],
     pch=c(17,19)[md.clim_GC_t2n$Year])
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown",'darkgreen','tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")

# ## PCA for spatial data ####
# # Extract coordinates matrix
# site_coords <- md.OCG[, c("Longitude", "Latitude")]
# site_coords <- site_coords[order(row.names(site_coords)),]
# 
# site_coords_scaled <- scale(site_coords)
# spatial_pca <- prcomp(site_coords_scaled, center = TRUE)
# summary(spatial_pca)
# biplot(spatial_pca, scale = 0, xlab = "PC1 (75.5%)", ylab = "PC2 (24.5%)")
# 
# # extract 2 PCs
# spatial_df <- as.data.frame(spatial_pca$x[, 1:2])

site_coords <- md.OCG[, c("Longitude", "Latitude")]
site_coords <- site_coords[order(row.names(site_coords)),]

# # rename "Plant Group" column to "Population"
# names(site_coords)[names(site_coords) == "Plant Group"] <- "Population"

# populations_coords <- site_coords %>%
#   group_by(Population) %>%
#   summarise(Longitude = mean(Longitude, na.rm = TRUE),
#             Latitude = mean(Latitude, na.rm = TRUE),
#             .groups = "drop")
# 
# D_km <- as.matrix(geodist(pops[, c("lon","lat")], measure = "geodesic")) / 1000
# rownames(D_km) <- colnames(D_km) <- pops$population



## Partial Mantel for spatial and climate data for 2012 GC - T_2n ####
# match row names for mantel test
# GC
gc12_t2n_geo_samples <- intersect(rownames(OCG_GC_2012_t2n_subset), rownames(site_coords)) # 50 samples

# Subset and reorder both
gc12_t2n_geo <- OCG_GC_2012_t2n_subset[gc12_t2n_geo_samples, ]

# subset spatial df to just t_2n 2012 samples
spatial_df_GC12_t2n <- subset(site_coords, row.names(site_coords) %in% row.names(OCG_GC_2012_t2n_subset))
# spatial_gc12_t2n <- spatial_df[gc12_t2n_geo_samples, ]
spatial_df_GC12_t2n_dist <- vegdist(spatial_df_GC12_t2n, method = "euclidean")

# # create geographic distance matrix from the first two PCs in a spatial PCA
# # library(geodist)
# geo_dist_gc12_t2n <- geodist(site_coords_gc12_t2n, measure = "geodesic")
# geo_dist_gc12_t2n <- as.dist(geo_dist_gc12_t2n)

gc12_t2n_dist <- vegdist(OCG_GC_2012_t2n_subset, method = "euclidean")

# subset climate metadata to just t_2n 2012 samples
md.clim_GC12_t2n <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC_2012_t2n_subset))
# clim_gc12_t2n <- subset(climate, rownames(climate) %in% md.clim_GC12_t2n$`Plant Group`) # plant group = population

# subset the climate data to just the climate related variables 
clim_gc12_t2n <- md.clim_GC12_t2n[, which(colnames(md.clim_GC12_t2n) == "elev"):which(colnames(md.clim_GC12_t2n) == "mapmtcm")]
# climate_vars_t_2n <- climate_vars_t_2n %>%
#   rownames_to_column(var = "Population")

row.names(clim_gc12_t2n) == row.names(OCG_GC_2012_t2n_subset) # 50 samples
row.names(clim_gc12_t2n) == row.names(spatial_df_GC12_t2n)

clim_gc12_t2n_dist <- vegdist(clim_gc12_t2n, method = "euclidean")

# partial mantel test for 2012 gc T_2n samples controlling for climate
partial_mantel_gc12_t2n_spatial <- mantel.partial(gc12_t2n_dist, spatial_df_GC12_t2n_dist, clim_gc12_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc12_t2n_spatial) # Mantel statistic r: ; not significant

# partial mantel test for 2012 gc T_2n samples controlling for spatial
partial_mantel_gc12_t2n_clim <- mantel.partial(gc12_t2n_dist, clim_gc12_t2n_dist, spatial_df_GC12_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc12_t2n_clim) # Mantel statistic r: -0.08712; not significant

###########################################
## Partial Mantel for spatial and climate data for 2021 GC - T_2n ####
# match row names for mantel test
# 2021 GC
gc21_t2n_geo_samples <- intersect(rownames(OCG_GC_2021_t2n_subset), rownames(site_coords)) # 24 samples

# Subset and reorder both
gc21_t2n_geo <- OCG_GC_2021_t2n_subset[gc21_t2n_geo_samples, ]

spatial_gc21_t2n <- subset(site_coords, row.names(site_coords) %in% row.names(OCG_GC_2021_t2n_subset))
# spatial_gc21_t2n <- site_coords[gc21_t2n_geo_samples, ]
spatial_gc21_t2n_dist <- vegdist(spatial_gc21_t2n, method = "euclidean")

gc21_t2n_dist <- vegdist(OCG_GC_2021_t2n_subset, method = "euclidean")

# subset metadata to just t2n 2021 samples
md.clim_GC21_t2n <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC_2021_t2n_subset))

clim_gc21_t2n <- md.clim_GC21_t2n[, which(colnames(md.clim_GC21_t2n) == "elev"):which(colnames(md.clim_GC21_t2n) == "mapmtcm")]

row.names(clim_gc21_t2n) == row.names(OCG_GC_2021_t2n_subset) # 50 samples
row.names(clim_gc21_t2n) == row.names(spatial_gc21_t2n)

clim_gc21_t2n_dist <- vegdist(clim_gc21_t2n, method = "euclidean")

# partial mantel test for 2021 gc t_2n samples controlling for climate
partial_mantel_gc21_t2n_spatial <- mantel.partial(gc21_t2n_dist, spatial_gc21_t2n_dist, clim_gc21_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc21_t2n_spatial) # Mantel statistic r: -0.07211; not significant

# partial mantel test for 2021 gc t2n samples controlling for spatial
partial_mantel_gc21_t2n_clim <- mantel.partial(gc21_t2n_dist, clim_gc21_t2n_dist, spatial_gc21_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc21_t2n_clim) # Mantel statistic r: -0.01105; not significant

###########################################
## Partial Mantel for spatial and climate data for 2012 GC - W_4n ####
# match row names for mantel test
# 2012 GC
gc12_w4n_geo_samples <- intersect(rownames(OCG_GC_2012_w4n_subset), rownames(site_coords)) # 39 samples

# Subset and reorder both
# gc12_w4n_geo <- OCG_GC_2012_w4n_subset[gc12_w4n_geo_samples, ]
# spatial_gc12_w4n <- spatial_df[gc12_w4n_geo_samples, ]

spatial_df_GC12_w4n <- subset(site_coords, row.names(site_coords) %in% row.names(OCG_GC_2012_w4n_subset))
spatial_gc12_w4n_dist <- vegdist(spatial_df_GC12_w4n, method = "euclidean")

# # create geographic distance matrix
# # library(geodist)
# geo_dist_gc12_w4n <- geodist(site_coords_gc12_w4n, measure = "geodesic")
# geo_dist_gc12_w4n <- as.dist(geo_dist_gc12_w4n)

gc12_w4n_dist <- vegdist(OCG_GC_2012_w4n_subset, method = "euclidean")

# subset metadata to just w4n 2012 samples
md.clim_GC12_w4n <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC_2012_w4n_subset))

clim_gc12_w4n <- md.clim_GC12_w4n[, which(colnames(md.clim_GC12_w4n) == "elev"):which(colnames(md.clim_GC12_w4n) == "mapmtcm")]

row.names(clim_gc12_w4n) == row.names(OCG_GC_2012_w4n_subset) # 50 samples
row.names(clim_gc12_w4n) == row.names(spatial_df_GC12_w4n)

clim_gc12_w4n_dist <- vegdist(clim_gc12_w4n, method = "euclidean")

# partial mantel test for 2012 gc w_4n samples controlling for climate
partial_mantel_gc12_w4n_spatial <- mantel.partial(gc12_w4n_dist, spatial_gc12_w4n_dist, clim_gc12_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc12_w4n_spatial) # significant

# partial mantel test for 2012 gc w4n samples controlling for spatial
partial_mantel_gc12_w4n_clim <- mantel.partial(gc12_w4n_dist, clim_gc12_w4n_dist, spatial_gc12_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc12_w4n_clim) # Mantel statistic r: -0.08712; not significant

###########################################
## Partial Mantel for spatial and climate data for 2021 GC - W_4n ####
# match row names for mantel test
# 2021 GC
gc21_w4n_geo_samples <- intersect(rownames(OCG_GC_2021_w4n_subset), rownames(site_coords)) # 26 samples

# Subset and reorder both
gc21_w4n_geo <- OCG_GC_2021_w4n_subset[gc21_w4n_geo_samples, ]
# spatial_gc21_w4n <- spatial_df[gc21_w4n_geo_samples, ]

spatial_df_GC21_w4n <- subset(site_coords, row.names(site_coords) %in% row.names(OCG_GC_2021_w4n_subset))
spatial_gc21_w4n_dist <- vegdist(spatial_df_GC21_w4n, method = "euclidean")

gc21_w4n_dist <- vegdist(gc21_w4n_geo, method = "euclidean")

# subset metadata to just w4n 2021 samples
md.clim_GC21_w4n <- subset(md_clim, row.names(md_clim) %in% row.names(OCG_GC_2021_w4n_subset))
clim_gc21_w4n <- md.clim_GC21_w4n[, which(colnames(md.clim_GC21_w4n) == "elev"):which(colnames(md.clim_GC21_w4n) == "mapmtcm")]

row.names(clim_gc21_w4n) == row.names(gc21_w4n_geo) # 50 samples
row.names(clim_gc21_w4n) == row.names(spatial_df_GC21_w4n)

clim_gc21_w4n_dist <- vegdist(clim_gc21_w4n, method = "euclidean")

# partial mantel test for 2021 gc w_4n samples controlling for climate
partial_mantel_gc21_w4n_spatial <- mantel.partial(gc21_w4n_dist, spatial_gc21_w4n_dist, clim_gc21_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc21_w4n_spatial) # not sig

# partial mantel test for 2021 gc w4n samples controlling for spatial
partial_mantel_gc21_w4n_clim <- mantel.partial(gc21_w4n_dist, clim_gc21_w4n_dist, spatial_gc21_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_gc21_w4n_clim) # not significant

###########################################
## RDA for spatial and climate data for LC-MS - T_2n ####
md.clim_LCMS_2012 <- subset(md_clim, md_clim$Year=="2012") # 153
md.clim_LCMS_t2n_2012 <- subset(md.clim_LCMS_2012, md.clim_LCMS_2012$Subsp_ploidy=="T_2n") # 54
OCG_LCMS_t2n_2012 <- subset(OCG_LCMS_subset, row.names(OCG_LCMS_subset) %in% row.names(md.clim_LCMS_t2n_2012)) # 30
md.clim_LCMS_t2n_2012 <- subset(md.clim_LCMS_t2n_2012, row.names(md.clim_LCMS_t2n_2012) %in% row.names(OCG_LCMS_t2n_2012)) # 30
OCG_LCMS_t2n_2012[is.na(OCG_LCMS_t2n_2012)] <- 0
OCG_LCMS_t2n_2012 <- OCG_LCMS_t2n_2012[order(row.names(OCG_LCMS_t2n_2012)),]
rownames(OCG_LCMS_t2n_2012) == rownames(md.clim_LCMS_t2n_2012)

# create geographic distance matrix for lcms T_2n samples
# library(geodist)
# geo_dist_lcms_t2n <- geodist(md.clim_LCMS_t2n[, c("Longitude", "Latitude")], measure = "geodesic")

# rda_lcms_clim_spatial<- vegan::rda(OCG_LCMS_t2n_2012 ~ PC1 + PC2 + PC3 + geo_dist_lcms_t2n, data = md.clim_LCMS_t2n)
# anova(rda_lcms_clim_spatial, by = "term")

# using the three climate variables from climate as predictors
# rda_lcms_clim_2012 <- vegan::rda(OCG_LCMS_t2n_2012 ~ long + mtwm + mtcm, data = md.clim_LCMS_t2n_2012)
# anova(rda_lcms_clim_2012, by = "term")
# 
# plot(rda_lcms_clim_2012, type = "n")
# points(rda_lcms_clim_2012, display = "sites", pch = 19)
# ordisurf(rda_lcms_clim_2012, md.clim_LCMS_t2n_2012$mtwm, 
#          col = "darkgreen", main = "mtwm", add = TRUE)
# 
# fit <- envfit(rda_lcms_clim_spatial2, 
#               md.clim_LCMS_t2n[, c("long", "mtwm", "mtcm")], 
#               permutations = 999)
# 
# # Plot significant vectors (e.g., p < 0.05)
# plot(fit, p.max = 0.05, col = "red")


# vegan::vif.cca(rda_gc_clim_spatial)

# scores_gc_clim_spatial <- scores(rda_gc_clim_spatial, display = "sites")
# scores_gc_clim_spatial <- as.data.frame(scores_gc_clim_spatial)
# 
# ordiplot(rda_gc_clim_spatial, type = "t",display = "sites",cex = .6)
# 
# plot(scores_gc_clim_spatial[, 1], scores_gc_clim_spatial[, 2],
#      xlab="axis 1", ylab="axis 2", 
#      col= c("pink","brown",'darkgreen','tan','lightblue')[md.clim_GC_t2n$Subsp_ploidy],
#      pch=c(17,19)[md.clim_GC_t2n$Year])
# legend("topleft", 
#        legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
#        col= c("pink","brown",'darkgreen','tan','lightblue'),
#        pch=19,
#        cex=0.8,
#        bty = "n")

## Partial Mantel for spatial and climate data for 2012 LC-MS - T_2n ####
# match row names for mantel test
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

# # create geographic distance matrix from the first two PCs in a spatial PCA
# # library(geodist)
# geo_dist_gc12_t2n <- geodist(site_coords_gc12_t2n, measure = "geodesic")
# geo_dist_gc12_t2n <- as.dist(geo_dist_gc12_t2n)

lcms12_t2n_dist <- vegdist(lcms12_t2n_geo, method = "euclidean")

# subset climate metadata to just t_2n 2012 samples
md.clim_lcms12_t2n <- subset(md_clim, row.names(md_clim) %in% row.names(lcms12_t2n_geo))
clim_lcms12_t2n <- md.clim_lcms12_t2n[, which(colnames(md.clim_lcms12_t2n) == "elev"):which(colnames(md.clim_lcms12_t2n) == "mapmtcm")]

row.names(clim_lcms12_t2n) == row.names(spatial_df_lcms12_t2n) 
row.names(clim_lcms12_t2n) == row.names(lcms12_t2n_geo)

clim_lcms12_t2n_dist <- vegdist(clim_lcms12_t2n, method = "euclidean")

# partial mantel test for 2012 lcms T_2n samples controlling for climate
partial_mantel_lcms12_t2n_spatial <- mantel.partial(lcms12_t2n_dist, spatial_df_lcms12_t2n_dist, clim_lcms12_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_lcms12_t2n_spatial) #significant

# partial mantel test for 2012 lcms T_2n samples controlling for spatial
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
md.clim_lcms21_t2n <- subset(md_clim, row.names(md_clim) %in% row.names(lcms21_t2n_geo))
clim_lcms21_t2n <- md.clim_lcms21_t2n[, which(colnames(md.clim_lcms21_t2n) == "elev"):which(colnames(md.clim_lcms21_t2n) == "mapmtcm")]

row.names(clim_lcms21_t2n) == row.names(spatial_df_lcms21_t2n) 
row.names(clim_lcms21_t2n) == row.names(lcms21_t2n_geo)

clim_lcms21_t2n_dist <- vegdist(clim_lcms21_t2n, method = "euclidean")

# partial mantel test for 2021 lcms T_2n samples controlling for climate
partial_mantel_lcms21_t2n_spatial <- mantel.partial(lcms21_t2n_dist, spatial_df_lcms21_t2n_dist, clim_lcms21_t2n_dist, method = "spearman", permutations = 999)
print(partial_mantel_lcms21_t2n_spatial) # sig

# partial mantel test for 2021 lcms T_2n samples controlling for spatial
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
md.clim_lcms21_w4n <- subset(md_clim, row.names(md_clim) %in% row.names(lcms21_w4n_geo))
clim_lcms21_w4n <- md.clim_lcms21_w4n[, which(colnames(md.clim_lcms21_w4n) == "elev"):which(colnames(md.clim_lcms21_w4n) == "mapmtcm")]

row.names(clim_lcms21_w4n) == row.names(spatial_df_lcms21_w4n) 
row.names(clim_lcms21_w4n) == row.names(lcms21_w4n_geo)

clim_lcms21_w4n_dist <- vegdist(clim_lcms21_w4n, method = "euclidean")

# partial mantel test for 2021 lcms w_4n samples controlling for climate
partial_mantel_lcms21_w4n_spatial <- mantel.partial(lcms21_w4n_dist, spatial_df_lcms21_w4n_dist, clim_lcms21_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_lcms21_w4n_spatial) 

# partial mantel test for 2021 w4n samples controlling for spatial
partial_mantel_lcms21_w4n_clim <- mantel.partial(lcms21_w4n_dist, clim_lcms21_w4n_dist, spatial_df_lcms21_w4n_dist, method = "spearman", permutations = 999)
print(partial_mantel_lcms21_w4n_clim) 

###########################################
## PERMANOVA for climate variables for LC-MS 2012 T_2n ####
row.names(OCG_LCMS_2012_t2n_subset) == row.names(md.clim_lcms12_t2n) 
permanova_lcms_2012_t2n <- adonis2(OCG_LCMS_2012_t2n_subset ~ tdiff + mtwm + mtcm, data = md.clim_lcms12_t2n, permutations = 999, by = "margin")
print(permanova_lcms_2012_t2n) 

# PCA on T_2n 2012 LCMS samples for ordisurf plots ####
OCG_LCMS_2012_t2n_subset_scaled <- scale(OCG_LCMS_2012_t2n_subset) 
OCG_LCMS_2012_t2n_subset_scaled[is.na(OCG_LCMS_2012_t2n_subset_scaled)] <- 0
pca_lcms_2012_t2n <- prcomp(OCG_LCMS_2012_t2n_subset_scaled)
summary(pca_lcms_2012_t2n) # 19.70 9.30
rownames(OCG_LCMS_2012_t2n_subset_scaled) == rownames(md.clim_LCMS_t2n_2012) # TRUE

#drop levels for ecoregion
md.clim_LCMS_t2n_2012$Ecoregion <- droplevels(md.clim_LCMS_t2n_2012$Ecoregion)

# palette for ecoregions
lcms_eco_t_2n_2012_palette <- c("#E79069","#E9AD6D",'#EAC87F','#EAE29C','#8BC982', '#1CA890')


# Longitude
ordi_lcms_12_t_2n_long <- ordisurf(pca_lcms_2012_t2n ~ md.clim_LCMS_t2n_2012$tdiff, 
                              col = "black", main = "TDIFF", add = TRUE)
summary(ordi_lcms_12_t_2n_long)

plot(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
     type = "n",  # don't draw points yet
     xlab = "PC 1 (19.7%)", ylab = "PC 2 (9.29%)",
     xlim = range(pca_lcms_2012_t2n$x[, 1], na.rm = TRUE),
     ylim = range(pca_lcms_2012_t2n$x[, 2], na.rm = TRUE))
ordisurf(pca_lcms_2012_t2n$x[, 1:2], md.clim_LCMS_t2n_2012$tdiff, 
         col = "black", add = TRUE)
points(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
       col = lcms_eco_t_2n_2012_palette[md.clim_LCMS_t2n_2012$Ecoregion],
       pch = 19)
legend("topright",
       legend=c("Central Basin & Range","Colorado Plateaus","Columbia Plateau","Northern Basin & Range","Snake River Plain", "Wasatch & Uinta Mountains"),
       col= lcms_eco_t_2n_2012_palette,
       pch=19,
       cex=0.6,
       bty = "n")

# mtwm
ordi_lcms_12_t_2n_mtwm <- ordisurf(pca_lcms_2012_t2n, md.clim_LCMS_t2n_2012$mtwm, 
                                   col = "black", main = "mtwm", add = TRUE)
summary(ordi_lcms_12_t_2n_mtwm)

plot(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
     type = "n",  # don't draw points yet
     xlab = "PC 1 (19.7%)", ylab = "PC 2 (9.29%)",
     xlim = range(pca_lcms_2012_t2n$x[, 1], na.rm = TRUE),
     ylim = range(pca_lcms_2012_t2n$x[, 2], na.rm = TRUE))
ordisurf(pca_lcms_2012_t2n$x[, 1:2], md.clim_LCMS_t2n_2012$mtwm, 
         col = "black", add = TRUE)
points(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
       col = lcms_eco_t_2n_2012_palette[md.clim_LCMS_t2n_2012$Ecoregion],
       pch = 19)
legend("topright",
       legend=c("Central Basin & Range","Colorado Plateaus","Columbia Plateau","Northern Basin & Range","Snake River Plain", "Wasatch & Uinta Mountains"),
       col= lcms_eco_t_2n_2012_palette,
       pch=19,
       cex=0.6,
       bty = "n")
# mtcm
ordi_lcms_12_t_2n_mtcm <- ordisurf(pca_lcms_2012_t2n, md.clim_LCMS_t2n_2012$mtcm, 
                                   col = "black", main = "mtcm", add = TRUE)
summary(ordi_lcms_12_t_2n_mtcm)

plot(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
     type = "n",  # don't draw points yet
     xlab = "PC 1 (19.7%)", ylab = "PC 2 (9.29%)",
     xlim = range(pca_lcms_2012_t2n$x[, 1], na.rm = TRUE),
     ylim = range(pca_lcms_2012_t2n$x[, 2], na.rm = TRUE))
ordisurf(pca_lcms_2012_t2n$x[, 1:2], md.clim_LCMS_t2n_2012$mtcm, 
         col = "black", add = TRUE)
points(pca_lcms_2012_t2n$x[, 1], pca_lcms_2012_t2n$x[, 2],
       col = lcms_eco_t_2n_2012_palette[md.clim_LCMS_t2n_2012$Ecoregion],
       pch = 19)
legend("topright",
       legend=c("Central Basin & Range","Colorado Plateaus","Columbia Plateau","Northern Basin & Range","Snake River Plain", "Wasatch & Uinta Mountains"),
       col= lcms_eco_t_2n_2012_palette,
       pch=19,
       cex=0.6,
       bty = "n")


###########################################
## PERMANOVA for climate variables for LC-MS 2021 T_2n ####

rownames(OCG_LCMS_2021_t2n_subset) == rownames(md.clim_lcms21_t2n) # TRUE

sum(is.na(md.clim_lcms21_t2n$tdiff)) # 0

rda_lcms_2021_t2n <- rda(OCG_LCMS_2021_t2n_subset ~ tdiff + mtwm + mtcm, data = md.clim_lcms21_t2n)
vif.cca(rda_lcms_2021_t2n) 

permanova_lcms_2021_t2n <- adonis2(OCG_LCMS_2021_t2n_subset ~ tdiff + mtwm + mtcm, data = md.clim_lcms21_t2n, permutations = 999, by = "margin")
print(permanova_lcms_2021_t2n)

cor(md.clim_lcms21_t2n[, c("tdiff","mtwm","mtcm")], use = "pairwise.complete.obs")

# write and save this as a csv 
write.csv(md.clim_lcms21_t2n, "md_clim_lcms21_t2n.csv", row.names = TRUE)

# PCA on T_2n 2021 LCMS samples for ordisurf plot ####
md.OCG.LCMS_clim_2021 <- subset(md_clim_LCMS, md_clim_LCMS$Year=="2021") #71
md.OCG.LCMS_clim_2021_t_2n <- subset(md.OCG.LCMS_clim_2021, md.OCG.LCMS_clim_2021$Subsp_ploidy=="T_2n") #24

# subset the lcms data to match the metadata
OCG_LCMS_2021_t2n_subset <- subset(OCG_LCMS_t2n_subset, row.names(OCG_LCMS_t2n_subset) %in% row.names(md.OCG.LCMS_clim_2021_t_2n)) # 24

OCG_LCMS_2021_t2n_subset_scaled <- scale(OCG_LCMS_2021_t2n_subset) #24
OCG_LCMS_2021_t2n_subset_scaled[is.na(OCG_LCMS_2021_t2n_subset_scaled)] <- 0

# remove columns with zero standard deviation
col_sds <- apply(OCG_LCMS_2021_t2n_subset_scaled, 2, sd, na.rm = TRUE)
OCG_LCMS_2021_t2n_subset_filtered <- OCG_LCMS_2021_t2n_subset_scaled[, col_sds > 0]

pca_lcms_2021_t2n <- prcomp(OCG_LCMS_2021_t2n_subset_filtered)
summary(pca_lcms_2021_t2n) 

rownames(OCG_LCMS_2021_t2n_subset_filtered) == rownames(md.OCG.LCMS_clim_2021_t_2n) # TRUE

#drop levels for ecoregion
md.OCG.LCMS_clim_2021_t_2n$Ecoregion <- droplevels(md.OCG.LCMS_clim_2021_t_2n$Ecoregion)

# palette for ecoregions
lcms_eco_t_2n_2021_palette <- c("#E79069",'#EAC87F','#EAE29C','#8BC982', '#1CA890', '#089392')


# MTWM
ordi_lcms_21_t_2n <- ordisurf(pca_lcms_2021_t2n, md.OCG.LCMS_clim_2021_t_2n$mtwm, 
         col = "black", add = TRUE)
summary(ordi_lcms_21_t_2n)

plot(pca_lcms_2021_t2n$x[, 1], pca_lcms_2021_t2n$x[, 2],
     type = "n",  # don't draw points yet
     xlab = "PC 1 (20.89%)", ylab = "PC 2 (10.85%)",
     xlim = range(pca_lcms_2021_t2n$x[, 1], na.rm = TRUE),
     ylim = range(pca_lcms_2021_t2n$x[, 2], na.rm = TRUE))
ordisurf(pca_lcms_2021_t2n$x[, 1:2], md.OCG.LCMS_clim_2021_t_2n$mtwm, 
         col = "black", add = TRUE)
points(pca_lcms_2021_t2n$x[, 1], pca_lcms_2021_t2n$x[, 2],
       col = lcms_eco_t_2n_2021_palette[md.OCG.LCMS_clim_2021_t_2n$Ecoregion],
       pch = 19)
legend("bottomleft",
       legend = c("Central Basin & Range", "Columbia Plateau", "Northern Basin & Range","Snake River Plain", "Wasatch & Uinta Mountains", "Wyoming Basin"),
       col = lcms_eco_t_2n_2021_palette,
       pch = 19,
       cex = 0.6,
       bty = "n")




###########################################
# ANCOM BC for LCMS tridentata 2n Ecoregion####
table(md.OCG.LCMS.tri_2n$Ecoregion)
md.OCG.LCMS.tri.abc <- md.OCG.LCMS.tri_2n[!md.OCG.LCMS.tri_2n$Ecoregion %in% c("Wyoming Basin", "Colorado Plateaus"), ] # ecoregions with too small of sample sizes
md.OCG.LCMS.tri.abc$Ecoregion <- droplevels(md.OCG.LCMS.tri.abc$Ecoregion)
levels(md.OCG.LCMS.tri.abc$Ecoregion) # 5 ecoregions
OCG_LCMS_tri.abc <- subset(OCG_LCMS_t2n, row.names(OCG_LCMS_t2n) %in% row.names(md.OCG.LCMS.tri.abc)) 

rownames(OCG_LCMS_tri.abc) == rownames(md.OCG.LCMS.tri.abc) # TRUE
rounded_matrix.lcmstri <- as.matrix(OCG_LCMS_tri.abc)
rounded_matrix.lcmstri <- round(rounded_matrix.lcmstri)
rounded_matrix.lcmstri<- rounded_matrix.lcmstri
rounded_matrix.lcmstri[is.na(rounded_matrix.lcmstri)] <- 0
rounded_matrix.lcmstri <- as.data.frame(rounded_matrix.lcmstri)


# Create the tse object
assays.lcmstri = S4Vectors::SimpleList(counts = t(rounded_matrix.lcmstri)) 
smd.lcmstri = S4Vectors::DataFrame(md.OCG.LCMS.tri.abc)
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

# save model
saveRDS(output.tri, "ancombc phytochemistry models/ancombc_lcms_loc_fdr.tri.RDS")
output.tri <- readRDS("ancombc phytochemistry models/ancombc_lcms_loc_fdr.tri.RDS")

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

abundance_sig.tri <- abundance_data_tri[significant_compounds.tri, , drop = FALSE] # subset abundane data to match significant compounds 

# log transform abundances 
abundance_log.tri <- log1p(abundance_sig.tri)  # log(1 + abundance)

# heat map with metadata call in 
annotation_col.tri <- data.frame(Ecoregion = md.OCG.LCMS.tri.abc$Ecoregion,
                                 row.names = rownames(md.OCG.LCMS.tri.abc))
rownames(annotation_col.tri) == colnames(abundance_log.tri)  # TRUE

Ecoregion_cluster <- hclust(dist(as.numeric(md.OCG.LCMS.tri.abc$Ecoregion)))

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
# log_abundance_data.tri <- log_abundance_data.tri[,1:length(md.OCG.LCMS.tri.abc$Location)] ## stops working here 
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

## Figure for ANCOMBC LOCATION for tridentata ####
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
# ANCOM BC 2012 GC subspecies ploidy ####
md.OCG.GC.2012.abc <- md.OCG.GC.2012
OCG_GC_2012.abc <- OCG_GC_2012

#adding subspecies ploidy number to row names
rownames(md.OCG.GC.2012.abc) <- paste(rownames(md.OCG.GC.2012.abc), md.OCG.GC.2012.abc$Subsp_ploidy, sep = "_") #sep argument controls the separator
rownames(OCG_GC_2012.abc) <- rownames(md.OCG.GC.2012.abc)

rownames(OCG_GC_2012.abc) == rownames(md.OCG.GC.2012.abc) # TRUE
# OCG_LCMS_tri.abc <- OCG_LCMS_tri.abc[,-1]
rounded_matrix.gc12 <- as.matrix(OCG_GC_2012.abc)
rounded_matrix.gc12 <- round(rounded_matrix.gc12)
rounded_matrix.gc12<- rounded_matrix.gc12
rounded_matrix.gc12[is.na(rounded_matrix.gc12)] <- 0
rounded_matrix.gc12 <- as.data.frame(rounded_matrix.gc12)

# Create the tse object
assays.gc12 = S4Vectors::SimpleList(counts = t(rounded_matrix.gc12)) 
smd.gc12 = S4Vectors::DataFrame(md.OCG.GC.2012.abc)
tse.gc12 = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays.gc12, colData = smd.gc12)

# run model 
output.gc12 = ancombc2(data = tse.gc12, assay_name = "counts", tax_level = NULL,
                       fix_formula = "Subsp_ploidy", rand_formula = NULL,
                       p_adj_method = "fdr", pseudo_sens = TRUE,
                       prv_cut = 0.20, lib_cut = 1000, s0_perc = 0.05,
                       group = "Subsp_ploidy", struc_zero = FALSE, neg_lb = FALSE,
                       alpha = 0.01, n_cl = 2, verbose = TRUE,
                       global = TRUE, pairwise = TRUE, 
                       dunnet = FALSE, trend = FALSE,
                       iter_control = list(tol = 1e-5, max_iter = 20, 
                                           verbose = FALSE),
                       em_control = list(tol = 1e-5, max_iter = 100),
                       lme_control = NULL, 
                       mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100), 
                       trend_control = NULL)

# save model
saveRDS(output.gc12, "ancombc phytochemistry models/ancombc_gc12_subsppl_fdr.RDS")
output.gc12 <- readRDS("ancombc phytochemistry models/ancombc_gc12_subsppl_fdr.RDS")

#subres#subpair#subet to the significant compounds between locations #43
res.gc12 <- output.gc12$res

# Rename columns containing "(Intercept)" to "Subsp_ploidyT_2n"
colnames(res.gc12) <- gsub("\\(Intercept\\)", "Subsp_ploidyT_2n", colnames(res.gc12))
pair.gc12 <- output.gc12$res_pair
dim(pair.gc12) #47 61

all_pair_sig.gc12 <- subset(res.gc12, diff_Subsp_ploidyT_2n == TRUE |diff_Subsp_ploidyT_4n == TRUE | diff_Subsp_ploidyV_2n == TRUE | diff_Subsp_ploidyV_4n == TRUE | diff_Subsp_ploidyW_4n == TRUE)

## Figure for ANCOMBC 2012 gc across subspecies and ploidy ####
lfc_long.gc12 <- all_pair_sig.gc12 %>%
  dplyr::select(taxon, starts_with("lfc_")) %>%  # Select compound and LFC columns
  tidyr::pivot_longer(cols = starts_with("lfc_"), 
                      names_to = "comparison", 
                      values_to = "lfc") %>%
  dplyr::mutate(comparison = gsub("lfc_", "", comparison))  %>%
  dplyr::mutate(comparison = factor(comparison, levels = c ("Subsp_ploidyT_2n", "Subsp_ploidyT_4n", "Subsp_ploidyV_2n", "Subsp_ploidyV_4n", "Subsp_ploidyW_4n")))# Keep rows without underscores

lfc_long.gc12 <- lfc_long.gc12 %>% 
  dplyr::mutate(lfc = ifelse(is.na(lfc), 0, lfc))

lfc_long.gc12$comparison <- sub("^Subsp_ploidy", "", lfc_long.gc12$comparison) #removing location from the start of all location names 

# Pivot the data to a matrix format
heatmap_data.gc12 <- reshape2::dcast(lfc_long.gc12, taxon ~ comparison, value.var = "lfc")
rownames(heatmap_data.gc12) <- heatmap_data.gc12$taxon
heatmap_matrix.gc12 <- as.matrix(heatmap_data.gc12[, -1])  # Remove the `taxon` column

pheatmap(heatmap_matrix.gc12,
         color = colorRampPalette(c("white", "mediumorchid2", "mediumorchid4"))(100),  # Custom color scale
         clustering_distance_rows = "euclidean",  # Distance metric for rows
         clustering_distance_cols = "euclidean",  # Distance metric for columns
         clustering_method = "complete",         # Clustering method
         scale = "none",                        # No scaling of data
         angle = 0)
####

significant_compounds.gc12 <- lfc_long.gc12 %>%
  dplyr::select(taxon,comparison) 

significant_abundance_data.gc12 <- OCG_GC_2012.abc[,colnames(OCG_GC_2012.abc) %in% significant_compounds.gc12$taxon]

significant_abundance_data.gc12 <- as.data.frame(significant_abundance_data.gc12)

significant_abundance_data.gc12[is.na(significant_abundance_data.gc12)] <- 0

abundance_log.gc12 <- log1p(significant_abundance_data.gc12)  # log(1 + abundance)

mean_abundance_gc12_df <- abundance_log.gc12 %>%
  # Extract subspecies ploidy
  mutate(Subsp_ploidy = str_sub(rownames(.), -4)) %>%
  # Group by subspecies ploidy
  group_by(Subsp_ploidy) %>%
  # Average log abundances for each compound
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  # Set subspecies ploidy as row names (optional, for pheatmap)
  column_to_rownames("Subsp_ploidy")

mean_abundance_gc12_df.t <- t(mean_abundance_gc12_df)

pheatmap(mean_abundance_gc12_df.t,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "mediumorchid2", "mediumorchid4"))(100), # Consider scaling by "column" or "none" as needed
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "2012 GC", 
         angle_col = 0)

# # Plot histogram of mean abundances
# hist(log10(mean_abundance_gc12_df.t + 1e-6),  # Adding a small constant to avoid log(0)
#      breaks = 50,
#      main = "Histogram of Mean Taxa Abundances (log10 scale)",
#      xlab = "Log10(Mean Relative Abundance + 1e-6)",
#      col = "skyblue", border = "black")

###########################################
# ANCOM BC 2021 GC subspecies ploidy ####
md.OCG.GC.2021.abc <- md.OCG.GC.2021[!md.OCG.GC.2021$Subsp_ploidy %in% c("V_2n", "V_4n"), ]
md.OCG.GC.2021.abc$Subsp_ploidy <- droplevels(as.factor(md.OCG.GC.2021.abc$Subsp_ploidy))
levels(md.OCG.GC.2021.abc$Subsp_ploidy)
OCG_GC_2021.abc <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(md.OCG.GC.2021.abc)) # 69 obs

#adding subspecies ploidy number to row names
rownames(md.OCG.GC.2021.abc) <- paste(rownames(md.OCG.GC.2021.abc), md.OCG.GC.2021.abc$Subsp_ploidy, sep = "_") #sep argument controls the separator

rownames(OCG_GC_2021.abc) <- rownames(md.OCG.GC.2021.abc) # add the ploidy argument to the end of the row names for ocg gc data 
rownames(OCG_GC_2021.abc) == rownames(md.OCG.GC.2021.abc) # TRUE

# OCG_LCMS_tri.abc <- OCG_LCMS_tri.abc[,-1]
rounded_matrix.gc21 <- as.matrix(OCG_GC_2021.abc)
rounded_matrix.gc21 <- round(rounded_matrix.gc21)
rounded_matrix.gc21<- rounded_matrix.gc21
rounded_matrix.gc21[is.na(rounded_matrix.gc21)] <- 0
rounded_matrix.gc21 <- as.data.frame(rounded_matrix.gc21)

# Create the tse object
assays.gc21 = S4Vectors::SimpleList(counts = t(rounded_matrix.gc21)) 
smd.gc21 = S4Vectors::DataFrame(md.OCG.GC.2021.abc)
tse.gc21 = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays.gc21, colData = smd.gc21)

# run model 
output.gc21 = ancombc2(data = tse.gc21, assay_name = "counts", tax_level = NULL,
                       fix_formula = "Subsp_ploidy", rand_formula = NULL,
                       p_adj_method = "fdr", pseudo_sens = TRUE,
                       prv_cut = 0.20, lib_cut = 1000, s0_perc = 0.05,
                       group = "Subsp_ploidy", struc_zero = FALSE, neg_lb = FALSE,
                       alpha = 0.05, n_cl = 2, verbose = TRUE,
                       global = TRUE, pairwise = TRUE, 
                       dunnet = FALSE, trend = FALSE,
                       iter_control = list(tol = 1e-5, max_iter = 20, 
                                           verbose = FALSE),
                       em_control = list(tol = 1e-5, max_iter = 100),
                       lme_control = NULL, 
                       mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100), 
                       trend_control = NULL)

# save model
saveRDS(output.gc21, "ancombc phytochemistry models/ancombc_gc21_subsppl_fdr.RDS")
output.gc21 <- readRDS("ancombc phytochemistry models/ancombc_gc21_subsppl_fdr.RDS")

res.gc21 <- output.gc21$res

# Rename columns containing "(Intercept)" to "Subsp_ploidyT_2n"
colnames(res.gc21) <- gsub("\\(Intercept\\)", "Subsp_ploidyT_2n", colnames(res.gc21))
pair.gc21 <- output.gc21$res_pair
dim(pair.gc21) 
#  19

all_pair_sig.gc21 <- subset(res.gc21, diff_Subsp_ploidyT_2n == TRUE |diff_Subsp_ploidyT_4n == TRUE | diff_Subsp_ploidyW_4n == TRUE )


## Figure for gc 21 across subspecies ploidy ####
lfc_long.gc21 <- all_pair_sig.gc21 %>%
  dplyr::select(taxon, starts_with("lfc_")) %>%  # Select compound and LFC columns
  tidyr::pivot_longer(cols = starts_with("lfc_"), 
                      names_to = "comparison", 
                      values_to = "lfc") %>%
  dplyr::mutate(comparison = gsub("lfc_", "", comparison))  %>%
  dplyr::mutate(comparison = factor(comparison, levels = c ("Subsp_ploidyT_2n", "Subsp_ploidyT_4n", "Subsp_ploidyW_4n")))# Keep rows without underscores

lfc_long.gc21 <- lfc_long.gc21 %>% 
  dplyr::mutate(lfc = ifelse(is.na(lfc), 0, lfc)) # any NAs to 0 

lfc_long.gc21$comparison <- sub("^Subsp_ploidy", "", lfc_long.gc21$comparison) #removing Subsp ploidy from the start of all names 

# Pivot the data to a matrix format
heatmap_data.gc21 <- reshape2::dcast(lfc_long.gc21, taxon ~ comparison, value.var = "lfc")
rownames(heatmap_data.gc21) <- heatmap_data.gc21$taxon
heatmap_matrix.gc21 <- as.matrix(heatmap_data.gc21[, -1])  # Remove the `taxon` column

pheatmap(heatmap_matrix.gc21,
         color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100),  # Custom color scale
         clustering_distance_rows = "euclidean",  # Distance metric for rows
         clustering_distance_cols = "euclidean",  # Distance metric for columns
         clustering_method = "complete",         # Clustering method
         scale = "none",                        # No scaling of data
         angle = 0)
####

significant_compounds.gc21 <- lfc_long.gc21 %>%
  dplyr::select(taxon,comparison) 

significant_abundance_data.gc21 <- OCG_GC_2021.abc[,colnames(OCG_GC_2021.abc) %in% significant_compounds.gc21$taxon]

significant_abundance_data.gc21 <- as.data.frame(significant_abundance_data.gc21)
significant_abundance_data.gc21[is.na(significant_abundance_data.gc21)] <- 0

abundance_log.gc21 <- log1p(significant_abundance_data.gc21)  # log(1 + abundance)

mean_abundance_gc21_df <- abundance_log.gc21 %>%
  # Extract subspecies ploidy
  mutate(Subsp_ploidy = str_sub(rownames(.), -4)) %>%
  # Group by subspecies ploidy
  group_by(Subsp_ploidy) %>%
  # Average log abundances for each compound
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  # Set subspecies ploidy as row names (optional, for pheatmap)
  column_to_rownames("Subsp_ploidy")

mean_abundance_gc21_df.t <- t(mean_abundance_gc21_df)

pheatmap(mean_abundance_gc21_df.t,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "mediumorchid2", "mediumorchid4"))(100), # Consider scaling by "column" or "none" as needed
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "2021 GC", 
         angle_col = 0)

# # Plot histogram of mean abundances
# hist(log10(mean_abundance_gc21_df.t + 1e-6),  # Adding a small constant to avoid log(0)
#      breaks = 50,
#      main = "Histogram of Mean Taxa Abundances (log10 scale)",
#      xlab = "Log10(Mean Relative Abundance + 1e-6)",
#      col = "skyblue", border = "black")
# 
# # Print summary statistics
# summary(log10(mean_abundance_gc21_df.t + 1e-6))

###########################################
# Procrustes for LCMS against GC data ####
# subset original OCG uncleaned data with original GC data 
OCG_GCp <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(OCG_LCMS_3uL)) #108 of 74 variables 
OCG_LCMS_3uLp <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(OCG_GC)) #108 of 307 variables 

rownames(OCG_GCp) == rownames(OCG_LCMS_3uLp) #TRUE

# Replace NA with 0
OCG_LCMS_3uLp[is.na(OCG_LCMS_3uLp)] <- 0
OCG_GCp[is.na(OCG_GCp)] <- 0

GC_v_LCMS <- protest(OCG_LCMS_3uLp, OCG_GCp, scale = TRUE) 
#technically this will run with LCMS being transfomed to GC 

# Calculate variance of each variable in OCG_LCMS_3uLp
variances <- apply(OCG_LCMS_3uLp, 2, var) 

# Select top N variables with highest variance (adjust N as needed)
n_top_vars <- ncol(OCG_GCp) # Example: Match the number of variables in OCG_GCp
top_vars <- names(sort(variances, decreasing = TRUE)[1:n_top_vars])

# Subset OCG_LCMS_3uLp to include only the top variables
OCG_LCMS_3uLp_reduced <- OCG_LCMS_3uLp[, top_vars]

GC_v_LCMS <- protest(OCG_LCMS_3uLp_reduced, OCG_GCp, scale = TRUE) # LCMS is being transformed to match with GC 
summary(GC_v_LCMS) 
plot(GC_v_LCMS)

# Create the data frame for plotting (Correctly using $X and $Yrot):
procrustes_df <- data.frame(
  x = GC_v_LCMS$X[, 1],       # LCMS data (the one being transformed) - Dimension 1
  y = GC_v_LCMS$X[, 2],       # LCMS data - Dimension 2
  xend = GC_v_LCMS$Yrot[, 1],  # GC data - Dimension 1
  yend = GC_v_LCMS$Yrot[, 2],  # GC data - Dimension 2
  Sample = rownames(GC_v_LCMS$X) # Sample names (use either dataset's rownames)
)

## procrustes plot #### 
ggplot(procrustes_df) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.6, color = "gray") +
  geom_point(aes(x = x, y = y, color = "LC-MS"), size = 2, shape = 19) +  # LCMS points
  geom_point(aes(x = xend, y = yend, color = "GC"), size = 2, shape = 17) + # Transformed GC points
  labs(x = "Procrustes axis 1", y = "Procrustes axis 2", color = "Data Source") +
  scale_color_manual(values = c("LC-MS" = "maroon", "GC" = "lightseagreen")) +
  theme_classic()+
  theme(legend.title = element_blank())

###########################################
# Procrustes using PCs ####
# using subsetted data that contains the same samples OCG_LCMSp and OCG_GCp 

# I have mutliple columns where compounds were never detected from the raw data so 
non_zero_vars <- names(OCG_GCp)[apply(OCG_GCp, 2, var) != 0]
OCG_GCp_filtered <- OCG_GCp[, non_zero_vars]

gc_pca <- prcomp(OCG_GCp_filtered, scale = TRUE)
lcms_pca <- prcomp(OCG_LCMS_3uLp, scale = TRUE)

GC_v_LCMS.pca <- protest(lcms_pca, gc_pca, scale = TRUE)
plot(GC_v_LCMS.pca)

procrustes_pca_df <- data.frame(
  x = GC_v_LCMS.pca$X[, 1],
  y = GC_v_LCMS.pca$X[, 2],
  xend = GC_v_LCMS.pca$Yrot[, 1],
  yend = GC_v_LCMS.pca$Yrot[, 2],
  Sample = rownames(lcms_pca$x) # Sample names
)

ggplot(procrustes_pca_df) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.6, color = "gray") +
  geom_point(aes(x = x, y = y, color = "LC-MS"), size = 2, shape = 19) +  # LCMS points
  geom_point(aes(x = xend, y = yend, color = "GC"), size = 2, shape = 17) + # Transformed GC points
  labs(x = "Procrustes axis 1", y = "Procrustes axis 2", color = "Data Source") +
  scale_color_manual(values = c("LC-MS" = "maroon", "GC" = "lightseagreen")) +
  theme_classic()+
  theme(legend.title = element_blank())


# # # # # # TRUE# # # # # # TRUETRUE
# rownames(abundance_sig.gc12) <- as.vector(abundance_sig.gc12[,1])
# sig_abundance_data.gc12 <- sig_abundance_data.gc12[,-1]
# 
# filtered_data <- abundance_data_gc12[sig_abundance_data.gc12, ]
# 
# # Rename columns containing "(Intercept)" to "AZ"
# colnames(res.gc12) <- gsub("\\(Intercept\\)", "Subsp_ploidyT_2n", colnames(res.gc12))
# pair.gc12 <- output.gc12$res_pair
# dim(pair.gc12) 
# #53 61
# 
# #subres#subpair#subet to the significant compounds between locations #43
# all_pair_sig.gc12 <- subset(res.gc12, diff_Subsp_ploidyT_2n == TRUE |diff_Subsp_ploidyT_4n == TRUE | diff_Subsp_ploidyV_2n == TRUE | diff_Subsp_ploidyV_4n == TRUE | diff_Subsp_ploidyW_4n == TRUE )
# 
# # Figure for GC 12 subspecies ploidy
# lfc_long.gc12 <- all_pair_sig.gc12 %>%
#   dplyr::select(taxon, starts_with("lfc_")) %>%  # Select compound and LFC columns
#   tidyr::pivot_longer(cols = starts_with("lfc_"), 
#                       names_to = "comparison", 
#                       values_to = "lfc")
#   # dplyr::mutate(comparison = gsub("lfc_", "", comparison)) %>%  # Clean comparison names
#   # dplyr::filter(!grepl("_", comparison)) %>%
#   # dplyr::mutate(comparison = factor(comparison, levels = c ("Subsp_ploidyT_2n", "Subsp_ploidyT_4n", "Subsp_ploidyV_2n", "Subsp_ploidyV_4n", "Subsp_ploidyW_4n")))# Keep rows without underscores
# 
# lfc_long.gc12 <- lfc_long.gc12 %>% 
#   dplyr::mutate(lfc = ifelse(is.na(lfc), 0, lfc))
# 
# lfc_long.gc12 <- lfc_long.gc12 %>%
#   mutate(comparison = recode(comparison,
#                              "lfc_Subsp_ploidyT_2n" = "T_2n",
#                              "lfc_Subsp_ploidyT_4n" = "T_4n",
#                              "lfc_Subsp_ploidyV_2n" = "V_2n",
#                              "lfc_Subsp_ploidyV_4n" = "V_4n",
#                              "lfc_Subsp_ploidyW_4n" = "W_4n"))
# 
# # Pivot the data to a matrix format
# heatmap_data <- reshape2::dcast(lfc_long.gc12, taxon ~ comparison, value.var = "lfc")
# rownames(heatmap_data) <- heatmap_data$taxon
# heatmap_matrix <- as.matrix(heatmap_data[, -1])  # Remove the `taxon` column
# 
# pheatmap(heatmap_matrix,
#          color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100),  # Custom color scale
#          clustering_distance_rows = "euclidean",  # Distance metric for rows
#          clustering_distance_cols = "euclidean",  # Distance metric for columns
#          clustering_method = "complete",         # Clustering method
#          scale = "none",                        # No scaling of data
#          angle = 0)

###########################################
#Creating a map of sites ####
# Read in shapefile
eco_shapefile <- st_read("data_csv/us_eco_l3_state_boundaries/us_eco_l3_state_boundaries.shp")

western_states <- c("Arizona", "New Mexico", "California", "Idaho", "Nevada", "Montana", "Utah", "Wyoming", "Colorado", "Oregon", "Washington")

eco_west <- eco_shapefile %>%
  filter(STATE_NAME %in% western_states)

# Turn site coordinates into an sf object
points_sf <- st_as_sf(md.OCG, coords = c("Longitude", "Latitude"), crs = 4326)

# Add the coordinates for the orchard
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

ggplot() +
  geom_sf(data = eco_west, aes(fill = highlight_group), color = "white", size = 0.3) +
  scale_fill_manual(values = final_colors, name = "Ecoregions") +
  geom_sf(data = points_sf, color = "darkslategrey", size = 1.5) +
  geom_sf(data = orchard_sf, shape = 7, color = "white", size = 3) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf()






















###########################################
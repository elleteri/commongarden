# Install and load necessary packages ####
if (!require("ANCOMBC")) {BiocManager::install("ANCOMBC"); require("ANCOMBC")} # ancom
if (!require("readr")) {BiocManager::install("readr"); require("readr")} # reading in data
if (!require("dplyr")) {BiocManager::install("dplyr"); require("dplyr")} #cleaning code
if (!require("tidyverse")) {BiocManager::install("tidyverse"); require("tidyverse")} #cleaning code
if (!require("ggplot2")) {BiocManager::install("ggplot2"); require("ggplot2")} #visualizing
if (!require("vegan")) {BiocManager::install("vegan"); require("vegan")} # community data analysis
if (!require("pairwiseAdonis")) {devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"); require("pairwiseAdonis")}

# Cleaned data read in####
#METADATA
md.OCG <- read.csv("data_csv/metadata_OCG.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #246 obs of 16 variables.
md.OCG <- md.OCG[order(row.names(md.OCG)),]

### Remove duplicates, negatve controls, and MTW.3.7.R_2012
rows_to_remove <- c('CAT.2.9_2012v1', 'CAV.2.7_2012v2','NVT.2.9_2012v2','ORT.2.10_2012v1','WAT.1.4_2012v2','WAT.1.9_2012v2','WAT.2.8_2012v1', 'ORT.1.5_2012', 'NEG_8-28-21', 'NEG_10-2-20', 'MTW.3.7.R_2012')
md.OCG <- md.OCG[!rownames(md.OCG) %in% rows_to_remove, ] #234 of 21 var

#make variables factor to plot and droplevels
md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Location","Description")] <- lapply(md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Location","Description")], as.factor)
md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Location","Description")] <- lapply(md.OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant","2020 STATUS","Location","Description")], droplevels)
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

# Alpha diversity ####
## GC alpha diversity ####
OCG.GC.shannon <- diversity(OCG_GC)
OCG.GC.ef <- exp(OCG.GC.shannon)
OCG.GC.ef.r <- round(OCG.GC.ef)

md.OCG.GC <- cbind(md.OCG.GC, effective_species = OCG.GC.ef.r)

glm.OCG.GC <- glm(effective_species ~ Year + Subsp_ploidy, family = Gamma(link = "log"), data = md.OCG.GC)
summary(glm.OCG.GC) #year and subspecies sig

plot(allEffects(glm.OCG.GC))

ggplot(data = md.OCG.GC, aes(Subsp_ploidy, effective_species, fill =Subsp_ploidy)) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "plasma")+
  labs(y = "Number of Compounds", x = "Subspecies + Ploidy")+
  theme_classic()+
  theme(legend.position = "none")

ggplot(md.OCG.GC, aes(Year, effective_species))+
  geom_boxplot(aes(group = Year, fill = Year))+
  scale_fill_viridis_d(option = "plasma")+
  labs(y = "Number of Compounds")+
  theme_classic()+
  theme(legend.position = "none")

## LCMS alpha diversity ####
OCG.LCMS.shannon <- diversity(OCG_LCMS_3uL)
OCG.LCMS.ef <- exp(OCG.LCMS.shannon)
OCG.LCMS.ef.r <- round(OCG.LCMS.ef)

md.OCG.LCMS.3 <- cbind(md.OCG.LCMS.3, Compounds = OCG.LCMS.ef.r)
glm.OCG.LCMS <- glm(compounds ~ Subsp_ploidy + Year, family = Gamma(link = "log"), data = md.OCG.LCMS.3)
summary(glm.OCG.LCMS)

plot(allEffects(glm.OCG.LCMS))

ggplot(md.OCG.LCMS.3, aes(Year, compounds))+
  geom_boxplot(aes(group = Year, fill = Year))+
  scale_fill_viridis_d(option = "plasma")+
  labs(y = "Number of Compounds")+
  theme_classic()+
  theme(legend.position = "none")

ggplot(data = md.OCG.LCMS.3, aes(Subsp_ploidy, compounds, fill = Subsp_ploidy)) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "plasma")+
  labs(y = "Number of Compounds", x = "Subspecies + Ploidy")+
  theme_classic()+
  theme(legend.position = "none")

ggplot(md.OCG.LCMS.3, aes(Ploidy, compounds, fill = Ploidy))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "plasma")+
  labs(y = "Number of Compounds")+
  theme_classic() +
  theme(legend.position = "none")


# Cleaning for PCAs####
## 2012 GC threshold defined ####
OCG_GC_2012[OCG_GC_2012 == 0] <- NA
colSums(is.na(OCG_GC_2012)) #checking for null values since they are used for filtering

#So now the data needs to be cleaned to only contain compounds that occur in more than 20% of the samples (plants). AKA columns need to be removed that don't have at least 20% of the rows containing a number=NA

#Calculate proportion of NA values that are in each column
na_proportion <- colMeans(is.na(OCG_GC_2012))
print(na_proportion)

#Now define the threshold of 10% - there are 73 compounds that remain after this
threshold <- 0.90

#identify which columns I need to keep
columns_to_keep <- na_proportion <= threshold

# Subset dataframe to keep only columns with NA proportion <= threshold
OCG_GC_2012_subset <- OCG_GC_2012[, columns_to_keep] #147 obs of 74 variables

# Replace NA values with zeroes
OCG_GC_2012_subset[is.na(OCG_GC_2012_subset)] <- 0
colSums(is.na(OCG_GC_2012_subset))

## 2021 GC threshold defined ####
OCG_GC_2021[OCG_GC_2021 == 0] <- NA
colSums(is.na(OCG_GC_2021))

na_proportion <- colMeans(is.na(OCG_GC_2021))
print(na_proportion)

threshold <- 0.90

columns_to_keep <- na_proportion <= threshold

OCG_GC_2021_subset <- OCG_GC_2021[, columns_to_keep] #70 obs of 43 var variables

OCG_GC_2021_subset[is.na(OCG_GC_2021_subset)] <- 0
colSums(is.na(OCG_GC_2021_subset)) 

## Full GC threshold defined #### 
OCG_GC[OCG_GC == 0] <- NA
colSums(is.na(OCG_GC)) 

na_proportion <- colMeans(is.na(OCG_GC))
print(na_proportion)

threshold <- 0.90

columns_to_keep <- na_proportion <= threshold

OCG_GC_subset <- OCG_GC[, columns_to_keep] #217 obs of 54 variables

OCG_GC_subset[is.na(OCG_GC_subset)] <- 0

colSums(is.na(OCG_GC_subset)) 

## LCMS threshold defined ####
OCG_LCMS_3uL[OCG_LCMS_3uL == 0] <- NA
colSums(is.na(OCG_LCMS_3uL))

na_proportion <- colMeans(is.na(OCG_LCMS_3uL))
print(na_proportion)

threshold <- 0.90

columns_to_keep <- na_proportion <= threshold

OCG_LCMS_3uL_subset <- OCG_LCMS_3uL[, columns_to_keep] #111 of 302 variables

OCG_LCMS_3uL_subset[is.na(OCG_LCMS_3uL_subset)] <- 0

colSums(is.na(OCG_LCMS_3uL_subset)) 

## 2012 GC scaling ####
#Compound
data_normalized_2012 <- scale(OCG_GC_2012_subset) #this will center the data 1:147, 1:47

#2012 plant ID
OCG_GC_2012_subset.t <- t(OCG_GC_2012_subset) #transpose the data first to put plant ID as columns and compound as row
colSums(is.na(OCG_GC_2012_subset.t)) #checking for null values
data_normalized_2012.ID <- scale(OCG_GC_2012_subset.t) #1:47, 1:147

## 2021 GC scaling ####
#Compound 
data_normalized_2021 <- scale(OCG_GC_2021_subset) #1:70, 1:42

#2021 plant ID
OCG_GC_2021_subset.t <- t(OCG_GC_2021_subset) 
colSums(is.na(OCG_GC_2021_subset.t))
data_normalized_2021.ID <- scale(OCG_GC_2021_subset.t) #1:42, 1:70

## Full GC scaling ####
#Compound
data_normalized_GC <- scale(OCG_GC_subset) #1:217, 1:54

#Plant ID
OCG_GC_subset.t <- t(OCG_GC_subset)
colSums(is.na(OCG_GC_subset.t))
data_normalized_GC.ID <- scale(OCG_GC_subset.t) #1:54, 1:217

## LCMS scaling ####
#Compound
data_LCMS3_normalized <- scale(OCG_LCMS_3uL_subset) #1:111, 1:302

#Plant ID
OCG_LCMS_3uL.t <- t(OCG_LCMS_3uL_subset)
colSums(is.na(OCG_LCMS_3uL.t)) 
data_LCMS3_normalized.ID <- scale(OCG_LCMS_3uL.t) #1:302, 1:111 



## Correlation ####
### 2012 GC CORRELATION
corr_matrix_2012 <- cor(data_normalized_2012) #47 compounds

### 2021 GC CORRELATION 
corr_matrix_2021 <- cor(data_normalized_2021) #42 compounds

### GC CORRELATION
corr_matrix_GC <- cor(data_normalized_GC) #53 compounds

### LCMS CORRELATION 
#Compound
corr_matrix_LCMS3 <- cor(data_LCMS3_normalized)#large matrix 91204 elements


# PCAs ####
set.seed(2425)
## 2012 GC PCA ####
#prcomp() has improved numerical accuracy, so is preferable to use this function.
data.pca_2012_ID <- prcomp(data_normalized_2012)
summary(data.pca_2012_ID)
autoplot(data.pca_2012_ID, label = TRUE)

rownames(data_normalized_2012) == rownames(md.OCG.GC.2012)

## PCA plot for 2012 GC subspecies ploidy ####
plot(data.pca_2012_ID$x[, 1], data.pca_2012_ID$x[, 2],
     xlab="PC 1 (19.81%)", ylab="PC 2 (12.46%)", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC.2012$Subsp_ploidy],
     pch=c(17),
     xlim = range(data.pca_2012_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2012_ID$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=17,
       cex=0.8,
       bty = "n")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(data.pca_2012_ID,groups = md.OCG.GC.2012$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

### PERMANOVA for 2012 GC subspecies ploidy ####
pca_scores_GC12 <- data.pca_2012_ID$x[]
GC_12_PCA_df <- as.data.frame(pca_scores_GC12)
pca_perm_GC_12_subsppl <- adonis2(pca_scores_GC12 ~ md.OCG.GC.2012$Subsp_ploidy, data = GC_12_PCA_df, method = "euclidean", by = "margin")
pca_perm_GC_12_subsppl # p = 0.001, R2 = 0.34638

OCG_GC_2012_subsp.pw.r <- pairwise.adonis(pca_scores_GC12, md.OCG.GC.2012$Subsp_ploidy, sim.method = "euclidean")
OCG_GC_2012_subsp.pw.r # sig between all subspecies

## 2021 GC PCA ####
data.pca_2021_ID <- prcomp(data_normalized_2021)
summary(data.pca_2021_ID)
autoplot(data.pca_2021_ID, label = TRUE)

rownames(data_normalized_2021) == rownames(md.OCG.GC.2021)

## PCA plot 2021 GC subspecies ploidy ####
plot(data.pca_2021_ID$x[, 1], data.pca_2021_ID$x[, 2],
     xlab="PC 1 (20.54%)", ylab="PC 2 (13.33%)",
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC.2021$Subsp_ploidy],
     pch=c(19),
     xlim = range(data.pca_2021_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_2021_ID$x[, 2], na.rm = TRUE))
legend("bottomright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(data.pca_2021_ID,groups = md.OCG.GC.2021$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")

### PERMANOVA for 2021 GC subspecies ploidy and location ####
pca_scores_GC21 <- data.pca_2021_ID$x[]
GC_21_PCA_df <- as.data.frame(pca_scores_GC21)
pca_perm_GC_21_subsppl <- adonis2(pca_scores_GC21 ~ md.OCG.GC.2021$Subsp_ploidy, data = GC_21_PCA_df, method = "euclidean", by = "margin")
pca_perm_GC_21_subsppl # p = 0.002, R2 = 0.23099

OCG_GC_2021_subsp.pw.r <- pairwise.adonis(pca_scores_GC21, md.OCG.GC.2021$Subsp_ploidy, sim.method = "euclidean")
OCG_GC_2021_subsp.pw.r # sig between all subspecies

## Full GC PCA ####
data.pca_GC_ID <- prcomp(data_normalized_GC)
summary(data.pca_GC_ID)
autoplot(data.pca_GC_ID, label = TRUE)

rownames(data_normalized_GC) == rownames(md.OCG.GC)

## PCA plot of full GC subspecies ploidy and year ####
plot(data.pca_GC_ID$x[, 1], data.pca_GC_ID$x[, 2],
     xlab="PC 1 (17.6%)", ylab="PC 2 (10.3%)", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.GC$Subsp_ploidy],
     pch=c(17,19)[md.OCG.GC$Year],
     xlim = range(data.pca_GC_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_GC_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown","darkgreen",'tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
legend("topleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")

### PERMANOVA for GC subspecies ploidy and location ####
pca_scores_GC <- data.pca_GC_ID$x[]
GC_PCA_df <- as.data.frame(pca_scores_GC)
pca_perm_GC_subsppl <- adonis2(pca_scores_GC ~ md.OCG.GC$Subsp_ploidy + md.OCG.GC$Year + md.OCG.GC$Location, data = GC_PCA_df, method = "euclidean", by = "margin")
pca_perm_GC_subsppl # p = 0.001, R2 = 0.15791 for subbspl R2 = 0.51987, p = 0.001 for year

## LCMS PCA ####
data.pca_LCMS3_ID <- prcomp(data_LCMS3_normalized)
summary(data.pca_LCMS3_ID)
autoplot(data.pca_LCMS3_ID, label = TRUE)

rownames(data_LCMS3_normalized) == rownames(md.OCG.LCMS.3)

## PCA plots of LCMS subspecies ploidy and year ####
plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1 (13.51%)", ylab="PC 2 (9.51%)", 
     col= c("pink","brown",'darkgreen','tan','lightblue')[md.OCG.LCMS.3$Subsp_ploidy],
     pch=c(17,19)[md.OCG.LCMS.3$Year],
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("pink","brown",'darkgreen','tan','lightblue'),
       pch=19,
       cex=0.8,
       bty = "n")
legend("topleft", 
       legend=c("2012","2021"),
       col= "black",
       pch=c(17,19),
       cex=0.8,
       bty = "n")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "T_2n", col = "pink")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "T_4n", col = "brown")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "V_2n", col = "darkgreen")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "V_4n", col = "tan")
ordispider(data.pca_LCMS3_ID,groups = md.OCG.LCMS.3$Subsp_ploidy, show.groups = "W_4n", col = "lightblue")


#"lightcoral","goldenrod",'olivedrab','lightseagreen','thistle'

## PCA plot for LCMS location and subspecies ####
plotcolor <- c("firebrick","cadetblue","sienna","skyblue","rosybrown","tomato","olivedrab","turquoise","burlywood","mediumaquamarine","darkseagreen")

plot(data.pca_LCMS3_ID$x[, 1], data.pca_LCMS3_ID$x[, 2],
     xlab="PC 1 (13.51%)", ylab="PC 2 (9.51%)", 
     col= plotcolor[md.OCG.LCMS.3$Location],
     pch=c(17,19,15)[md.OCG.LCMS.3$Subspecies],
     xlim = range(data.pca_LCMS3_ID$x[, 1], na.rm = TRUE),
     ylim = range(data.pca_LCMS3_ID$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("AZ","CA", "CO", "ID", "MT", "NM", "NV", "OR", "UT", "WA", "WY"),
       col= plotcolor,
       pch=19,
       cex=0.8,
       bty = "n")
legend("topleft", 
       legend=c("T","V","W"),
       col= "black",
       pch=c(17,19,15),
       cex=0.8,
       bty = "n")

### PERMANOVA for LCMS subspecies, and year ####
pca_scores_LCMS <- data.pca_LCMS3_ID$x[,1:2]
LCMS_PCA_df <- as.data.frame(pca_scores_LCMS)

pca_perm_lcms_subs <- adonis2(pca_scores_LCMS ~ md.OCG.LCMS.3$Subsp_ploidy + md.OCG.LCMS.3$Year, data = LCMS_PCA_df, by = "margin", method = "euclidean")
pca_perm_lcms_subs

lcms_subsp.pw.r <- pairwise.adonis(pca_scores_LCMS, md.OCG.LCMS.3$Subsp_ploidy, sim.method = "euclidean")
lcms_subsp.pw.r # sig between all subspecies

### PERMANOVA for LCMS Location, subspecies, and year ####
pca_perm_lcms_loc <- adonis2(pca_scores_LCMS ~ md.OCG.LCMS.3$Location + md.OCG.LCMS.3$Subsp_ploidy + md.OCG.LCMS.3$Year, data = LCMS_PCA_df, by = "margin", method = "euclidean")
pca_perm_lcms_loc 

##LCMS PCA of just tridentata ####
#subset LCMS to just tridentata
OCG_LCMS_tri <- subset(OCG_LCMS_3uL, md.OCG.LCMS.3$Subspecies=="T") #80 of 308 variables
md.OCG.LCMS.tri <- subset(md.OCG.LCMS.3, row.names(md.OCG.LCMS.3) %in% row.names(OCG_LCMS_tri)) #80

#make 0 NA to redefine threshold
OCG_LCMS_tri[OCG_LCMS_tri == 0] <- NA
colSums(is.na(OCG_LCMS_tri))
na_proportion <- colMeans(is.na(OCG_LCMS_tri))
print(na_proportion)
threshold <- 0.90
columns_to_keep <- na_proportion <= threshold
OCG_LCMS_tri_subset <- OCG_LCMS_tri[, columns_to_keep] #80 0f 297 var
OCG_LCMS_tri_subset[is.na(OCG_LCMS_tri_subset)] <- 0
colSums(is.na(OCG_LCMS_tri_subset)) 

#SCALING
#Compound
data_LCMS_tri_normalized <- scale(OCG_LCMS_tri_subset) #1:80, 1:297

#Plant ID
OCG_LCMS_tri.t <- t(OCG_LCMS_tri_subset)
colSums(is.na(OCG_LCMS_tri.t)) 
data_LCMS_tri_normalized.ID <- scale(OCG_LCMS_tri.t) #1:302, 1:111

#PCA
pca_LCMS.tri <- prcomp(data_LCMS_tri_normalized)
summary(pca_LCMS.tri)
autoplot(pca_LCMS.tri, label = TRUE)

rownames(data_LCMS_tri_normalized) == rownames(md.OCG.LCMS.tri)

# PCA plot for tridentata by location####
levels(md.OCG.LCMS.tri$Location)
plotcolor.tri <- c("firebrick","cadetblue","rosybrown","skyblue", "darkkhaki","coral","burlywood","goldenrod", "forestgreen")

plot(pca_LCMS.tri$x[, 1], pca_LCMS.tri$x[, 2],
     xlab="PC 1 (17.1%)", ylab="PC 2 (9.3%)", 
     col= plotcolor.tri[md.OCG.LCMS.tri$Location],
     pch=c(19),
     xlim = range(pca_LCMS.tri$x[, 1], na.rm = TRUE),
     ylim = range(pca_LCMS.tri$x[, 2], na.rm = TRUE))
legend("topleft", 
       legend=c("AZ","CA", "ID", "MT", "NM", "NV", "OR", "UT", "WA"),
       col= plotcolor.tri,
       pch=19,
       cex=0.8,
       bty = "n")

rownames(data_LCMS_tri_normalized) == rownames(md.OCG.LCMS.tri)

## PCA latitude plot ####
lat_colors <- viridis(n = 100)
latitude_scaled <- scales::rescale(md.OCG.LCMS.tri$Latitude, to = c(1, 100))  # Scale Latitude to match color range
color_by_lat <- lat_colors[latitude_scaled]  # Assign colors based on scaled Lat

plot(pca_LCMS.tri$x[, 1], pca_LCMS.tri$x[, 2],
     xlab="PC 1 (17.06%)", ylab="PC 2 (9.3%)", 
     col=color_by_lat,
     pch=c(19)[md.OCG.LCMS.tri$Subspecies],
     xlim=range(pca_LCMS.tri$x[, 1], na.rm=TRUE),
     ylim=range(pca_LCMS.tri$x[, 2], na.rm=TRUE))
image.plot(legend.only = TRUE, 
           zlim = range(md.OCG.LCMS.tri$Latitude, na.rm = TRUE), 
           col = viridis(100),  # Match legend with plot color scale
           legend.args = list(text = "Latitude", side = 4, line = 1, cex = 0.8),
           smallplot = c(0.945, 0.965, 0.3, 0.8),
           legend.mar = 1)

## PERMANOVA for wyomingensis location ####
pca_scores_lcms.tri <- pca_LCMS.tri$x[]
LCMS_PCA_tri.df <- as.data.frame(pca_scores_lcms.tri)

pca_perm_lcms_loc.tri <- adonis2(pca_scores_lcms.tri ~ md.OCG.LCMS.tri$Location, data = LCMS_PCA_tri.df, method = "euclidean")

pca_perm_lcms_loc.tri

## PERMANOVA of tridentata lat, long, elev ####
pca_perm_lcms_lat.tri <- adonis2(pca_scores_lcms.tri ~ md.OCG.LCMS.tri$Latitude + md.OCG.LCMS.tri$Longitude + md.OCG.LCMS.tri$`Elevation (m)`, data = LCMS_PCA_tri.df, by = "margin", method = "euclidean")

pca_perm_lcms_lat.tri # p = 0.001, R2 = 0.59972 for location R2 = 0.02828, p = 0.002 for year

## LCMS PCA of just wyomingensis ####
#subset LCMS to just wyomingensis
OCG_LCMS_wy <- subset(OCG_LCMS_3uL, md.OCG.LCMS.3$Subspecies=="W") 
md.OCG.LCMS.wy <- subset(md.OCG.LCMS.3, row.names(md.OCG.LCMS.3) %in% row.names(OCG_LCMS_wy)) 

#make 0 NA to redefine threshold
OCG_LCMS_wy[OCG_LCMS_wy == 0] <- NA
colSums(is.na(OCG_LCMS_wy))
na_proportion <- colMeans(is.na(OCG_LCMS_wy))
print(na_proportion)
threshold <- 0.90
columns_to_keep <- na_proportion <= threshold
OCG_LCMS_wy_subset <- OCG_LCMS_wy[, columns_to_keep] 
OCG_LCMS_wy_subset[is.na(OCG_LCMS_wy_subset)] <- 0
colSums(is.na(OCG_LCMS_wy_subset)) 

#SCALING
#Compound
data_LCMS_wy_normalized <- scale(OCG_LCMS_wy_subset) #1:26, 1:299

#Plant ID
OCG_LCMS_wy.t <- t(OCG_LCMS_wy_subset)
colSums(is.na(OCG_LCMS_wy.t)) 
data_LCMS_wy_normalized.ID <- scale(OCG_LCMS_wy.t) #1:299, 1:26

#PCA
pca_LCMS.wy <- prcomp(data_LCMS_wy_normalized)
summary(pca_LCMS.wy)
autoplot(pca_LCMS.wy, label = TRUE)

rownames(data_LCMS_wy_normalized) == rownames(md.OCG.LCMS.wy)

levels(md.OCG.LCMS.wy$Location)
plotcolor.wy <- c("rosybrown","skyblue", "darkkhaki","mediumaquamarine","goldenrod", "forestgreen")

plot(pca_LCMS.wy$x[, 1], pca_LCMS.wy$x[, 2],
     xlab="PC 1 (17.1%)", ylab="PC 2 (9.3%)", 
     main="PCA of LCMS tridentata data by location", 
     col= plotcolor.wy[md.OCG.LCMS.wy$Location],
     pch=c(19),
     xlim = range(pca_LCMS.wy$x[, 1], na.rm = TRUE),
     ylim = range(pca_LCMS.wy$x[, 2], na.rm = TRUE))
legend("topright", 
       legend=c("CO", "ID", "MT", "OR", "UT", "WA"),
       col= plotcolor.wy,
       pch=19,
       cex=0.8,
       bty = "n")

## PERMANOVA for wyomingensis location ####
pca_scores_lcms.wy <- pca_LCMS.wy$x[]
LCMS_PCA_wy.df <- as.data.frame(pca_scores_lcms.wy)

pca_perm_lcms_loc.wy <- adonis2(pca_scores_lcms.wy ~ md.OCG.LCMS.wy$Location, data = LCMS_PCA_wy.df, method = "euclidean")

pca_perm_lcms_loc.wy # p = 0.001, R2 = 0.59972 for location R2 = 0.02828, p = 0.002 for year

## PERMANOVA of wyomingensis lat, long, elev ####
pca_perm_lcms_lat.wy <- adonis2(pca_scores_lcms.wy ~ md.OCG.LCMS.wy$Latitude + md.OCG.LCMS.wy$Longitude + md.OCG.LCMS.wy$`Elevation (m)`, data = LCMS_PCA_wy.df, by = "margin", method = "euclidean")

pca_perm_lcms_lat.wy # p = 0.001, R2 = 0.59972 for location R2 = 0.02828, p = 0.002 for year


## wy PCA with longitude ####
long_colors <- viridis(n = 100, option = "plasma")
long_scaled <- scales::rescale(md.OCG.LCMS.wy$Longitude, to = c(1, 100))
color_by_long <- long_colors[long_scaled]  

plot(pca_LCMS.wy$x[, 1], pca_LCMS.wy$x[, 2],
     xlab="PC 1 (19.14%)", ylab="PC 2 (15.41%)", 
     col=color_by_long,
     pch= 19,
     xlim=range(pca_LCMS.wy$x[, 1], na.rm=TRUE),
     ylim=range(pca_LCMS.wy$x[, 2], na.rm=TRUE))
image.plot(legend.only = TRUE, 
           zlim = range(md.OCG.LCMS.wy$Longitude, na.rm = TRUE), 
           col = viridis(100, option = "plasma"),
           smallplot = c(0.85, 0.9, 0.3, 0.8),
           legend.args = list(text = "Longitude", side = 4, line = 2.5, cex = 0.8))


# ANCOM BC for tridentata location####
md.OCG.LCMS.tri.abc <- md.OCG.LCMS.tri[!md.OCG.LCMS.tri$Location %in% c("MT", "NV"), ]
md.OCG.LCMS.tri.abc$Location <- droplevels(md.OCG.LCMS.tri.abc$Location)
levels(md.OCG.LCMS.tri.abc$Location)
OCG_LCMS_tri.abc <- subset(OCG_LCMS_tri, row.names(OCG_LCMS_tri) %in% row.names(md.OCG.LCMS.tri.abc)) 

rownames(OCG_LCMS_tri.abc) == rownames(md.OCG.LCMS.tri.abc) # TRUE
# OCG_LCMS_tri.abc <- OCG_LCMS_tri.abc[,-1]
rounded_matrix <- as.matrix(OCG_LCMS_tri.abc)
rounded_matrix <- round(rounded_matrix)
rounded_matrix1<- rounded_matrix
rounded_matrix1[is.na(rounded_matrix1)] <- 0
rounded_matrix1 <- as.data.frame(rounded_matrix1)


# Create the tse object
assays = S4Vectors::SimpleList(counts = t(rounded_matrix1)) 
smd = S4Vectors::DataFrame(md.OCG.LCMS.tri.abc)
tse.tri = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)

# run model 
output.tri = ancombc2(data = tse.tri, assay_name = "counts", tax_level = NULL,
                  fix_formula = "Location", rand_formula = NULL,
                  p_adj_method = "fdr", pseudo_sens = TRUE,
                  prv_cut = 0.40, lib_cut = 1000, s0_perc = 0.05,
                  group = "Location", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.0001, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100), 
                  trend_control = NULL)

# save model
saveRDS(output.tri, "ancombc_lcms_loc_fdr.tri.RDS")
output.tri <- readRDS("ancombc_lcms_loc_fdr.tri.RDS")

# Extract raw abundance data (counts)
abundance_data_tri <- assay(tse.tri, "counts")
View(abundance_data_tri)

# Extract p-values from ANCOM-BC 2 results
p_vals.tri <- as.data.frame(output.tri$res_global$q_val)

abundance_sig.gc12 <- abundance_data_gc12[sig_compounds.gc12, , drop = FALSE]

# identify sig compounds = 22
significant_compounds.tri <- output.tri$res_global %>%
  dplyr::filter(q_val < 0.05) %>%
  dplyr::select(taxon) %>%
  dplyr::pull()

abundance_sig.tri <- abundance_data_tri[significant_compounds.tri, , drop = FALSE] # subset abundane data to match significant compounds 

# log transform abundances 
abundance_log.tri <- log1p(abundance_sig.tri)  # log(1 + abundance)

# heat map with metadata call in 
annotation_col.tri <- data.frame(Location = md.OCG.LCMS.tri.abc$Location)
rownames(annotation_col.tri) == colnames(abundance_log.tri)  # TRUE

Location_cluster <- hclust(dist(as.numeric(md.OCG.LCMS.tri.abc$Location)))

pheatmap(abundance_log.tri, 
         cluster_rows = TRUE, 
         cluster_cols = Location_cluster,
         show_colnames = TRUE,
         annotation_legend =  TRUE,  # Add metadata annotation
         annotation_names_col = FALSE,
         scale = "row",
         color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100),
         main = "Differentially Abundant Compounds - Tri location")


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
colnames(res.tri) <- gsub("\\(Intercept\\)", "LocationAZ", colnames(res.tri))
pair.tri <- output.tri$res_pair
dim(pair.tri) 
#238 127

str(res.tri)

#subres#subpair#subet to the significant compounds between locations #43
all_pair_sig.tri <- subset(res.tri, diff_LocationAZ | diff_LocationCA == TRUE |diff_LocationID == TRUE | diff_LocationNM == TRUE | diff_LocationOR == TRUE | diff_LocationUT == TRUE | diff_LocationWA == TRUE)

## Figure for ANCOMBC LOCATION for tridentata ####
lfc_long.tri <- all_pair_sig.tri %>%
  dplyr::select(taxon, starts_with("lfc_")) %>%  # Select compound and LFC columns
  tidyr::pivot_longer(cols = starts_with("lfc_"), 
                      names_to = "comparison", 
                      values_to = "lfc") %>%
  dplyr::mutate(comparison = gsub("lfc_", "", comparison)) %>%  # Clean comparison names
  dplyr::filter(!grepl("_", comparison)) %>%
  dplyr::mutate(comparison = factor(comparison, levels = c ("LocationAZ", "LocationCA", "LocationID","LocationNM", "LocationOR", "LocationUT", "LocationWA")))# Keep rows without underscores

lfc_long.tri <- lfc_long.tri %>% 
  dplyr::mutate(lfc = ifelse(is.na(lfc), 0, lfc))

lfc_long.tri$comparison <- sub("^Location", "", lfc_long.tri$comparison) #removing location from the start of all location names 
lfc_long.tri <- lfc_long.tri %>%
  mutate(comparison = recode(comparison,
                             "AZ" = "Arizona",
                             "CA" = "California",
                             "ID" = "Idaho",
                             "NM" = "New Mexico",
                             "OR" = "Oregon",
                             "UT" = "Utah",
                             "WA" = "Washington"))

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

state_abbrs <- substr(colnames(significant_abundance_data.tri), 1, 2)
significant_abundance_data.tri[is.na(significant_abundance_data.tri)] <- 0

abundance_log.tri <- log1p(significant_abundance_data.tri)  # log(1 + abundance)

# Step 2: Take the mean for columns grouped by state abbreviation
mean_abundance_tri <- abundance_log.tri %>%
  # Transpose and bind row names for easier manipulation
  t() %>%
  as.data.frame() %>%
  mutate(group = state_abbrs) %>%
  group_by(group) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "group")  # Set group names as row names

mean_abundance_tri_t <- t(mean_abundance_tri)

state_names <- c(AZ = "Arizona", 
                 CA = "California",
                 ID = "Idaho",
                 NM = "New Mexico",
                 OR = "Oregon", 
                 UT = "Utah",
                 WA = "Washington")

new_colnames <- colnames(mean_abundance_tri_t)
for (abbr in names(state_names)) {
  new_colnames <- gsub(paste0("^", abbr), state_names[abbr], new_colnames)
}
colnames(mean_abundance_tri_t) <- new_colnames

pheatmap(mean_abundance_tri_t,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100),  # Optional: Viridis color palette,  # Scale abundances by row
         show_rownames = TRUE,
         show_colnames = TRUE,
         angle_col = 0,
         fontsize_row = 8)

# ANCOM BC for wyomingensis location####
OCG_LCMS_wy <- subset(OCG_LCMS_3uL, md.OCG.LCMS.3$Subspecies=="W") 
md.OCG.LCMS.wy <- subset(md.OCG.LCMS.3, row.names(md.OCG.LCMS.3) %in% row.names(OCG_LCMS_wy)) 
md.OCG.LCMS.wy.abc <- md.OCG.LCMS.wy[!md.OCG.LCMS.wy$Location %in% c("AZ", "CA", "CO","MT", "NM","NV", "OR", "WA", "WY"), ]
md.OCG.LCMS.wy.abc$Location <- droplevels(md.OCG.LCMS.wy.abc$Location)
levels(md.OCG.LCMS.wy.abc$Location)
OCG_LCMS_wy.abc <- subset(OCG_LCMS_wy, row.names(OCG_LCMS_wy) %in% row.names(md.OCG.LCMS.wy.abc)) 

rownames(OCG_LCMS_wy.abc) == rownames(md.OCG.LCMS.wy.abc) # TRUE
# OCG_LCMS_tri.abc <- OCG_LCMS_tri.abc[,-1]
rounded_matrix.wy <- as.matrix(OCG_LCMS_wy.abc)
rounded_matrix.wy <- round(rounded_matrix.wy)
rounded_matrix.wy1<- rounded_matrix.wy
rounded_matrix.wy1[is.na(rounded_matrix.wy1)] <- 0
rounded_matrix.wy1 <- as.data.frame(rounded_matrix.wy1)


# Create the tse object
assays.wy = S4Vectors::SimpleList(counts = t(rounded_matrix.wy1)) 
smd.wy = S4Vectors::DataFrame(md.OCG.LCMS.wy.abc)
tse.wy = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays.wy, colData = smd.wy)

# run model 
output.wy = ancombc2(data = tse.wy, assay_name = "counts", tax_level = NULL,
                      fix_formula = "Location", rand_formula = NULL,
                      p_adj_method = "fdr", pseudo_sens = TRUE,
                      prv_cut = 0.40, lib_cut = 1000, s0_perc = 0.05,
                      group = "Location", struc_zero = FALSE, neg_lb = FALSE,
                      alpha = 0.0001, n_cl = 2, verbose = TRUE,
                      global = FALSE, pairwise = FALSE, 
                      dunnet = FALSE, trend = FALSE,
                      iter_control = list(tol = 1e-5, max_iter = 20, 
                                          verbose = FALSE),
                      em_control = list(tol = 1e-5, max_iter = 100),
                      lme_control = NULL, 
                      mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100), 
                      trend_control = NULL)


## ANCOMBC for LOCATION ####
output.tri = ancombc2(data = tse.tri, assay_name = "counts", tax_level = NULL,
                      fix_formula = "Location", rand_formula = NULL,
                      p_adj_method = "fdr", pseudo_sens = TRUE,
                      prv_cut = 0.40, lib_cut = 1000, s0_perc = 0.05,
                      group = "Location", struc_zero = FALSE, neg_lb = FALSE,
                      alpha = 0.0001, n_cl = 2, verbose = TRUE,
                      global = TRUE, pairwise = TRUE, 
                      dunnet = FALSE, trend = FALSE,
                      iter_control = list(tol = 1e-5, max_iter = 20, 
                                          verbose = FALSE),
                      em_control = list(tol = 1e-5, max_iter = 100),
                      lme_control = NULL, 
                      mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100), 
                      trend_control = NULL)

saveRDS(output.tri, "ancombc_lcms_loc_fdr.tri.RDS")
output.tri <- readRDS("ancombc_lcms_loc_fdr.tri.RDS")

res.tri <- output.tri$res

abundance_data <- output.tri$res_global

sig_abundance_data.tri <- abundance_data[output.tri$res_global$p_val < 0.05,]
sig_abundance_data.tri <- abundance_data[output.tri$res_global$diff_abn == "TRUE",]
rownames(sig_abundance_data.tri) <- as.vector(sig_abundance_data.tri[,1])
sig_abundance_data.tri <- sig_abundance_data.tri[,-1]
sig_abundance_data.tri <- sig_abundance_data.tri[,-4]

log_abundance_data.tri <- log(sig_abundance_data.tri + 1)

log_abundance_data.tri <- as.matrix(log_abundance_data.tri)  

log_abundance_data.tri <- log_abundance_data.tri[, 1:length(md.OCG.LCMS.tri$Location)] ## stops working here 

pheatmap(log_abundance_data.tri, cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = data.frame(Location = tse$Location),
         scale = "row", color = viridis::viridis(100))

# Rename columns containing "(Intercept)" to "AZ"
colnames(res.tri) <- gsub("\\(Intercept\\)", "LocationAZ", colnames(res.tri))
pair.tri <- output.tri$res_pair
dim(pair.tri) 
#238 127

str(res)

#subres#subpair#subet to the significant compounds between locations #43
all_pair_sig.tri <- subset(res.tri, diff_LocationCA == TRUE |diff_LocationCO == TRUE | diff_LocationID == TRUE | diff_LocationMT == TRUE | diff_LocationNM == TRUE | diff_LocationNV == TRUE | diff_LocationOR == TRUE | diff_LocationUT == TRUE | diff_LocationWA == TRUE)

## Figure for ANCOMBC LOCATION ####
lfc_long <- all_pair_sig %>%
  dplyr::select(taxon, starts_with("lfc_")) %>%  # Select compound and LFC columns
  tidyr::pivot_longer(cols = starts_with("lfc_"), 
                      names_to = "comparison", 
                      values_to = "lfc") %>%
  dplyr::mutate(comparison = gsub("lfc_", "", comparison)) %>%  # Clean comparison names
  dplyr::filter(!grepl("_", comparison)) %>%
  dplyr::mutate(comparison = factor(comparison, levels = c ("LocationAZ", "LocationNM", "LocationUT", "LocationCO", "LocationCA", "LocationNV", "LocationMT", "LocationID", "LocationOR", "LocationWA")))# Keep rows without underscores

lfc_long <- lfc_long %>% 
  dplyr::mutate(lfc = ifelse(is.na(lfc), 0, lfc))

lfc_long$comparison <- sub("^Location", "", lfc_long$comparison) #removing location from the start of all location names 
lfc_long <- lfc_long %>%
  mutate(comparison = recode(comparison,
                             "AZ" = "Arizona",
                             "CA" = "California",
                             "CO" = "Colorado",
                             "ID" = "Idaho",
                             "MT" = "Montana",
                             "NM" = "New Mexico",
                             "NV" = "Nevada",
                             "OR" = "Oregon",
                             "UT" = "Utah",
                             "WA" = "Washington"))

install.packages("pheatmap")
library(pheatmap)

# Pivot the data to a matrix format
heatmap_data <- reshape2::dcast(lfc_long, taxon ~ comparison, value.var = "lfc")
rownames(heatmap_data) <- heatmap_data$taxon
heatmap_matrix <- as.matrix(heatmap_data[, -1])  # Remove the `taxon` column

pheatmap(heatmap_matrix,
         color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100),  # Custom color scale
         clustering_distance_rows = "euclidean",  # Distance metric for rows
         clustering_distance_cols = "euclidean",  # Distance metric for columns
         clustering_method = "complete",         # Clustering method
         scale = "none",                        # No scaling of data
         angle = 0)
####

significant_compounds <- lfc_long %>%
  dplyr::select(taxon,comparison) 

significant_abundance_data <- OCG_LCMS_3uL[,colnames(OCG_LCMS_3uL) %in% significant_compounds$taxon]

significant_abundance_data <- as.data.frame(t(significant_abundance_data))

state_abbrs <- substr(colnames(significant_abundance_data), 1, 2)
significant_abundance_data[is.na(significant_abundance_data)] <- 0

# Step 2: Take the mean for columns grouped by state abbreviation
mean_abundance_df <- significant_abundance_data %>%
  # Transpose and bind row names for easier manipulation
  t() %>%
  as.data.frame() %>%
  mutate(group = state_abbrs) %>%
  group_by(group) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "group")  # Set group names as row names

mean_abundance_df_t <- t(mean_abundance_df)

state_names <- c(AZ = "Arizona", 
                 NM = "New Mexico", 
                 UT = "Utah", 
                 CO = "Colorado", 
                 CA = "California", 
                 NV = "Nevada", 
                 MT = "Montana", 
                 ID = "Idaho", 
                 OR = "Oregon", 
                 WA = "Washington")
new_colnames <- colnames(mean_abundance_df_t)
for (abbr in names(state_names)) {
  new_colnames <- gsub(paste0("^", abbr), state_names[abbr], new_colnames)
}
colnames(mean_abundance_df_t) <- new_colnames

pheatmap(mean_abundance_df_t,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100),  # Optional: Viridis color palette,
         scale = "row",  # Scale abundances by row
         show_rownames = TRUE,
         show_colnames = TRUE,
         angle_col = 0,
         fontsize_row = 8)

# how many plants died ####
md.OCG.2012$`2020 STATUS` <- factor(md.OCG.2012$`2020 STATUS`, levels = c("A", "D"))
levels(md.OCG.2012$`2020 STATUS`)


sum(md.OCG.2012$`2020 STATUS` == "D", na.rm = TRUE) # 75 had a 2020 dead status and there was 82 total that were lost between 2012 and 2021 





table(md.OCG.GC.2021$Subspecies)

# ANCOM BC 2012 GC subspecies ploidy ####
md.OCG.GC.2012.abc <- md.OCG.GC.2012
OCG_GC_2012.abc <- OCG_GC_2012
#adding subspecies ploidy number to row names
rownames(md.OCG.GC.2012.abc) <- paste(rownames(md.OCG.GC.2012.abc), md.OCG.GC.2012.abc$Subsp_ploidy, sep = "_") #sep argument controls the separator
levels(md.OCG.GC.2012.abc$Subsp_ploidy)
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
                      prv_cut = 0.06, lib_cut = 1000, s0_perc = 0.05,
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

#subres#subpair#subet to the significant compounds between locations #43
res.gc12 <- output.gc12$res

# Rename columns containing "(Intercept)" to "Subsp_ploidyT_2n"
colnames(res.gc12) <- gsub("\\(Intercept\\)", "Subsp_ploidyT_2n", colnames(res.gc12))
pair.gc12 <- output.gc12$res_pair
dim(pair.gc12) #53 61

all_pair_sig.gc12 <- subset(res.gc12, diff_Subsp_ploidyT_2n == TRUE |diff_Subsp_ploidyT_4n == TRUE | diff_Subsp_ploidyV_2n == TRUE | diff_Subsp_ploidyV_4n == TRUE | diff_Subsp_ploidyW_4n == TRUE)

## Figure for ANCOMBC LOCATION for tridentata ####
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
         color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100),  # Custom color scale
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
         color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100), # Consider scaling by "column" or "none" as needed
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "2012 GC", 
         angle_col = 0)

# Plot histogram of mean abundances
hist(log10(mean_abundance_gc12_df.t + 1e-6),  # Adding a small constant to avoid log(0)
     breaks = 50,
     main = "Histogram of Mean Taxa Abundances (log10 scale)",
     xlab = "Log10(Mean Relative Abundance + 1e-6)",
     col = "skyblue", border = "black")

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
                       prv_cut = 0.01, lib_cut = 1000, s0_perc = 0.05,
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

res.gc21 <- output.gc21$res

# Rename columns containing "(Intercept)" to "Subsp_ploidyT_2n"
colnames(res.gc21) <- gsub("\\(Intercept\\)", "Subsp_ploidyT_2n", colnames(res.gc21))
pair.gc21 <- output.gc21$res_pair
dim(pair.gc21) 
# 50 19

all_pair_sig.gc21 <- subset(res.gc21, diff_Subsp_ploidyT_2n == TRUE |diff_Subsp_ploidyT_4n == TRUE | diff_Subsp_ploidyW_4n == TRUE )


## Figure for ANCOMBC LOCATION for tridentata ####
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
         color = colorRampPalette(c("darkolivegreen", "white", "darkmagenta"))(100), # Consider scaling by "column" or "none" as needed
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "2021 GC", 
         angle_col = 0)

# Plot histogram of mean abundances
hist(log10(mean_abundance_gc21_df.t + 1e-6),  # Adding a small constant to avoid log(0)
     breaks = 50,
     main = "Histogram of Mean Taxa Abundances (log10 scale)",
     xlab = "Log10(Mean Relative Abundance + 1e-6)",
     col = "skyblue", border = "black")

# Print summary statistics
summary(log10(mean_abundance_gc21_df.t + 1e-6))


# Procrustes for LCMS against GC data ####
# subset original OCG uncleaned data with original GC data 
OCG_GCp <- subset(OCG_GC, row.names(OCG_GC) %in% row.names(OCG_LCMS_3uL)) #107 of 74 variables 
OCG_LCMS_3uLp <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(OCG_GC)) #107 of 307 variables 

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


# procrustes using pca matrices ####
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

# ####
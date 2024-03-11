#Read in cleaned data ####
#metadata read in
mdITS_OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #metadata read in #154 of 16 variables
head(mdITS_OCG)

#asv read in
asvITS.OCG <- read.csv("data_csv/asvITS.OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #asv table read in 154 obs of 2135 variables

summary(rowSums(asvITS.OCG)) #507
summary(colSums(asvITS.OCG)) #2

row.names(asvITS.OCG) == row.names(mdITS_OCG) # TRUE

#LCMS 1uL read in and subset with asv
OCG_LCMS_1uL <- read.csv("data_csv/OCG_LCMS_1uL_cleaned.csv", row.names = 1) #110 obs of 309 var
OCG_LCMS_1uL <- OCG_LCMS_1uL[order(row.names(OCG_LCMS_1uL)),] # order samples alphabetically

OCG_LCMS_1uL <- subset(OCG_LCMS_1uL, row.names(OCG_LCMS_1uL) %in% row.names(asvITS.OCG)) #58 of 309 variables

asvITS.OCG.LC1 <- subset(asvITS.OCG, row.names(asvITS.OCG) %in% row.names(OCG_LCMS_1uL)) #58 of 2135 variables

#LCMS 3uL read in and subset with asv
OCG_LCMS_3uL <- read.csv("data_csv/OCG_LCMS_3uL_cleaned.csv", row.names = 1) #110 obs of 309 var
OCG_LCMS_3uL <- OCG_LCMS_3uL[order(row.names(OCG_LCMS_3uL)),] # order samples alphabetically

OCG_LCMS_3uL <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(mdITS_OCG)) #59 of 309 variables

asvITS.OCG.LC3 <- subset(asvITS.OCG, row.names(asvITS.OCG) %in% row.names(OCG_LCMS_3uL)) #59 of 2135 variables

#2012 GC data read in and subset with asv
OCG_AUC_2012 <- read.csv("data_csv/OCG_AUC_2012_cleaned.csv", row.names = 1) 
#157 obs of 73 variables
OCG_AUC_2012 <- OCG_AUC_2012[order(row.names(OCG_AUC_2012)),] # order samples alphabetically

OCG_AUC_2012 <- subset(OCG_AUC_2012, row.names(OCG_AUC_2012) %in% row.names(asvITS.OCG)) #97 of 73 var
asvITS.OCG.GC12 <- subset(asvITS.OCG, row.names(asvITS.OCG) %in% row.names(OCG_AUC_2012)) #97 of 2135 var

#2021 GC data read in and subset with asv
OCG_AUC_2021 <- read.csv("data_csv/OCG_AUC_2021.csv", row.names = 1) 
#70 obs of 75 variables
OCG_AUC_2021 <- OCG_AUC_2021[order(row.names(OCG_AUC_2021)),] # order samples alphabetically

OCG_AUC_2021 <- subset(OCG_AUC_2021, row.names(OCG_AUC_2021) %in% row.names(asvITS.OCG)) #43 of 74 var
asvITS.OCG.GC21 <- subset(asvITS.OCG, row.names(asvITS.OCG) %in% row.names(OCG_AUC_2021)) #43 of 2135 var

##Read in nmds files
set.seed(41)
#rarefy
#LCMS 1
#summary(rowSums(asvITS.OCG.LC1)) #507
#summary(colSums(asvITS.OCG.LC1)) #0.00
#asvITS.LCMS_1.r <- rrarefy(asvITS.OCG.LC1,507) ## rarefy. warning message

#LCMS 3
#summary(rowSums(asvITS.OCG.LC3)) #507
#summary(colSums(asvITS.OCG.LC3)) #0.00
#asvITS.LCMS_3.r <- rrarefy(asvITS.OCG.LC3,507) ## rarefy. warning message

#2012 GC
#summary(rowSums(asvITS.OCG.GC12)) #507
#summary(colSums(asvITS.OCG.GC12)) #0.00
#asvITS.OCG.GC12.r <- rrarefy(asvITS.OCG.GC12,507) ## rarefy. warning message

#2021 GC
#summary(rowSums(asvITS.OCG.GC21)) #783
#summary(colSums(asvITS.OCG.GC21)) #0.00
#asvITS.OCG.GC21.r <- rrarefy(asvITS.OCG.GC21,783) ## rarefy. warning message

# Replace NA with 0
OCG_LCMS_1uL[is.na(OCG_LCMS_1uL)] <- 0
OCG_LCMS_3uL[is.na(OCG_LCMS_3uL)] <- 0
OCG_AUC_2012[is.na(OCG_AUC_2012)] <- 0
OCG_AUC_2021[is.na(OCG_AUC_2021)] <- 0

#asv nmds to match LCMS 1
#set.seed(7)
#asvITS_OCG_LCMS1.nmds <- metaMDS(asvITS.LCMS_1.r, trymax=500) ###solution reached! warning message
#save(asvITS_OCG_LCMS1.nmds, file = "nmds/asvITS_OCG_LCMS1.nmds.rda")

#LCMS 1 to match asv
#set.seed(9)
#OCG_LCMS1_pro.nmds <- metaMDS(OCG_LCMS_1uL, trymax=500) ###solution reached! 
#save(OCG_LCMS1_pro.nmds, file = "nmds/OCG_LCMS1_pro.nmds.rda")

#asv nmds to match LCMS 3
#set.seed(87)
#asvITS_OCG_LCMS3.nmds <- metaMDS(asvITS.LCMS_3.r, trymax=500) ###solution reached! warning message
#save(asvITS_OCG_LCMS3.nmds, file = "nmds/asvITS_OCG_LCMS3.nmds.rda")

#LCMS 3 to match asv
# set.seed(56)
# OCG_LCMS3_pro.nmds <- metaMDS(OCG_LCMS_3uL, trymax=500) ###solution reached!
# save(OCG_LCMS3_pro.nmds, file = "nmds/OCG_LCMS3_pro.nmds.rda")

#asv nmds to match 2012 GC
#set.seed(53)
#asvITS.OCG.GC12.nmds <- metaMDS(asvITS.OCG.GC12.r, trymax=500) ###solution reached! warning message
#save(asvITS.OCG.GC12.nmds, file = "nmds/asvITS.OCG.GC12.nmds.rda")

#2012 GC to match asv
#set.seed(85)
#OCG_GC12_pro.nmds <- metaMDS(OCG_AUC_2012, trymax=500) ###solution reached! 
#save(OCG_GC12_pro.nmds, file = "nmds/OCG_GC12_pro.nmds.rda")

#asv nmds to match 2021 GC
#set.seed(98)
#asvITS_OCG_GC21.nmds <- metaMDS(asvITS.OCG.GC21.r, trymax=500) ###solution reached! warning message
#save(asvITS_OCG_GC21.nmds, file = "nmds/asvITS_OCG_GC21.nmds.rda")

#2021 GC to match asv
#set.seed(5)
#OCG_GC21_pro.nmds <- metaMDS(OCG_AUC_2021, trymax=500) ###solution reached!
#save(OCG_GC21_pro.nmds, file = "nmds/OCG_GC21_pro.nmds.rda")

#load all nmds files
load("nmds/asvITS_OCG_LCMS1.nmds.rda") #asv file for LCMS 1uL data = asvITS.OCG
load("nmds/asvITS_OCG_LCMS3.nmds.rda") #asv file for LCMS 3uL data = asvITS.OCG
load("nmds/OCG_LCMS1_pro.nmds.rda") #LCMS 1uL file for asv data = OCG_LCMS_1uL
load("nmds/OCG_LCMS3_pro.nmds.rda") #LCMS 3uL file for asv data = OCG_LCMS_3uL

load("nmds/asvITS.OCG.GC12.nmds.rda") #asv file for GC 2012 data = asvITS.OCG
load("nmds/asvITS_OCG_GC21.nmds.rda") #asv file for GC 2021 data = asvITS.OCG
load("nmds/OCG_GC12_pro.nmds.rda") #2012 GC file for asv data = OCG_AUC_2012
load("nmds/OCG_GC21_pro.nmds.rda") #2021 GC file file for asv data = OCG_AUC_2021

#Procrustes analyses#
library(vegan)
library(ggplot2)

#LCMS 1uL####
asv_LCMS_1.pro <- protest(asvITS_OCG_LCMS1.nmds, OCG_LCMS1_pro.nmds, symmetric=T) 
asv_LCMS_1.pro ## Correlation in a symmetric Procrustes rotation: 0.1718, Significance: 0.323
summary(asv_LCMS_1.pro) 
plot(asv_LCMS_1.pro)
## Re-plotting the procrustes 
asv_LCMS_1_prodat <- as.data.frame(asv_LCMS_1.pro$X)
asv_LCMS_1_prodat <- cbind(asv_LCMS_1_prodat,asv_LCMS_1.pro$Yrot)
colnames(asv_LCMS_1_prodat)[colnames(asv_LCMS_1_prodat)=="1"] <- "Xend"
colnames(asv_LCMS_1_prodat)[colnames(asv_LCMS_1_prodat)=="2"] <- "Yend"

ggplot() + 
  geom_segment(data=asv_LCMS_1_prodat, mapping=aes(x=NMDS1, y=NMDS2, xend=Xend, yend=Yend), size=0.8, color="gray") + 
  geom_point(data=asv_LCMS_1_prodat, mapping=aes(x=NMDS1, y=NMDS2), size=2, shape=19, color = "maroon") +
  geom_point(data=asv_LCMS_1_prodat, mapping=aes(x=Xend, y=Yend), size=2, shape=17, color = "lightseagreen") +
  labs(x="Procrustes axis 1", y="Procrustes axis 2",title = "LCMS1 vs asv procrustes plot") +
  theme_classic() 


#LCMS 3uL####
asv_LCMS_3.pro <- protest(asvITS_OCG_LCMS3.nmds, OCG_LCMS3_pro.nmds, symmetric=T) 
asv_LCMS_3.pro ## Correlation in a symmetric Procrustes rotation: 0.1575, Significance: 0.434
summary(asv_LCMS_3.pro) 
plot(asv_LCMS_3.pro)

## Re-plotting the procrustes 
asv_LCMS_3_prodat <- as.data.frame(asv_LCMS_3.pro$X)
asv_LCMS_3_prodat <- cbind(asv_LCMS_3_prodat,asv_LCMS_3.pro$Yrot)
colnames(asv_LCMS_3_prodat)[colnames(asv_LCMS_3_prodat)=="1"] <- "Xend"
colnames(asv_LCMS_3_prodat)[colnames(asv_LCMS_3_prodat)=="2"] <- "Yend"

ggplot() + 
  geom_segment(data=asv_LCMS_3_prodat, mapping=aes(x=NMDS1, y=NMDS2, xend=Xend, yend=Yend), size=0.8, color="gray") + 
  geom_point(data=asv_LCMS_3_prodat, mapping=aes(x=NMDS1, y=NMDS2), size=2, shape=19, color = "maroon") +
  geom_point(data=asv_LCMS_3_prodat, mapping=aes(x=Xend, y=Yend), size=2, shape=17, color = "lightseagreen") +
  labs(x="Procrustes axis 1", y="Procrustes axis 2",title = "LCMS3 vs asv procrustes plot") +
  theme_classic() 

#2012 GC data####
asv_GC12.pro <- protest(asvITS.OCG.GC12.nmds, OCG_GC12_pro.nmds, symmetric=T) 
asv_GC12.pro ## Correlation in a symmetric Procrustes rotation: 0.1555, Significance: 0.184
summary(asv_GC12.pro) 
plot(asv_GC12.pro)

## Re-plotting the procrustes 
asv_GC12.pro_prodat <- as.data.frame(asv_GC12.pro$X)
asv_GC12.pro_prodat <- cbind(asv_GC12.pro_prodat,asv_GC12.pro$Yrot)
colnames(asv_GC12.pro_prodat)[colnames(asv_GC12.pro_prodat)=="1"] <- "Xend"
colnames(asv_GC12.pro_prodat)[colnames(asv_GC12.pro_prodat)=="2"] <- "Yend"

ggplot() + 
  geom_segment(data=asv_GC12.pro_prodat, mapping=aes(x=NMDS1, y=NMDS2, xend=Xend, yend=Yend), size=0.8, color="gray") + 
  geom_point(data=asv_GC12.pro_prodat, mapping=aes(x=NMDS1, y=NMDS2), size=2, shape=19, color = "maroon") +
  geom_point(data=asv_GC12.pro_prodat, mapping=aes(x=Xend, y=Yend), size=2, shape=17, color = "lightseagreen") +
  labs(x="Procrustes axis 1", y="Procrustes axis 2",title = "GC 2012 vs asv procrustes plot") +
  theme_classic() 

#2021 GC data
asv_GC21.pro <- protest(asvITS_OCG_GC21.nmds, OCG_GC21_pro.nmds, symmetric=T) 
asv_GC21.pro ## Correlation in a symmetric Procrustes rotation: 0.2263, Significance: 0.179
summary(asv_GC21.pro) 
plot(asv_GC21.pro)

## Re-plotting the procrustes 
asv_GC21.pro_prodat <- as.data.frame(asv_GC21.pro$X)
asv_GC21.pro_prodat <- cbind(asv_GC21.pro_prodat,asv_GC21.pro$Yrot)
colnames(asv_GC21.pro_prodat)[colnames(asv_GC21.pro_prodat)=="1"] <- "Xend"
colnames(asv_GC21.pro_prodat)[colnames(asv_GC21.pro_prodat)=="2"] <- "Yend"

ggplot() + 
  geom_segment(data=asv_GC21.pro_prodat, mapping=aes(x=NMDS1, y=NMDS2, xend=Xend, yend=Yend), size=0.8, color="gray") + 
  geom_point(data=asv_GC21.pro_prodat, mapping=aes(x=NMDS1, y=NMDS2), size=2, shape=19, color = "maroon") +
  geom_point(data=asv_GC21.pro_prodat, mapping=aes(x=Xend, y=Yend), size=2, shape=17, color = "lightseagreen") +
  labs(x="Procrustes axis 1", y="Procrustes axis 2",title = "GC 2021 vs asv procrustes plot") +
  theme_classic() 

#Mantel tests measuring correlation between the distance matrices##
#GC 2012
OCG_AUC_2012.dist <- vegdist(OCG_AUC_2012)
asvITS.OCG.GC12.dist <- vegdist(asvITS.OCG.GC12)
GC12.asv.mant <- mantel(OCG_AUC_2012.dist,asvITS.OCG.GC12.dist, permutations = 9999) 
GC12.asv.mant ## r = -0.02619, P =0.6681

#GC 2021
OCG_AUC_2021.dist <- vegdist(OCG_AUC_2021)
asvITS.OCG.GC21.dist <- vegdist(asvITS.OCG.GC21)
GC21.asv.mant <- mantel(OCG_AUC_2021.dist,asvITS.OCG.GC21.dist, permutations = 9999) 
GC21.asv.mant ## r = -0.01556, P =0.5612

#LCMS1
OCG_LCMS_1uL.dist <- vegdist(OCG_LCMS_1uL)
asvITS.OCG.LC1.dist <- vegdist(asvITS.OCG.LC1)
LCMS1.asv.mant <- mantel(OCG_LCMS_1uL.dist,asvITS.OCG.LC1.dist, permutations = 9999) 
LCMS1.asv.mant ## r = -0.1146, P =0.9645

#LCMS3
OCG_LCMS_3uL.dist <- vegdist(OCG_LCMS_3uL)
asvITS.OCG.LC3.dist <- vegdist(asvITS.OCG.LC3)
LCMS3.asv.mant <- mantel(OCG_LCMS_3uL.dist,asvITS.OCG.LC3.dist, permutations = 9999) 
LCMS3.asv.mant ## r = -0.1345, P =0.9879


#LCMS against 1ul and 3ul procrustes to check #### 
OCG_LCMS_3uL. <- subset(OCG_LCMS_3uL, row.names(OCG_LCMS_3uL) %in% row.names(OCG_LCMS_1uL)) #58
OCG_LCMS_1uL. <- subset(OCG_LCMS_1uL, row.names(OCG_LCMS_1uL) %in% row.names(OCG_LCMS_3uL)) #58

set.seed(72)
OCG_LCMS3_prot.nmds <- metaMDS(OCG_LCMS_3uL., trymax=500) ###solution reached!


set.seed(92)
OCG_LCMS1_prot.nmds <- metaMDS(OCG_LCMS_1uL., trymax=500) ###solution reached!

LCMS_3_1.pro <- protest(OCG_LCMS1_prot.nmds ,OCG_LCMS3_prot.nmds, symmetric=T) 
LCMS_3_1.pro ## Correlation in a symmetric Procrustes rotation: 0.7926, Significance: 0.001
summary(LCMS_3_1.pro) 
plot(LCMS_3_1.pro)


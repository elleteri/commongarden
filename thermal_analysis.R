# Stable isotope analysis#### 
#2012 
stable_iso_data_2012<-read.csv("data_csv/Stable_isotope_2012.csv")
str(stable_iso_data_2012) #96
# Extracting numbers from Plant_ID column
stable_iso_data_2012$Sample.ID <- as.numeric(gsub("[^0-9]", "", stable_iso_data_2012$ID))
stable_iso_data_2012$Sample.ID <- as.character(stable_iso_data_2012$Sample.ID)
stable_iso_data_2012$Subspecies <- recode(stable_iso_data_2012$Subspecies, 
                                          "Tridentata" = "T",
                                          "Vaseyana" = "V",
                                          "Wyomingensis" = "W")
levels(stable_iso_data_2012$Subspecies)
stable_iso_data_2012$Subspecies <- as.factor(stable_iso_data_2012$Subspecies)
str(stable_iso_data_2012)

#2021
stable_iso_data_2021<-read.csv("data_csv/Stable_isotope_2021.csv") #74
str(stable_iso_data_2021)
## Remove V-134677 (not sure what subspecies since it says "T/W"). look into later
stable_iso_data_2021 <- stable_iso_data_2021[!(row.names(stable_iso_data_2021) == "31"),] #73
# need to check what subspecies 134677 is 
levels(stable_iso_data_2021$Subspecies)<- list(Tridentata="T", Vaseyana="V", Wyomingensis="W")
levels(stable_iso_data_2021$Subspecies)
stable_iso_data_2021$Subspecies <- as.factor(stable_iso_data_2021$Subspecies)

#renaming the columns to combine the data together
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "Ampl..28"] <- "Ampl28"
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "Ampl..44"] <- "Ampl44"
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "d15N"] <- "Delta15N"
names(stable_iso_data_2021)[names(stable_iso_data_2021) == "d13C"] <- "Delta13C"

#FULL STABLE ISOTOPE DATA
stable_iso_df <- full_join(stable_iso_data_2012,stable_iso_data_2021) #169 of 15 var
stable_iso_df$Subspecies <- droplevels(stable_iso_df$Subspecies)
str(stable_iso_df)

# 2012 BETADISPERSION
# Delta 15 N
D15Ndist <- vegdist(stable_iso_data_2012$Delta15N, method = "bray")
D15Nbetadisper <- betadisper(D15Ndist, group = stable_iso_data_2012$Subspecies)
permutest(D15Nbetadisper) #0.343

D15N12fit<-aov(Delta15N~Subspecies,data = stable_iso_data_2012)
summary(D15N12fit) #sig between subspecies

D15Nsubsp.pw <- pairwise.adonis(D15Ndist, as.factor(stable_iso_data_2012$Subspecies))


#Delta 13 C
D13Cdist <- vegdist(stable_iso_data_2012$Delta13C, method = "euclidean")
D13Cbetadisper <- betadisper(D13Cdist, group = stable_iso_data_2012$Subspecies)
permutest(D13Cbetadisper) #0.036
#pairwiseadonis
D13Csubsp.pw.12 <- pairwise.adonis(D13Cdist,as.factor(stable_iso_data_2012$Subspecies))
D13Csubsp.pw.12 #T vs V sig

D13C12fit<-aov(Delta13C~Subspecies,data = stable_iso_data_2012)
summary(D13C12fit)

#2012 VISUALIZATION
# Grouped Scatter plot with marginal density plots
ggscatterhist(
  stable_iso_data_2012, x = "Delta13C", y = "Delta15N", group = "Subspecies",
  color = "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.params = list(fill = "Subspecies", color = "black", size = 0.2)
)

ggscatterhist(
  stable_iso_data_2012, x = "Delta13C", y = "Delta15N", group="Subspecies",
  color = "Subspecies", fill= "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.plot = "boxplot",
  margin.params = list(fill = "Subspecies", color = c("olivedrab", "cadetblue", "goldenrod"), size = 0.2),
  ggtheme = theme_bw()
)

#2021 BETADISPERSION
#DELTA 15 N
D15N21fit<-aov(Delta15N~Subspecies,data = stable_iso_data_2021)
summary(D15N21fit) #p -value= 0.1812

D15Ndist <- vegdist(stable_iso_data_2021$Delta15N, method = "euclidean")
D15Nbetadisper <- betadisper(D15Ndist, group = stable_iso_data_2021$Subspecies)
permutest(D15Nbetadisper) #0.79
#pairwiseadonis
D15Nsubsp.pw.21 <- pairwise.adonis(D15Ndist,as.factor(stable_iso_data_2021$Subspecies))
D15Nsubsp.pw.21 
#DELTA 13 C
D13Cdist <- vegdist(stable_iso_data_2021$Delta13C, method = "euclidean")
D13Cbetadisper <- betadisper(D13Cdist, group = stable_iso_data_2021$Subspecies)
permutest(D13Cbetadisper) #0.491
#pairwiseadonis
D13Csubsp.pw.21 <- pairwise.adonis(D13Cdist,as.factor(stable_iso_data_2021$Subspecies))
D13Csubsp.pw.21 

D13C21fit<-aov(Delta13C~Subspecies,data = stable_iso_data_2021)
summary(D13C21fit) #p-value = 0.037

D15Nsubsp.pw.21 <- pairwise.adonis(stable_iso_data_2021$Delta13C, stable_iso_data_2021$Subspecies)

#2021 VISUALIZATION
ggscatterhist(
  stable_iso_data_2021, "Delta13C", y = "Delta15N", group="Subspecies",
  color = "Subspecies", fill= "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.plot = "boxplot",
  margin.params = list(fill = "Subspecies", color = c("olivedrab", "cadetblue", "goldenrod"), size = 0.2),
  ggtheme = theme_bw()
)

#FULL VISUALIZATION
ggscatterhist(
  stable_iso_df, "Delta13C", y = "Delta15N", group="Subspecies",
  color = "Subspecies", fill= "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.plot = "boxplot",
  margin.params = list(fill = "Subspecies", color = c("olivedrab", "cadetblue", "goldenrod"), size = 0.2),
  ggtheme = theme_bw()
)

ggscatterhist(
  stable_iso_df, x = "Delta13C", y = "Delta15N", group = "Subspecies",
  color = "Subspecies", size = 3, alpha = 0.6,
  palette = c("olivedrab", "cadetblue", "goldenrod"),
  margin.params = list(fill = "Subspecies", color = "black", size = 0.2),
  ggtheme = theme_bw()
)

D15N_model <- adonis2(stable_iso_df$Delta15N ~ stable_iso_df$Subspecies) # not sig
D15Ndist <- vegdist(stable_iso_df$Delta15N)
D15Nbetadisper <- betadisper(D15Ndist, group = stable_iso_df$Subspecies)
permutest(D15Nbetadisper) #0.339
D15Nsubsp.pw <- pairwise.adonis(D15Ndist,as.factor(stable_iso_df$Subspecies))
D15Nsubsp.pw 

D13C_model <- adonis2(stable_iso_df$Delta13C ~ stable_iso_df$Subspecies, method = "euclidean") # not sig
D13Cdist <- vegdist(stable_iso_df$Delta13C, method = "euclidean")
D13Cbetadisper <- betadisper(D13Cdist, group = stable_iso_df$Subspecies)
permutest(D13Cbetadisper) #0.401
D13Csubsp.pw <- pairwise.adonis(D13Cdist,as.factor(stable_iso_df$Subspecies))
D13Csubsp.pw 

p1<-ggplot(stable_iso_data_2012,aes(Subspecies,Delta15N))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+ ggtitle("2012")+ theme(legend.position = "none") + ylab("Delta15N")
p2<-ggplot(stable_iso_data_2021,aes(Subspecies,Delta15N))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+ ggtitle("2021")+theme(legend.position = "none")
p3<-ggplot(stable_iso_data_2012,aes(Subspecies,Delta13C))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+ggtitle("2012")+ theme(legend.position = "none")
p4<-ggplot(stable_iso_data_2021,aes(Subspecies,Delta13C))+
  geom_boxplot(aes(fill=Subspecies))+theme_classic()+
  scale_fill_manual(values=c("olivedrab","cadetblue","goldenrod"))+ggtitle("2021")+theme(legend.position = "none")
grid.arrange(p1, p2, p3, p4, nrow = 2)

stable_iso_df$Year <- as.factor(stable_iso_df$Year)
D15N_year_fit <- aov(stable_iso_df$Delta15N ~ stable_iso_df$Year)
summary(D15N_year_fit)

D13C_year_fit <- aov(stable_iso_df$Delta13C ~ stable_iso_df$Year)
summary(D13C_year_fit)

stablefit <- lm(stable_iso_df$Delta15N ~ stable_iso_df$Delta13C)
summary(stablefit)
plot(stable_iso_df$Delta15N ~ stable_iso_df$Delta13C)

# mantel test between stable iso data and 2012 GC ####
md_gc_2012_cluster <- read.csv("data_csv/md.OCG.GC.2012_cluster.csv")
md_gc_2021_cluster <- read.csv("data_csv/md.OCG.GC.2021_cluster.csv")
old_name <- "Sample.ID"
new_name <- "Garden.Plant.ID"

# Rename the column
# 2012
names(stable_iso_data_2012)[names(stable_iso_data_2012) == old_name] <- new_name
md2012 <- merge(md_gc_2012_cluster, stable_iso_data_2012, by = "Garden.Plant.ID")

# 2021 
names(stable_iso_data_2021)[names(stable_iso_data_2021) == old_name] <- new_name
md2021 <- merge(md_gc_2021_cluster, stable_iso_data_2021, by = "Garden.Plant.ID")

md.w.iso <- merge(md2012, md2021, all = TRUE)
rownames(md.w.iso) <- md.w.iso[, 2]
md.w.iso <- md.w.iso[, -2]
md.w.iso <- md.w.iso[order(row.names(md.w.iso)),]

#GC with D15N
OCG_GC_mantel <- subset(OCG_GC_subset, row.names(OCG_GC_subset) %in% row.names(md.w.iso)) #140
OCG_GC_mantel <- OCG_GC_mantel[order(row.names(OCG_GC_mantel)),]
rownames(OCG_GC_mantel) == rownames(md.w.iso) #TRUE

OCG_GC_mantel.dist <- vegdist(OCG_GC_mantel)
D15dist <- vegdist(md.w.iso$Delta15N)

GC.D15N.mant <- mantel(OCG_GC_mantel.dist,D15dist, permutations = 9999) 
GC.D15N.mant ## r = 0.1926, sig < 0.001

# GC with D13C
D13dist <- vegdist(md.w.iso$Delta13C, method = "euclidean")

GC.D13C.mant <- mantel(OCG_GC_mantel.dist,D13dist, permutations = 9999) 
GC.D13C.mant ## r = 0.04033, sig = 0.0278

# LCMS 
OCG_LCMS_mantel <- subset(OCG_LCMS_3uL_subset, row.names(OCG_LCMS_3uL_subset) %in% row.names(md.w.iso))#65
md.w.iso.lcms <- subset(md.w.iso, row.names(md.w.iso) %in% row.names(OCG_LCMS_3uL_subset))
OCG_LCMS_mantel <- OCG_LCMS_mantel[order(row.names(OCG_LCMS_mantel)),]
rownames(OCG_LCMS_mantel) == rownames(md.w.iso.lcms) #TRUE

# D15N versus LCMS mantel test
OCG_LCMS_mantel.dist <- vegdist(OCG_LCMS_mantel)
D15dist.lcms <- vegdist(md.w.iso.lcms$Delta15N)

LCMS.D15N.mant <- mantel(OCG_LCMS_mantel.dist,D15dist.lcms, permutations = 9999) 
LCMS.D15N.mant ## -0.06531, P = 0.8764

# D13C
D13dist.lcms <- vegdist(md.w.iso.lcms$Delta13C, method = "euclidean")

LCMS.D13C.mant <- mantel(OCG_LCMS_mantel.dist,D13dist.lcms, permutations = 9999) 
LCMS.D13C.mant ## 0.1032, P = 0.0077

# D15N versus D13C 
D15N.D13C.mant <- mantel(D15dist,D13dist.lcms, permutations = 9999) 
D15N.D13C.mant ## 0.01514, P = 0.3566








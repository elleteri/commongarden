# GC plots
 
#housekeeping #### 
setwd("/Users/ellehorwath/Documents/Orchard_Common_Garden/commongarden")

library(ggplot2)
library(viridis)
library(ggthemes)
# detach("package:DBI", unload = TRUE)

#METADATA
mdITS <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) # 505 obs of 16 var
mdITS <- mdITS[order(row.names(mdITS)),] 
mdITS.OCG <- subset(mdITS, mdITS$Project == "OCG") #244 obs of 16 var   
mdITS.OCG <- subset(mdITS.OCG, mdITS.OCG$Location != "NEG")  

#GC
OCG_GC_w_RT <- read.csv("data_csv/OCG_GC_w_RT_full_clean.csv", head=T, check.names = F,stringsAsFactors = T) #226
rownames(OCG_GC_w_RT) <- OCG_GC_w_RT[,1]
mdITS.OCG.GC <- subset(mdITS.OCG, row.names(mdITS.OCG) %in% row.names(OCG_GC_w_RT)) #217
OCG_GC_w_RT <- OCG_GC_w_RT[order(row.names(OCG_GC_w_RT)),] 
OCG_GC_w_RT <- subset(OCG_GC_w_RT, row.names(OCG_GC_w_RT) %in% row.names(mdITS.OCG)) #217

row.names(OCG_GC_w_RT) == row.names(mdITS.OCG.GC) # sanity check: TRUE

OCG_GC_w_RT <- merge(mdITS.OCG.GC, OCG_GC_w_RT, by = "row.names")

#subset to just 2012 
OCG_GC_w_RT_2012 <- subset(OCG_GC_w_RT, OCG_GC_w_RT$Year=="2012") #147

#subset to just 2021
OCG_GC_w_RT_2021 <- subset(OCG_GC_w_RT, OCG_GC_w_RT$Year=="2021") #70

#subset to just tridentata
OCG_GC_w_RT_t <- subset(OCG_GC_w_RT, OCG_GC_w_RT$Subspecies=="T") #102

#subset to just wyomingensis
OCG_GC_w_RT_w <- subset(OCG_GC_w_RT, OCG_GC_w_RT$Subspecies=="W") #65

#subset to just vaseyana
OCG_GC_w_RT_v <- subset(OCG_GC_w_RT, OCG_GC_w_RT$Subspecies=="V") #50

#remove md from OCG_GC_w_RT to just have RT and Compounds 
rownames(OCG_GC_w_RT) <- OCG_GC_w_RT[,1]
OCG_GC_w_RT <- OCG_GC_w_RT[,-c(1:18)]

#re-organize GC_w_RT to be one column for RT and one for peak area correspinding to compound # ####
df_long <- OCG_GC_w_RT %>%
  pivot_longer(cols = starts_with("RT"), names_to = "Retention_Time_Column", values_to = "Retention_Time") %>%
  pivot_longer(cols = starts_with("C"), names_to = "Peak_Area_Column", values_to = "Peak_Area") %>%
  mutate(
    Retention_Time_Column = gsub("RT", "", Retention_Time_Column),  # Extract numeric part from RT column name
    Peak_Area_Column = gsub("C", "", Peak_Area_Column),  # Extract numeric part from C column name
    Compound = paste0("C", Peak_Area_Column)
  ) %>%
  filter(as.numeric(Retention_Time_Column) == as.numeric(Peak_Area_Column))  # Filter to ensure correspondence

# Remove unnecessary columns
df_long <- df_long[, c("Compound", "Retention_Time", "Peak_Area")]

#re-organize GC_w_RT 2012 to be one column for RT and one for peak area correspinding to compound # ####
#remove md from OCG_GC_w_RT to just have RT and Compounds 
OCG_GC_w_RT_2012 <- OCG_GC_w_RT_2012[,-c(1:18)]

df_long_2012 <- OCG_GC_w_RT_2012 %>%
  pivot_longer(cols = starts_with("RT"), names_to = "Retention_Time_Column", values_to = "Retention_Time") %>%
  pivot_longer(cols = starts_with("C"), names_to = "Peak_Area_Column", values_to = "Peak_Area") %>%
  mutate(
    Retention_Time_Column = gsub("RT", "", Retention_Time_Column),  # Extract numeric part from RT column name
    Peak_Area_Column = gsub("C", "", Peak_Area_Column),  # Extract numeric part from C column name
    Compound = paste0("C", Peak_Area_Column)
  ) %>%
  filter(as.numeric(Retention_Time_Column) == as.numeric(Peak_Area_Column))  # Filter to ensure correspondence

# Remove unnecessary columns
df_long_2012 <- df_long_2012[, c("Compound", "Retention_Time", "Peak_Area")]

#re-organize GC_w_RT 2021 to be one column for RT and one for peak area correspinding to compound # ####
#remove md from OCG_GC_w_RT to just have RT and Compounds 
OCG_GC_w_RT_2021 <- OCG_GC_w_RT_2021[,-c(1:16)]

df_long_2021 <- OCG_GC_w_RT_2021 %>%
  pivot_longer(cols = starts_with("RT"), names_to = "Retention_Time_Column", values_to = "Retention_Time") %>%
  pivot_longer(cols = starts_with("C"), names_to = "Peak_Area_Column", values_to = "Peak_Area") %>%
  mutate(
    Retention_Time_Column = gsub("RT", "", Retention_Time_Column),  # Extract numeric part from RT column name
    Peak_Area_Column = gsub("C", "", Peak_Area_Column),  # Extract numeric part from C column name
    Compound = paste0("C", Peak_Area_Column)
  ) %>%
  filter(as.numeric(Retention_Time_Column) == as.numeric(Peak_Area_Column))  # Filter to ensure correspondence

# Remove unnecessary columns
df_long_2021 <- df_long_2021[, c("Compound", "Retention_Time", "Peak_Area")]

#GC_w_RT for vaseyana #### 
OCG_GC_w_RT_v <- OCG_GC_w_RT_v[,-c(1:16)]

df_long_v <- OCG_GC_w_RT_v %>%
  pivot_longer(cols = starts_with("RT"), names_to = "Retention_Time_Column", values_to = "Retention_Time") %>%
  pivot_longer(cols = starts_with("C"), names_to = "Peak_Area_Column", values_to = "Peak_Area") %>%
  mutate(
    Retention_Time_Column = gsub("RT", "", Retention_Time_Column),  # Extract numeric part from RT column name
    Peak_Area_Column = gsub("C", "", Peak_Area_Column),  # Extract numeric part from C column name
    Compound = paste0("C", Peak_Area_Column)
  ) %>%
  filter(as.numeric(Retention_Time_Column) == as.numeric(Peak_Area_Column))  # Filter to ensure correspondence

# Remove unnecessary columns
df_long_v <- df_long_v[, c("Compound", "Retention_Time", "Peak_Area")]

#GC_w_RT for wyomingensis ####
OCG_GC_w_RT_w <- OCG_GC_w_RT_w[,-c(1:16)]

df_long_w <- OCG_GC_w_RT_w %>%
  pivot_longer(cols = starts_with("RT"), names_to = "Retention_Time_Column", values_to = "Retention_Time") %>%
  pivot_longer(cols = starts_with("C"), names_to = "Peak_Area_Column", values_to = "Peak_Area") %>%
  mutate(
    Retention_Time_Column = gsub("RT", "", Retention_Time_Column),  # Extract numeric part from RT column name
    Peak_Area_Column = gsub("C", "", Peak_Area_Column),  # Extract numeric part from C column name
    Compound = paste0("C", Peak_Area_Column)
  ) %>%
  filter(as.numeric(Retention_Time_Column) == as.numeric(Peak_Area_Column))  # Filter to ensure correspondence

# Remove unnecessary columns
df_long_w <- df_long_w[, c("Compound", "Retention_Time", "Peak_Area")]

#GC_w_RT for tridentata ####
OCG_GC_w_RT_t <- OCG_GC_w_RT_t[,-c(1:16)]

df_long_t <- OCG_GC_w_RT_t %>%
  pivot_longer(cols = starts_with("RT"), names_to = "Retention_Time_Column", values_to = "Retention_Time") %>%
  pivot_longer(cols = starts_with("C"), names_to = "Peak_Area_Column", values_to = "Peak_Area") %>%
  mutate(
    Retention_Time_Column = gsub("RT", "", Retention_Time_Column),  # Extract numeric part from RT column name
    Peak_Area_Column = gsub("C", "", Peak_Area_Column),  # Extract numeric part from C column name
    Compound = paste0("C", Peak_Area_Column)
  ) %>%
  filter(as.numeric(Retention_Time_Column) == as.numeric(Peak_Area_Column))  # Filter to ensure correspondence

# Remove unnecessary columns
df_long_t <- df_long_t[, c("Compound", "Retention_Time", "Peak_Area")]

#Plotting ####

#FULL GC
plot(df_long$Retention_Time, df_long$Peak_Area, type = "b", 
     xlab = "Retention Time (RT)", ylab = "Peak Area",
     main = "Peaks at Retention Times for all GC")

plot1 <- ggplot(df_long, aes(Retention_Time, Peak_Area, color = Compound))+
  geom_line(size = 1)+
  theme_clean()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(1,30, by = 2))+
  labs(x = "Retention time (min)", y = "Peak Area", title = "Full GC peak area at retention times") 

#2012 GC
plot(df_long_2012$Retention_Time, df_long_2012$Peak_Area, type = "b", 
     xlab = "Retention Time (RT)", ylab = "Peak Area",
     main = "Peaks at Retention Times for 2012 GC")

plot2 <- ggplot(df_long_2012, aes(Retention_Time, Peak_Area, color = Compound))+
  geom_line(size = 1)+
  theme_clean()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(1,30, by = 2))+
  labs(x = "Retention time (min)", y = "Peak Area", title = "2012 GC peak area at retention times")

#2021 GC
plot(df_long_2021$Retention_Time, df_long_2021$Peak_Area, type = "b", 
     xlab = "Retention Time (RT)", ylab = "Peak Area",
     main = "Peaks at Retention Times for 2021 GC")

plot3 <- ggplot(df_long_2021, aes(Retention_Time, Peak_Area, color = Compound))+
  geom_line(size = 1)+
  theme_clean()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(1,30, by = 2))+
  labs(x = "Retention time (min)", y = "Peak Area", title = "2021 GC peak area at retention times") 

grid.arrange(plot2, plot3, plot1, ncol = 1)

# WYOMINGENSIS
plot(df_long_w$Retention_Time, df_long_w$Peak_Area, type = "b", 
     xlab = "Retention Time (RT)", ylab = "Peak Area",
     main = "Peaks at Retention Times for tridentata GC")

plot4 <- ggplot(df_long_w, aes(Retention_Time, Peak_Area, color = Compound))+
  geom_line(size = 1)+
  theme_clean()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(1,30, by = 2))+
  labs(x = "Retention time (min)", y = "Peak Area", title = "GC data for wyomingensis")

# VASEYANA
plot(df_long_v$Retention_Time, df_long_v$Peak_Area, type = "b", 
     xlab = "Retention Time (RT)", ylab = "Peak Area",
     main = "Peaks at Retention Times for vaseyana GC")

plot5 <- ggplot(df_long_v, aes(Retention_Time, Peak_Area, color = Compound))+
  geom_line(size = 1)+
  theme_clean()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(1,30, by = 2))+
  labs(x = "Retention time (min)", y = "Peak Area", title = "GC data for vaseyana")

# TRIDENTATA
plot(df_long_t$Retention_Time, df_long_t$Peak_Area, type = "b", 
     xlab = "Retention Time (RT)", ylab = "Peak Area",
     main = "Peaks at Retention Times for tridentata GC")

plot6 <- ggplot(df_long_t, aes(Retention_Time, Peak_Area, color = Compound))+
  geom_line(size = 1)+
  theme_clean()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(1,30, by = 2))+
  labs(x = "Retention time (min)", y = "Peak Area", title = "GC data for tridentata")

grid.arrange(plot4, plot5, plot6, ncol = 1)

# COMPARING 2012 TO 2021 FOR INDIVIDUAL COMPOUNDS
plot(OCG_GC_w_RT_2012$RT003, OCG_GC_w_RT_2012$C003, type = "b", 
     xlab = "Retention Time (RT)", ylab = "Peak Area",
     main = "Peaks at Retention Times")

plot(OCG_GC_w_RT_2021$RT003, OCG_GC_w_RT_2021$C003, type = "b", 
     xlab = "Retention Time (RT)", ylab = "Peak Area",
     main = "Peaks at Retention Times")

#both of these RT are around 4.2 

plot(OCG_GC_w_RT_2012$RT002, OCG_GC_w_RT_2012$C002, type = "b", 
     xlab = "Retention Time (RT)", ylab = "Peak Area",
     main = "Peaks at Retention Times")

plot(OCG_GC_w_RT_2021$RT002, OCG_GC_w_RT_2021$C002, type = "b", 
     xlab = "Retention Time (RT)", ylab = "Peak Area",
     main = "Peaks at Retention Times")

#these ones are different

# LCMS RT plot
## 3UL LCMS RAW DATA READ IN
OCG_LCMS_3uL <- read.csv("data_csv/OCG_LCMS_3uL_cleaned_w_RT.csv", head=T, check.names = F,stringsAsFactors = T, row.names = 1) #120 of 929 variables

df_long_lcms <- OCG_LCMS_3uL %>%
  pivot_longer(cols = starts_with("RT"), names_to = "Retention_Time_Column", values_to = "Retention_Time") %>%
  pivot_longer(cols = starts_with("C"), names_to = "Peak_Area_Column", values_to = "Peak_Area") %>%
  mutate(
    Retention_Time_Column = gsub("RT", "", Retention_Time_Column),  # Extract numeric part from RT column name
    Peak_Area_Column = gsub("C", "", Peak_Area_Column),  # Extract numeric part from C column name
    Compound = paste0("C", Peak_Area_Column)
  ) %>%
  filter(as.numeric(Retention_Time_Column) == as.numeric(Peak_Area_Column))  # Filter to ensure correspondence

ggplot(df_long_lcms, aes(Retention_Time, Peak_Area, color = Compound))+
  geom_line()+
  theme_clean()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(1,30, by = 1))+
  labs(x = "Retention time (min)", y = "Peak Area", title = "LCMS")

# compounds that occur before 5 minutes (coumarins)
# compounds that occur after 5 minutes (flavonoids)

df_long_lcms <- df_long_lcms %>%
  mutate(compound_class = case_when(
    Retention_Time >= 0 & Retention_Time < 3 ~ "Coumarin",
    Retention_Time >= 3 & Retention_Time < 6 ~ "Glycosylated Coumarin",
    Retention_Time >= 6 & Retention_Time < 8 ~ "Flavonoid",
    Retention_Time >= 8 ~ "Glycosylated Flavonoid",
    TRUE ~ "Unknown"
  ))

df_long_lcms <- df_long_lcms %>%
  filter(compound_class != "Unknown")

df_long_lcms <- df_long_lcms[,-c(1:4)]

df_long_lcms$compound_class <- as.factor(df_long_lcms$compound_class)
str(df_long_lcms)

df_long_lcms <- df_long_lcms %>%
  distinct(Compound, .keep_all = TRUE)

write.csv(df_long_lcms, file = "data_csv/lcms_compound_class.csv",row.names = FALSE)



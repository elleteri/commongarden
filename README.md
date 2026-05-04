# commongarden
This repository contains the script for analyzing, and visualizing the LC-MS and GC phytochemistry data of sagebrush leaves from a commmon garden in Orchard, ID.

The goal of this study is to investigate plant chemistry in a common garden setting with co-occurring subspecies. We sampled sagebrush leaves from over 70 plants in both 2012 and 2021. For leaf chemical analyses we conducted Liquid Column Mass Spectrometry (LCMS) and Gas Chromatography (GC).

Script: 
1. orchard_phytochemistry_analysis_&_figures.R
   First an analysis of compound richness across subspecies:cytotype levels, ecoregion, and year while controlling for plant ID using a       generalized linear mixed model for both GC and LC-MS. Figures depicting compound richness across subspecies:cytotype groups and year for    both datasets are generated here.
   Then we prepared data for PCAs with scaling and threshold definition (see paper for details and tested betadispersion for GC and LC-MS.    Our PCA figures for GC and LC-MS data were then generated and PERMANOVAs to look at chemical variation across subspecies:cytotype and       ecoregion.
   We subsetted our data to subspecies:cytotypes for each chemistry dataset for partial mantel tests across sptial and climate data from       the origin populations the seeds were collected. We performed PCAs using selected climate data for the LC-MS T_2n subspecies:cytotype       group and plotted these figures using ordisurf.
   We used ANCOM-BC2 to test differential abundance of compounds for our T_2n group across ecoregions.
   We performed a Procrustes analysis to look for correlation between LC-MS and GC data.
   Lastly used shape files from https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states to created a map of our seed-origin populations in the Western US colored by Level III ecoregions. 
  
Folders: 
1. data_csv_format
   - contains all of the csv files needed for the above script
   a. us_eco_l3_state_boundaries
      - contains the necessary shape files to create the ecoregion map (Figure 1)

   



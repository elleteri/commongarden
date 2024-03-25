# commongarden
This repository contains scripts on cleaning, analyzing, and visualizing microbial and chemical 
ecology data. 

The goal of this study is to investigate the relationship between plant chemistry and 
leaf microbiome in a common garden setting with co-occurring subspecies. We sampled sagebrush leaves 
from over 70 plants in both 2012 and 2021. For leaf chemical analyses, we measured nitrogen and 
carbon stable isotopes, Liquid Column Mass Spectrometry (LCMS), as well as Gas Chromatography (GC). We did fungal
metabarcoding to examine microbial community composition across host subspecies, ploidy, and 
chemistry. Ongoing analyses will connect distinct chemical peaks with differences in microbial taxa. 

Script: 
1. OCG_data_cleaning_analyses_&_visualization.R
   - The first part of this script is cleaning a fungal amplicon sequencing table compiled from the plants sampled in the common garden. Analyses techniques include permanovas, pairwise adonis, shannon diversity, effective species number calculatins, and linear modeling. Visualization includes NMDS plots, and box plots. Pairwise adonis and PERMANOVAS were done to look at differences between different subspecies, ploidy number, and subspecies ploidy combined in their microbial community composition. To delve into fungal community composition seen on the plants, Metacoder and ANCOM are performed based on the ASV and taxonomy data.
   - After microbial analyses, this script contains chemistry data (LCMS and GC). LCMS (1uL & 3uL) and GC samples prepared from 2012 and 2021 are cleaned and subsetted to the existing metadata. PCA, NMDS. PcoA, Binary jaccard, permanovas, pairwise adonis are all performed to explore the chemical compostion between subspecies, ploidy, subspecies ploidy, and year.
   - This script then contains a procrustes analysis on the chemistry data and microbial data of the scripts above. A mantel test is run at the end. These analyses are performed to assess relationships hypothesized between the chemical and microbial data.
   - Stable isotope analyses between years and subspecies is done at the end of this script.
  
Folders: 
1. nmds
   - contains nmds rda files from the code to prevent having to rerun metaMDS.
2. jnmds
   - contains rda files from binary jaccards run from the code to prevent having to rerun. 
3. data_csv
   - contains all of the csv files needed for the above script
4. RAWDATAFILES
   - included until publication for cleaning and organizing raw data files. 


   



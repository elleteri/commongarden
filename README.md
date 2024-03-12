# commongarden
This repository contains scripts on cleaning, analyzing, and visualizing microbial and chemical 
ecology data. 

The goal of this study is to investigate the relationship between plant chemistry and 
leaf microbiome in a common garden setting with co-occurring subspecies. We sampled sagebrush leaves 
from over 70 plants in both 2012 and 2021. For leaf chemical analyses, we measured nitrogen and 
carbon stable isotopes, Liquid Column Mass Spectrometry (LCMS), as well as Gas Chromatography (GC). We did fungal
metabarcoding to examine microbial community composition across host subspecies, ploidy, and 
chemistry. Ongoing analyses will connect distinct chemical peaks with differences in microbial taxa. 

Scripts: 
1. Microbial_data_cleaning_analyses_&_visualization.R
   - This script is cleaning a fungal amplicon sequencing table compiled from the plants sampled in the common garden. Analyses techniques include permanovas, pairwise adonis, shannon diversity, effective species number calculatins, and linear modeling. Visualization includes NMDS plots, and box plots. Pairwise adonis and PERMANOVAS were done to look at differences between different subspecies, ploidy number, and subspecies ploidy combined in their microbial community composition.
2. Chem_data_cleaning_&_analysis.R
   - This script contains chemistry data (LCMS and GC). LCMS (1uL & 3uL) and GC samples prepared from 2012 and 2021 are cleaned and subsetted to th existing metadata. PCA, NMDS. PcoA, Binary jaccard, permanovas, pairwise adonis are all performed to explore the chemical compostion between subspecies, ploidy, subspecies ploidy, and year.
3. procrustes.R
   - This script contains a procrustes analysis on the chemistry data and microbial data of the scripts above. A mantel test is run at the end. 



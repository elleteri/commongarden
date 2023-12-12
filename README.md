# commongarden
This repository contains scripts on cleaning, analyzing, and visualizing microbial and chemical 
ecology data. 

The goal of this study is to investigate the relationship between plant chemistry and 
leaf microbiome in a common garden setting with co-occurring subspecies. We sampled sagebrush leaves 
from over 70 plants in both 2012 and 2021. For leaf chemical analyses, we measured nitrogen and 
carbon stable isotopes, as well as Gas Chromatography (GC). We did fungal and limited bacterial 
metabarcoding to examine microbial community composition across host subspecies, ploidy, and 
chemistry. Ongoing analyses will connect distinct chemical peaks with differences in microbial taxa. 

Scripts: 
1. Microbial_cleaning_OCG.Rmd
   - This script is cleaning a fungal amplicon sequencing table compiled from the plants sampled in the common garden.
2. Microbial analyses_&_visualization.Rmd
   - Using the cleaned and subsetted asv table and metadata, statistical analyses (Pairwise adonis and PERMANOVAS) were done to look at differences between different subspecies, ploidy number, and subspecies ploidy combined in their microbial community composition. 
3. Chemistry_2012_data_cleaning_OCG.Rmd
   - GC samples prepared from 2012 are cleaned and subsetted to just the plants that were part of the fungal amplicon sequencing collection.
4. Chemistry_2012_OCG_data_analyses_&_visualization.Rmd
   - Using the cleaned GC data, PCA plots, correlation tests and plots, PERMANOVAS, Pairwise adonis, and NMDS plots were done to look at the differences found between subspecies, ploidy number, and subspecies ploidy combined in their GC chemistry composition.
5. References

# O'Donnell et al. (2022) R Scripts

This repository contains R scripts used for figures and analysis in O'Donnell et al. (2022) Molecular and morphological analyses support recognition of Prostanthera volucris (Lamiaceae), a new species from the Central Tablelands of New South Wales 

MAP - This script produces a species distribution map with cartographic features, seen as Fig. 1 in O'Donnell et al. (2022)

DART - This script outlines the processing and analysis pipeline that was used on genomic SNP data delivered by DArTseq. This includes:
  - Data import and filtering
  - Identification and removal of potential clones
  - PCA ordination
  - Compute a distance matrix and export for Neighbour-Network analysis in SplitsTree
  - Calculate species and population level statistics
  - SNMF clustering analysis
  - SNMF plotting (bar plot and map with pie charts)
  - Create data subset for coalescent analysis
  - Export SNP matrix to SVDQuartets
  - Import and plot resultant phylogeny
  

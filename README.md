# Code and data for: 'Mapping shifts in spatial synchrony in grassland birds to inform conservation planning' (Allen & Lockwood, Conservation Biology 2020)
This code (i.e., the GLBsynchrony.Rmd file along with associated scripts "gridcollapse.func.R", "run.gamyet.func.R", and "mapsync.func.R") can be used to reproduce the analyses and figures from the associated manuscript (Allen & Lockwood 2020). It processes route-level North American Breeding Bird Survey count data summarized into 2x2 degree grid-cells to analyze spatial synchrony for 19 species of grassland birds. It relies heavily on data and code provided within the bbsBayes package (Edwards & Smith 2020; version 2.3.3.2020) and a Bayesian hierarchical model implemented within that package and described in Smith and Edwards (2020). 

IMPORTANT NOTE: this code uses North American Breeding Bird Survey data (Pardieck et al. 2020) that must be accessed through the BBSbayes package (Edwards & Smith 2020). The code was most recently tested on 2020-10-27 with R version 4.0.2 (2020-06-22), RStudio version 1.3.1056, and bbsBayes version 2.3.3.2020. Version information for other packages are provided in the comments of Section 1 and full session info is provided in Section 8 below.

ANOTHER IMPORTANT NOTE: Many of the computations performed by this code can take a long time to complete depending on your computer (e.g., up to 15 hours per species - time period combination). Therefore, we provide the option to skip Sections 1-3 below by pre-loading synchrony measurements calculated from each of 3000 data sets resampled from the BBS index posterior values (i.e., the values used in the analysis presented in Allen & Lockwood 2020). Alternatively, you can skip all four 'processing' sections (1-4) and proceed directly to the data visualization section (Sections 5 & 6). This is possible because we also pre-load the final summarized spatial synchrony data in Section 1 as an Rdata file (file name: sync.summaries.Rdata). 

## Minimum list of files required for all code to work properly
File names with relative file paths. Note this does not include files that are included as a means of skipping sections as described above.

"GLBsynchrony.Rmd" - this document, which acts as a guide for all analyses.

"scripts/gridcollapse.func.R" - function to collapse data into 2x2 grid cells.

"scripts/run.gamyet.func.R" - function to run JAGS models to generate BBS indices.

"scripts/mapsync.func.R" - a function to create many resampled data sets from BBS index posteriors and compute synchrony metrics for each.

"scripts/gamyet.bug" - the heavy-tailed GAMYE JAGS model script from bbsBayes package (Edwards & Smith 2020, Smith & Edwards 2020).

"data/jags/db.csv" - a list of 1x1 latitude/longitude grid cells with area information for JAGS model.

"data/shapefiles/BCR_Terrestrial_master_International_clip.shp" - a shapefile for mapping Bird Conservation Regions (also need associated .shx, .dbf, .prj files).

"data/shapefiles/province.shp" - a shapefile for mapping Canadian provinces (also need associated .shx, .dbf, .prj files).

## Table of contents for the main .Rmd script
Sections

1. Load libraries, functions, data

2. Setup and run JAGS model to generate BBS indices

3. Computing mean spatial synchrony by species and grid cell

4. Summarize spatial synchrony data for mapping and analysis

5. Create figures in the manuscript body

6. Create figures and recreate analysis from Supporting Information

7. References

8. Package versions and session info

## References

Allen MC, Lockwood JL. 2020. Mapping shifts in spatial synchrony in grassland birds to inform conservation planning. Conservation Biology. doi: 10.1111/cobi.13662.

Edwards BPM, Smith AC. 2020. bbsBayes: An R Package for Hierarchical
Bayesian Analysis of North American Breeding Bird Survey Data. bioRxiv
2020.05.27.118901. https://doi.org/10.1101/2020.05.27.118901

Pardieck KL, Ziolkowski Jr. DJ, Lutmerding M, Aponte VI, Hudson M-AR. 2020. North American Breeding Bird Survey Dataset 1966 - 2019:
U.S. Geological Survey data release, https://doi.org/10.5066/P9J6QUF6.
  
Smith AC, Edward BPM. 2020. Improved status and trend estimates from the North American Breeding Bird Survey using a Bayesian hierarchical generalized additive model. bioRxiv:2020.03.26.010215. 

The version used in the publication (v1.0.0) is archived on Zenodo: [![DOI](https://zenodo.org/badge/161780440.svg)](https://zenodo.org/badge/latestdoi/161780440) and also on OSF: https://osf.io/zrdcy/ . The OSF version has additional output files (csv files containing the full posterior samples).

# Trans-equatorial migration links oceanic frontal habitats across the Pacific Ocean: year-round movements and foraging activity of a small gadfly petrel 

Thomas A. Clay, Michael de L. Brooke

## Overview and data

The following repository contains codes to derive spatial locations based on light-based geolocation for Stejneger's petrels Pterodroma longirostris, in revision in Marine Biology. DOI: X

We based our process on the manual published by Lisovski et al. (2020) Light-Level Geolocator Analyses: A user’s guide. Journal of Animal Ecology 89:221-236. DOI: 10.1111/1365-2656.13036. See https://geolocationmanual.vogelwarte.ch/ for details. 

Raw light data and twilights, as well as processed geolocator locations, can be viewed and/or downloaded in the Movebank Data Repository: Clay TA, Brooke M de L. (2024). Data from: Trans-equatorial migration links oceanic frontal habitats across the Pacific Ocean: year-round movements and foraging activity of a small gadfly petrel. Movebank Data Repository. https://doi.org/10.5441/001/1.301

The timing of geolocator calibration periods are available in the "Data_inputs" folder. We suggest you paste the ".lux" files from the Movebank Data Repository into the "./Data_inputs/Raw Light" folder and twilight ".csv" files into the "./Data_outputs/Full twilights" folder.

## General statement (please read before using data)

The archived files (in the repostory listed above) contain data derived from fieldwork on Isla Alejandro Selkirk, Juan Fernández Islands, Chile. We request that you let Tommy Clay (tommy.clay@outlook.com) or Michael Brooke (mb10005@cam.ac.uk) know if you use them, for the following reasons:
- The data are complex and workers who do not know the study species are likely to benefit from advice on interpretation.
- People within the project or collaborators may be analyzing data from this project and it is desirable to prevent duplication of effort.

## Codes

Here are the main codes in R to reproduce the locations using the TwGeos and SGAT R packages in the aforementioned paper.

- **1_calibration_zeniths.r**: Code to crop light files to get calibration periods, identify twilights and obtain sun elevation angles (zeniths).
- **2_deriving_annotating_twilights_manual.r**: Code to manually derive twilights from the whole track and manually check and correct based on visualizing light curves.
- **3_deriving_locations_SGAT.r**: Code to load in twilights and zeniths, create land mask and movement prior and run SGAT iteratively for each individual. 

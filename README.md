# Species Distribution Models of the Marmosini species of Colombia

## Table of contents
* [Additional datasets] Additional datasets
* [How to run this workflow] How to run this workflow
* [Other software involved] Other software involved
* [Authors and participations] Authors and participations
* [Future work] Future work

This is the workflow for the paper entitled 'Distribution and conservation of the species of Marmosini (Didelphimorphia, Didelphidae) from Colombia'. It is an **R** workflow with some minor additions and modification from QGIS. This work uses MaxEnt models to delimit the distribution of 16 species of Marmosini (Didelphidae) from Colombia, which include the genera *Marmosa* and *Monodelphis*. The workflow part from data preparation and gathering and ends in a full reproducible scripts of our main results and supporting information. Please see below and read **METHODS** section of our manuscript (pending acceptance) for a detailed explanation of methods and requirements. 

A link to the preprint version of the manuscript can be found at: [](www.asda.com)

### **IMPORTANT:** THIS SCRIPT IS DESING TO FUNCTION WITH PACKAGE `ENMeval` < 2.0.0 

## Additional datasets

Following is the folders that the user should create in the base directory of the project, the content each of them should have and some comments of our election. Only format is required to be mantained but the data could vary according to user necessities. 

|Folder name | Content | Comments |
| --- | --- | ---- |
| LayersBank/ | Predictors rasters (.grd) layers. | This layes should be aling in the same origin and extent, and their resolution should be match at ~1 km² |
| MODIS/Terra/Vegetation_Index/ | MODIS Terra images of MOD13A3 1km Monthly | For this data set, we used the year 2000 but people replicating this work could choose any year |
|wmap/ | World map shapefile (.shp) from https://www.naturalearthdata.com/ | We used the 10 m level 0 (country level) data set |
| HFP/ | Human Index Data set in raster (.grd) format | In this case we used a Colombian HFP index for the year 2015 by [Correa Ayram et al. (2020)](https://www.sciencedirect.com/science/article/abs/pii/S1470160X20305677) |
| BiomeSHP/ | Ecoregion map in shapefile (.shp) format | We used the dataset by [Dinerstein et al. (2017)](https://academic.oup.com/bioscience/article/67/6/534/3102935) |
| WDPA/ | Word Database on Protected Areas in .shp format | We used the dataset from February 2021, available at https://www.protectedplanet.net/ |


## How to run this workflow
- You have two options, go one by one of the scripts or use the 'R/3_fit_models.R' archive to source automatically the previous scrpts. Either way, all other scripts should be run separately since they are not inter-dependent as the preivous ones. 

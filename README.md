# Species Distribution Models of the Marmosini species of Colombia

## Table of contents
* [Additional datasets](#additional-datasets)
* [How to run this workflow](#how-to-run-this-workflow)
* [Other software involved](#other-software-involved)
* [Authors and contributions](#authors-and-contributions)
* [Future work](#future-work)

---

This is the workflow for the paper entitled 'Distribution and conservation of the species of Marmosini (Didelphimorphia, Didelphidae) from Colombia'. It is an **R** workflow with some minor additions and modification from QGIS. This work uses MaxEnt models to delimit the distribution of 16 species of Marmosini (Didelphidae) from Colombia, which include the genera *Marmosa* and *Monodelphis*. The workflow part from data preparation and gathering and ends in a full reproducible scripts of our main results and supporting information. Please see below and read **METHODS** section of our manuscript (pending acceptance) for a detailed explanation of methods and requirements. 

A link to the preprint version of the manuscript can be found at: [here](https://doi.org/10.21203/rs.3.rs-557895/v1)

---

### **IMPORTANT:** THIS SCRIPT IS DESING TO FUNCTION WITH PACKAGE `ENMeval` < 2.0.0 

---

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

---

## How to run this workflow
You have two options, go one by one of the scripts or use the 'R/2_fit_models.R' archive to source automatically the previous scrpts. Either way, all other scripts should be run separately since they are not inter-dependent as the preivous ones. This script is design to function at the later stages with Colombia as the spatial target. However, this and other parameters can be tweak easely to replicate to other regions. 

---

## Other software involved
We used **QGIS 3.16** software to modified manually final models according to geographical barriers proposed (see [here](https://doi.org/10.21203/rs.3.rs-557895/v1) preprint text for details). Importantly, this script was run with **MaxEnt v. 3.4.0** through the `dismo` package of R. Similarly, we used **gdal 3.3.0** from R to make some calculations faster (compared to the package `raster` from R, but note that `terra` package has overcome this caveats too). 

---


## Authors and contributions
* This repository was desing, organized and commented by **Baltazar González** (baltazargch@gmail.com)
* Gabriel Martin and Federico Brook (co-authors of the manuscript associated with this workflow) participated in all the methodological and conceptual design.
* **NatureMap Argentina** provided the computer power necessary for this study

---


## Future work
This workflow will be apply with improvement to all other marsupial species of Colombia. All insights, suggestions and bugs or improvements are welcome at the first author email (baltazargch@gmail.com). 

---

### **PLEASE, FELL FREE TO COMMENT, SHARE AND IMPROVE THIS WORK**

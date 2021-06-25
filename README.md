![](images/readme_image.pdf)

# Species Distribution Models of the Marmosini species of Colombia

## Table of contents
* [Additional datasets](#additional-datasets)
* [How to run this workflow](#how-to-run-this-workflow)
* [Download the range data](#download-the-range-data)
* [Other software involved](#other-software-involved)
* [Authors and contributions](#authors-and-contributions)
* [Future work](#future-work)
* [How to cite](#how-to-cite)

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
| LayersBank/ | Predictors rasters (.grd) layers. | These layers should be aligned in the same origin and extent, and their resolution should match at ~1 km² |
| MODIS/Terra/Vegetation_Index/ | MODIS Terra images of MOD13A3 1km Monthly | For this data set we used the year 2000, but people replicating this work could choose any year |
|wmap/ | World map shapefile (.shp) from https://www.naturalearthdata.com/ | We used the 10 m level 0 (country level) data set |
| HFP/ | Human Index Data set in raster (.grd) format | In this case we used a Colombian HFP index for the year 2015 by [Correa Ayram et al. (2020)](https://www.sciencedirect.com/science/article/abs/pii/S1470160X20305677) |
| BiomeSHP/ | Ecoregion map in shapefile (.shp) format | We used the dataset by [Dinerstein et al. (2017)](https://academic.oup.com/bioscience/article/67/6/534/3102935) |
| WDPA/ | Word Database on Protected Areas in .shp format | We used the dataset from February 2021, available at https://www.protectedplanet.net/ |

---

## How to run this workflow
You have two options: go one by one with each of the scripts, or use the 'R/5_make_maps_for_revision.R' to source automatically the previous scripts. Either way, all other scripts should be ran separately since they are not inter-dependent as the previous ones. Specially the 'R/6_make_shapes_final_modif.R' where the object `goodmodels` should be chose manually by the author according to its expert criteria of maps and model performance. All model results will be placed in the folder *output* upong script excecution. Some figures are saved in the folder *results* upon runing the 'R/reults.R' script. These scripts are design to function at the later stages with Colombia as the spatial target. However, this and other parameters can be easily tweaked to replicate the analyses in other regions.

---

## Download the range data
The range maps that resulted from this work are available for download at: https://zenodo.org/record/4813016#.YMoPxXVKj0o

---

## Other software involved
We used **QGIS 3.16** software to manually modify the final models according to geographical barriers proposed (see here preprint text for details). Importantly, this script was ran with **MaxEnt v. 3.4.0** through the `dismo` package of R. Similarly, we used **gdal 3.3.0** from R to make some faster calculations (compared to the package `raster` from R, but note that `terra` package has overcome these caveats too).

---


## Authors and contributions
* This repository was designed, organized and commented by **Baltazar González** (baltazargch@gmail.com)
* **Gabriel Martin** and **Federico Brook** (co-authors of the manuscript associated with this workflow) participated in all the methodological and conceptual design
* NatureMap Argentina provided the computer power needed for this study

---


## Future work
This workflow will be applied with improvement to all other marsupial species of Colombia. All insights, suggestions and improvements (or bugs found) are welcomed at the first author's email (baltazargch@gmail.com).

---

## How to cite
Please, when using any part or the whole code give corresponding credits both for us and the packages used. Suggested cite to this recurse:

González, B., Brook, F., and Martin, G. (2021). Distribution and conservation of the species of Marmosini (Didelphimorphia, Didelphidae) from Colombia. PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-557895/v1].

**PLEASE, FELL FREE TO COMMENT, SHARE AND IMPROVE THIS WORK**

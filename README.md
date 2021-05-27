# Species Distribution Models of the Marmosin species of Colombia

This is the workflow for the paper entitled 'Distribution and conservation of the species of Marmosini (Didelphimorphia, Didelphidae) from Colombia'. It is an **R** workflow with some minor additions and modification from QGIS. This work uses MaxEnt models to delimit the distribution of 16 species of Marmosini (Didelphidae) from Colombia, which include the genera *Marmosa* and *Monodelphis*. The workflow part from data preparation and gathering and ends in a full reproducible scripts of our main results and supporting information. Please see below and read **METHODS** section of our manuscript (pending acceptance) for a detailed explanation of methods and requirements. 

A link to the preprint version of the manuscript can be found at: ![](www.asda.com)

## **IMPORTANT:** THIS SCRIPT IS DESING TO FUNCTION WITH PACKAGE <ENMeval> < 2.0.0 

## Additional datasets required to run this repository

- Layers: people who are interested in running this workflow should create a folder name *LayersBank* with the datasets declared in the manuscript. 
- MODIS data: additionally, to create the MSAVI layer, see the directory in script 'R/0_prepare_data.R' in the *MAKE MSAVI PREDICTOR LAYER* section.
- World Map: people interested in replicating this script should include the Natural Earth world map at 10 m level 0 in a folder name *wmap*. Map available for donwload at https://www.naturalearthdata.com/
- Humand Index and WDPA: similarly, those who will replicate this workflow should add two folders name *HFP* and *WDPA* with the Human Index and Protected Areas dataset declared in the manuscript named above (preprint). 
- Biomes: finally, a they should create a folder name *BiomeSHP* and include there the ecoregion dataset mentioned in the manuscript's methods.

## How to run this workflow
- You have two options, go one by one of the scripts or use the 'R/3_fit_models.R' archive to source automatically the previous scrpts. Either way, all other scripts should be run separately since they are not inter-dependent as the preivous ones. 

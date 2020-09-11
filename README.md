# Stier-Coral-Morphometrics-2020
Coral Morphometrics using Metashape

This project compares traditional field measurements of *Pocillopora* coral heads with 3D models constructed from digital images of the same coral using Agisoft Metashape software. We are interested in (1) creating a workflow utilizing best practices for photo acquistion and photogrammetry models of individual corals (2) compare information gained or lost from photogrammetry vs. in-field manual measurements and (3) the ecological relevance as a predictor for invertebrate communities within the coral. Data was collected in Moorea, French Polynesia on August 2019 for a project funded by the National Science Foundation (NSF).

This repo is maintained by Stier Lab Lead Technician Joseph Curtis (GitHub: [@jlscurtis512](https://github.com/jlscurtis512) and undergraduate student Journ Galvan (GitHub: [@journgalvan](https://github.com/journgalvan)) at the University of California, Santa Barbara in the Department of Ecology, Evolution, & Marine Biology.


# Morphometric Data
*/photogrammetry_data* Folder containing coral photogrammetry and field data
file name | description 
---|-----------
branchwidth_data.csv | Measurements of distances between coral branches including the average distance for each coral
photogrammetry_data_v3_2020_9_10.csv | Measurements estimated by the photogrammetry software
field_experiment_colony_measurements_moorea_summer2019.csv | Measurements taken in the field

# CAFI Data
*/cafi_data* Folder containing coral associated fish and invertebrate data
file name | description 
---|-----------
cafi_data_w_taxonomy_summer2019_2020_5_21.csv | contains most up to date and thorough taxonomy data
/old | subfolder with original invertebrate data sheets
/prelim_cafi_counts_moorea_summer2019.csv | Intitial counts of invertebrates 
/revised_cafi_data_moorea_summer2019_11_27.csv | Updated species level data on invertebrate communities

# Code
*/code* 
file name | description 
---|-----------
photogrammetry_CAFI_outline_JC_9_8_2020.Rmd | Synthesis and cleaning of data sets as well as code for figures and stats analysis comparing photogrammetry and traditional measurements. 
photogrammetry_CAFI_outline_invertonly_JC_9_8_2020.csv | Contains annotated suggestions and edits of the original code
/figures | Subfolder with code containing preliminary figures comparing software and manual measurements and how they relate to invertebrate communities



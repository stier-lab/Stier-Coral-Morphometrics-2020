# Stier-Coral-Morphometrics-2020
Coral Morphometrics using Metashape

This project compares traditional field measurements of *Pocillopora* coral heads with 3D models constructed from digital images of the same coral using Agisoft Metashape software. We are interested in (1) creating a workflow utilizing best practices for photo acquistion and photogrammetry models of individual corals (2) compare information gained or lost from photogrammetry vs. in-field manual measurements and (3) the ecological relevance as a predictor for invertebrate communities within the coral. Data was collected in Moorea, French Polynesia on August 2019 for a project funded by the National Science Foundation (NSF).

![](https://user-images.githubusercontent.com/47797235/113909750-e3eb5b00-978c-11eb-8c23-51aebe8c6d28.mp4)

*Note*: Project models and photos are stored in a GoogleDrive folder and are available upon request.

This repository is maintained by Ocean Recoveries Lab Manager Joseph Curtis (GitHub: [@jlscurtis512](https://github.com/jlscurtis512) and UCSB Alumni Journ Galvan (GitHub: [@journgalvan](https://github.com/journgalvan)) affiliated with the University of California, Santa Barbara in the Department of Ecology, Evolution, & Marine Biology.


# Morphometric Data
*/morphometric_data* Folder containing coral photogrammetry and field data
file name | description 
---|-----------
branchwidth_data.csv | Measurements of distances between coral branches including the average distance for each coral
photogrammetry_data_v3_2020_9_10.csv | Photogrammetric measurements from colonies imaged in August 2019
field_experiment_colony_measurements_moorea_summer2019.csv | Measurements taken in the field in August 2019
field_experiment_colony_measurements_moorea_december2019_v2.csv | Measurements taken in the field in December
maatea_experiment_photo_measurements_december_JC_2020_7_8.csv | Photogrammetric measurements from colonies imaged in December 2019

# CAFI Data
*/cafi_data* Folder containing coral associated fish and invertebrate data
file name | description 
---|-----------
cafi_data_w_taxonomy_summer2019_2020_5_21.csv | contains most up to date and thorough taxonomy data

# Code
*/code* 
file name | description 
---|-----------
photogrammetry_CAFI_figs_v3_2020_12_8.Rmd | Synthesis and cleaning of data sets as well as code for figures and stats analysis comparing photogrammetry and traditional measurements

photogrammetry_stats_aggregated_2021_4_7.R | Additional statistical analyses comparing photogrammetry and traditional measurements

/figures | Subfolder with figures generated by "photogrammetry_CAFI_figs_v3_2020_12_8.Rmd"

# Notes
*/notes* 
file name | description 
---|-----------
ORL_Photogrammetry_Protocol_Working.docx | Detailed manual explaining methodology and process for generating photogrammetric measurements from still images using Metashape

# Python Code for Batch Analysis in Metashape

Matthew Gottleib, a UCSB alumni, developed Python code for the Ocean Recoveries lab that can facilitate photogrammetry processing of multiple files by enabling Alignment and 3D model creation of batches of files without user input. Those files are stored and maintained on Matthew's GitHub: https://github.com/Mgla96/OceanRecoveryLabScripts

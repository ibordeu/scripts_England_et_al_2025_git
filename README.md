# Dysregulated stem cell co-option of the lung regeneration program drives tumor initiation. 
## England et al. Cell Stem Cell (2025)

[![DOI](https://zenodo.org/badge/914657790.svg)](https://doi.org/10.5281/zenodo.14673088)

This repository contains:

1. **data** : Folder with the clonal data and corresponding postprocessed data to reproduce the main and supplementary figures in the paper.

2. **plots** : Folder with the plots for the cumulative distribution of clone size, used in sim_two_pop_model scipt.

3. **data_analysis** : the scripts to reproduce all corresponding panels, as indicated. 

4. **sim_two_pop_model** : The script to run the simulations of the two-population model described in the Supplementary Note. This script produces plots comparing the empirical cumulative distribution of clone sizes with the one obtained numerically, as shown in the main and supplementary figures of the paper.

5. **utils** : Folder containing third-party plotting functions used in **data_analysis** to generate:\
 (i) density scatter plots [Dave B (2023). densityScatterChart (https://github.com/MATLAB-Graphics-and-App-Building/densityScatterChart/releases/tag/v1.2, GitHub. Retrieved 4 April 2023.]\
(ii) violin plots [bastibe (2023). Violinplot-Matlab https://github.com/bastibe/Violinplot-Matlab), GitHub. Retrieved 29 January 2023]

### To run the scripts:
-Download the scripts, together with **data** and **plots** folders, and save in a local directory.

-Open in **MATLAB** (tested in version 202a), change diretory to the location of the scripts and press **RUN**.

## Repository for Commonness as a reliable surrogacy strategy for the conservation planning of rare tree species in the subtropical Atlantic Forest
[ADD DOI LATER]

<!-- badges: start -->

[![License: GPL (\>=
2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)
[![Dependencies](https://img.shields.io/badge/dependencies-2/94-green?style=flat)](#)

<!-- badges: end -->

<br />

### Summary
This repository contains all the codes needed to replicate the analyses of the aforementioned paper. At the moment, data can only be made available upon reasonable request. Most modeling codes, which bear the prefix ```number_modeling_```roughly follow [Valavi et al. (2021)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1486), with modifications. It must be noted that some ```.rds``` filenames that may appear throughout the codes have **enfa** in their names. We did not proceed with enfa analyses but left that name there to avoid rewriting some parts of the codes. All codes work just fine and we commented pretty much everything, so replication should not be a problem (but if a issue is found, [please contact me](mailto:ggrittz@usp.br)).

<br />

### Workflow
Codes are structured as follows:  
1. **Preparing data**: scripts ```preparing_survey_data.R``` is used to produce site-by-species matrices for both common and rare species, as defined in the paper. ```preparing_env_data.R``` is used for geospatial analyses (crop and mask the study region etc.) and preparing all environmental variables through PCA analyses;  
2. **Loading data**: script ```loading_data.R``` is sourced to load everything that is needed for each modeling approach, that is, surveys, environmental variables, and a ```data.frame``` to save the results.
3. **Modeling common species**: scripts ending with one of the following suffixes were used to model common species: ```rf```, ```brt```, ```gam```, ```glm```, ```mars```, ```bart```, ```svm```, or ```xgboost```;
4. **Preparing matrices for Zonation and modeling performance** ```stack_and_performance.R``` is used to create stacked species distribution models (we did not use them but this part is available anyway), save each species distribution model as a matrix for later use in Zonation and evaluate each model performance through observed x predicted plots;
5. **Ensemble of Small Models**: all modeling related to rare species was performed within the script ```esm.R```;
6. **Zonation**: script ```zonation_analyses.R``` was used to prepare the data in the right format for the Zonation app and to analyze a bit the geospatial distribution of our data.

<br />

### How to cite

Please cite this GitHub as:

> **{{ ADD CITATION LATER }}**

<br />

### References
>Valavi, R., G. Guillera-Arroita, J. J. Lahoz-Monfort, and J. Elith. 2021. Predictive performance of presence-only species distribution models: a benchmark study with reproducible code. Ecological Monographs >00(00):e01486. 10.1002/ecm.1486

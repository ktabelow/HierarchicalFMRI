[![DOI](https://zenodo.org/badge/650516625.svg)](https://doi.org/10.5281/zenodo.16677811)

# HierarchicalFMRI

This repository contains the R scripts 

- realDataFigure2.R	
- simulationTable2.R

for the paper 

Neumann, A., Peitek, N., Brechmann, A., Tabelow, K., Dickhaus, T., 2021. Utilizing anatomical information for signal detection in functional magnetic resonance imaging. Preprint: Weierstraß-Institut für Angewandte Analysis und Stochastik. https://doi.org/10.20347/WIAS.PREPRINT.2806

## Required software

The scripts require an installation of R (https://cran.r-project.org)n and additional R packages that are available from CRAN:

- library(mutoss) - https://doi.org/10.32614/CRAN.package.mutoss
- library(ggplot2) - https://doi.org/10.32614/CRAN.package.ggplot2
- library(neuRosim) - https://doi.org/10.32614/CRAN.package.neuRosim
- library(fmri) - https://doi.org/10.32614/CRAN.package.fmri
- library(oro.nifti) - https://doi.org/10.32614/CRAN.package.oro.nifti

## Required data

The script realDataFigure2.R requires data that is available from 

Peitek, N., Brechmann, A., Tabelow, K., Dickhaus, T., 2024. Utilizing anatomical information for signal detection in functional magnetic resonance imaging - Data. https://doi.org/10.20347/wias.data.9

The script should automatically download the data from the resource. 

## Results

The script realDataFigure2.R creates the p-value in Figure 2 of the above mentioned paper. 
The script simulationTable2.R creates type-I and type-II error rates using different methods in three vectors: Aparc_rejected_over_MC (method from this paper), Cluster_rejected_over_MC (cluster-corrected inference), Mask_rejected_over_MC (voxel-wise inference).

## Tutorial

Both scripts, the simulation and real data analysis, can be adjusted to your own needs. Please refer to the [tutorial](tutorial.md) for more details.

## Authors

- Tabelow, Karsten, karsten.tabelow@wias-berlin.de, https://github.com/ktabelow 
- Brechmann, André 
- Dickhaus, Thorsten
- Peitek, Norman 


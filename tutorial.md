# Tutorial

## How to use the scripts on expertimental data

### Data preparation 

The script expects input data as csv file with the following columns (use the same names for the script to work). The data used in the paper is published at https://doi.org/10.20347/wias.data.9. If you want to use your own data, create the input data along the lines and use the script realDataFigure2.R to process it, see below.

#### Participants

Label for the subject. Usually, a study contains data from multiple subjects, which should be uniquely labelled. In our example data the labels are "subj1", "subj2" etc.

#### X, Y, Z

X, Y, and Z coordinate of the voxel under consideration in this line of the table.

#### Aparc

Label of spatial region from the corresponding atlas, an example is "ctx_lh_G_and_S_transv_frontopol". As the scripts where developed using Aparc labels, the column is named as such. However, any regional atlas can be used. (For compatibility reasons with the published data https://doi.org/10.20347/wias.data.9 we leave the column name as it is though.)

#### PC2_Maske_STAT

??? (This column is discarded by the script and only kept to be compatible with the published dataset at https://doi.org/10.20347/wias.data.9)

#### PC2_Maske_p

p-value ...

#### PC2_Compr_STAT

??? (This column is discarded by the script and only kept to be compatible with the published dataset at https://doi.org/10.20347/wias.data.9)

#### PC2_Compr_p

p-value ...

### Using the script realDataFigure2.R

To run the script on the data run this script with potential adjustments for your own data: Line 25 contains the location of the csv file described above and reads the data using R's read.csv() function. This can be an online ressource or some local file. If you generate your data already without column 6 and 8 (PC2_Maske_STAT, PC2_Compr_STAT) because they are not needed anyway, adjust the line 25 to
```
data <- read.csv("filelocation") 
```
so discard the indexing at the end.

The script will automatically detect the number of regional labels and the number of subjects. 

It will finally save a figure similar to the data part of Figure 4 of the paper https://doi.org/10.12688/f1000research.166549.1 to disk. Change the file name in line 162. (Figure 4 has been manually extended by images of the respective local regions)

A dump of output of the script is written to disk as file for expert examination, the file name is set in line 23.


## How to use simulation script

The simulation script first defines a function generate_fMRI_data() which generates data in lines 9-141. If you want to adjust this to other signal distributions and geometries, this function has to be rewritten. In its current form it takes the number of subjects of the simulation study and the signal SNR as arguments. The function creates simulated fmri data using the R package neuRosim and analyses the data using the R package fmri. As output it creates a matrix very similar to the experimental data described above. The difference is, that the columns 6, 7, 8, 9 contain the calculated voxelwise p-values for the first and second contrast (column 6 and 7) and the p-values for both contrast using cluster inference (column 8 and 9).

The main loop of the simulation is set up in lines 163-166, where you can choose the number of simulations, the number of spatial regions, the number of subjects, and the SNR for the data.

The output is given in the variables, from which the Type-I and Type-II errors such as in Table 2 of https://doi.org/10.12688/f1000research.166549.1 can be immediately calculated.
```
REGION_rejected_over_MC
CLUSTER_rejected_over_MC
VOXEL_rejected_over_MC
```
These vectors give, for all 8 considered regions in the simulation data, the number of simulation runs, where a regional hypotheses has been rejected. The labels REGION (method from this paper), CLUSTER (cluster-based detection), and VOXEL (voxelwise inference) in the variable name represent the utilized methods.

If you want to adjust the data generation, you can re-write the function generate_fMRI_data(). Minor changes of the dimension of the data (line 14) or the temporal design (lines 17-28) can be easily done adjusting the values. The original function just uses the R package neuRosim to create the fmri data, adjusting the geometry of the activation can be done by changing the parameters of the respective functions. The fmri data is then analysed by R package fmri.  


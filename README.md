# Transformed-Linear Prediction for Extremes

# Author Contributions Checklist Form

## Data

### Abstract
This manuscript analyzes daily levels of $\text{NO}_2$ measured in parts per billion (ppb) by the US Environmental Protection Agency (EPA).
The air pollution dataset contains geographical information for gauging stations in Washington DC.
Additionally, we analyze the daily UK precipitation of 30 stations from 1960 to 2024 across the Cumbria and Lancashire regions.

### Availability
To access the air pollution data, users can visit the EPA website and download the daily data at (https://www.epa.gov/outdoor-air-quality-data/download-daily-data). Similarly, the UK precipitation data are publicly available at Met Office (https://catalogue.ceda.ac.uk/uuid/8ddfd4dd5af443f9ad382cd77366d877).

### Description
The preprocessing of both datasets is described in Section 6 of the manuscript, and the preprocessed data can be found in the 'Data' folder. The EPA and the Met Office websites provide data dictionaries for both datasets.

## Code

### Abstract
The provided files offer guidance on how to apply transformed-linear prediction in both simulation studies and applications. By following the script files, users can generate optimized transformed-linear predictions, uncertainty quantification, and assess the coverage rate, etc.

### Description
The provided files contain comprehensive instructions and functions for applying transformed-linear prediction in simulation studies and applications, along with reproducing results in Sections 5 and 6 of the manuscript:

### 1. Script and its source files for simulation studies.

* A_Script_simulation.R: A script file that outlines the steps to reproduce results in the simulation study.
* functions.R: all source codes together.

* TransformedOperations.R: Functions that define transformed linear operations.
* PredictionError.R: A function that calculates the tail ratio of prediction error. 
* genDataParams.R: A function that simulates a regularly varying random vector X constructed from a matrix multiplication of A and returns a data matrix 'X_t' and relevant true quantities. 
* estimateParams.R: A function that estimates the TPDM and obtains the 'prediction' TPDM estimate.
* CPfactor.R: A function that implements completely positive decomposition for a 2x2 prediction (estimated) TPDM.
* JointRegion.R: A function that creates the 95% 'joint polar region' from normalized angular measure derived from the CP-factorization.
    * Note: To create the 95% joint polar region, Xhats in the test set must be calculated first. See line 65-66 in the script file.
* crossValidate.R: A function that finds the bandwidth for a kernel density estimation via k-fold cross-validation.
* KDE_w.R: A function that performs kernel density estimation for the angular density 'h' using angular components and masses from the CP-factorization with either the default bandwidth selection or the cross-validated bandwidth. 
    * Requires the "VGAM" and "ks" packages to perform kernel density estimation.
* condDensity.R: A function that provides an approximate conditional density given $\hat{X}$ is large using the kernel density estimate.
* coverageRate.R: A function that creates the 95% 'prediction interval' given large predicted values and assess the coverage rate.
* Mix_ECDF_GPD.R: A function that returns a mixture distribution of an empirical CDF below the threshold and a generalized Pareto distribution above the threshold.

### 2. Script and its source files for applications.
* A1_DataCleaningNO2.R: A script that creates a new data frame with daily $\text{NO}_2$ data from five stations, after detrending.
    * Requires the "readr" and "dplyr" packages to load and process multiple CSV files, and the "ismev" package to estimate parameters.
* A2_Mov_Avg.R: A function that calculates the moving average and moving standard deviation of a time series.
* A3_ScriptNO2Mixture.R: A script that reproduces the results of the air pollution application.
* A3_ScriptNO2MixtureMissing.R: A script that reproduces the results of the air pollution application when values at four stations are observed, but none are taken at Alexandria.
* A4_Script_UKprec.R: A script that reproduces the results of the UK precipitation application.

## Output files

### Rdata files for simulation study:

* CPout.RData: contains angular components and angular point masses obtained from the completely positive decomposition of the estimated prediction TPDM.

### Rdata files for applications:

####  Air pollution

* NewData.Rdata: contains preprocessed NO2 data for five stations.
* Train_NO2_Alx.Rdata: contains a training set with Pareto margins added.
* Test_NO2_Alx.Rdata: contains a test set with Pareto margins added.
* recentData.Rdata: contains a NO2 dataset for 2016-2020 when four stations are observed and no observation is taken at Alexandria.
* NewDat.Rdata: contains a dataset from 'NewData.Rdata' with Pareto margins added.

####  UK precipitation

* cumbriaLanc.RData: contains cleaned daily precipitation data across the Cumbria and Lancahsire regions.
* CP_prcp_final.Rdata: contains saved angles and angular masses obtained from the decomposition of the prediction TPDM estimate in a training set.

CPfactor.R implements completely positive decomposition to obtain multiple $2 \times q$ nonnegative matrices $C$. The computational time depends on the dimension of q* and the number of iterations in a for loop. By default, we choose q* to be 10 and set the number of iterations to 5000 until the algorithm converges. This function takes approximately 2-3 minutes to obtain about 500 angles and angular masses. If users want to skip this function, you can directly use the corresponding Rdata files in the **output** folder to reproduce the results.





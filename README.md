# NAM Model Validation for Tropical Cyclone Precipitation

North American Mesoscale (NAM) model validation based on forecasted precipitation for over 50 tropical storms and hurricanes since 2004. Precipitation forecasts compared with Stage IV (ST4) observational data product. Both the NAM forecast and the ST4 observational data are provided by NCEP/NOAA. We propose a Bayesian hierarchical framework to upgrade these forecasts from deterministic to probabilistic by leveraging information on the errors (difference between forecasted and observed precipitation) from previous TCs.

## Data structure 

As of August 2022, the raw forecast data and observational data products are available for download at https://www.ncei.noaa.gov/products/weather-climate-models/north-american-mesoscale and https://data.eol.ucar.edu/dataset/21.093, respectively. 

Upon downloading, tropical cyclone (TC) data will be organized within the directory, **NAMandST4**. Within this directory, there will be 47 subdirectories: one for each of the TCs from 2004 to 2017. Each of these storm directories contains two folders: one containing the NAM files for the storm, the other containing the ST4 files. All of these files are in gridded binary (GRIB) format, with most being in GRIB2 and some earlier storms (2006 or earlier) being GRIB (GRIB1). For example, Tropical Cyclone Arlene (2005) will have a directory `NAMandST4/2005arlene`, and contains subdirectories `NAMandST4/2005arlene/2005arlene` for the ST4 data and `NAMandST4/2005arlene/nam_218_2005061000arlene` for the NAM forecasts.

## Guidelines to reproduce results

### 1) SqrtPrecipVariograms.R (for data preprocessing)
This script reads in the downloaded data, based on the folder hierarchy described above. NAM and ST4 are on different resolutions (approximately 12km and 4km, respectively), so nearest neighbor interpolation is used upon the ST4 data to match the spatial resolution of the two datasets. This script also takes the square root transformation of the NAM and ST4, before subtracting the NAM from the ST4 to obtain an error field for each TC (which are written in .csv format). Plots are created for each of the 47 training storms (in .pdf format), and a simple variogram fit is also done. Note, this script should be run with `pred=F` as well as `pred=T` for the training and test storms, respectively.

This script can also calculate an average error for each grid point across the training storms (`makePWmean=T`, `subtractPWmean=T`), see Appendix of Walsh et al. (accepted for publication in Annals of Applied Statistics, October 2022).

### 2) MLE_storm_expntl_oneHess.R (calculate MLEs and Hessians)
Here, each error field is read in and the MLEs for the spatial parameters of the exponential covariance function are found, as well as the corresponding Hessian matrix (containing uncertainty information related to the MLEs). This script can be looped in a bash environment to run MLE calculations for each storm simultaneously. These are saved as .csvs. (For the nonspatial case, use `MLEs_nonsp_var.R`.)

### Optional) aggregate_covmtxs_expntl.R, combine_covmtxs.R and PWmean_post_expntl.R (Posterior calculations related to the pointwise mean)
We found assuming zero mean outperforms this estimation for this application. This is also computationally intensive, given there are over 26000 grid points when you aggregate all TC landfall locations.)

### 3) GibbsSamplerHurrSE_LM2.R (Gibbs sampler\*)
This takes MLEs and Hessian matrices as inputs and obtains posterior draws via MCMC for each theta_i, as well as mu_theta and Sigma_theta. Details on the modeling are provided in Walsh et al. Output saved as RData.

### 4) predictionNA_2018_SE_LM2.R (Prediction algorithm\*)
Using the Gibbs sampler output, uncertainty quantification for each test storm is computed by creating 1000 synthetic rainfall scenarios. This can also be looped in a bash environment to run simultaneously for the six test storms from 2018-2019 (therefore, this job should be run for `ste` values 1 through 6). 

### 5) laplace-metropolis.R (Initial model comparison)
To compare different versions of the modeling framework (a constant mean across all TCs, landfall-specific means for Atlantic, Florida and Gulf storms, etc), there are multiple versions of the codes marked with a \* above. These will have "", "LM2", "LM3" appended on the end to illustrate a different model was utilized in that script. LM2 performed the best, so we encourage that as a default in this application; see the manuscript for details.

### 6) SC_basins_mm_only.R (Obtain UQ on a basin-wide level)
For each test storm, find a distribution for the total precipitation over a 24 period. Different U.S. states are analyzed depending on the TC landfall location. Log scores are computed. This script should be ran for `ste` values 1 through 6, but storms 3 and 4 also have two different states for which it can be ran, therefore there is a total of 8 storm/state combinations to replicate the results in the manuscript.

### 7) scores_by_basin.R (Calculate log scores)
To obtain model comparisons, we look at the log scores for each proposed model, and look at the results overall, as well as broken up by region (ATL, FL, GULF) and by whether the basins are coastal or not.

### 8) NMV_plots.R (plots)
Reproduce the plots found in Walsh et al. (2022).

# NAM Model Validation for Hurricane Precipitation

North American Mesoscale (NAM) model validation based on precipitation observed over 40 tropical storms and hurricanes since 2005. Precipitation forecasts compared with Stage IV (ST4) observations. Both the NAM forecast and the ST4 observational data are provided by NCEP/NOAA.

Hurricane Data in .grib(2) format will be organized within the directory, **NAMandST4**. Within this directory, there will be 47 subdirectories: one for each of the tropical storms/hurricanes from 2004 to 2017. Each of these storm directories contains two folders: one containing the NAM files for the storm, the other containing the ST4 files. All of these files are in gridded binary (GRIB) format, with most being in GRIB2 and some earlier storms (2006 or earlier) being GRIB (GRIB1).

The NAM and ST4 are on different resolutions (approximately 12km and 4km, respectively), so nearest neighbor interpolation is used upon the ST4 data to match the spatial resolution of the two datasets.
 Goals include studying the uncertainty and errors of the NAM forecast by estimating spatial parameters of a Matern covariance function, utilizing a Gibbs sampler to generate posterior samples for each of these parameters, and use these samples for prediction purposes on storms making landfall in 2018 and 2019.


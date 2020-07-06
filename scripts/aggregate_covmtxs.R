# Work on aggregating precision data related to 
# the error fields

library(raster)

# 11 is default since fastest to calculate; loop through 1 to 47 in sbatch/linux
storm <- 11 

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

source("~/NAM-Model-Validation/scripts/mosaicList.R")

# load, modify and project land/sea mask 
mask <- raster("~/NAM-Model-Validation/lsmask.nc")
mask[mask==-1]  <- NA #changed from 0 to NA because mismatch rows due to off-coast pts
extent(mask)[1] <- extent(mask)[1]-360
extent(mask)[2] <- extent(mask)[2]-360
mask.regrid <- resample(mask, projectRaster(raster(
  "~/NAM-Model-Validation/nam_218_20050829_1200_f012.grib"),
  crs = "+proj=longlat +datum=WGS84"), method='ngb')  #/Volumes/LACIEHD/

# Read in error counts
count_files <- list.files("~/NAM-Model-Validation/error_rasters_counts/", pattern = ".grd", full.names = T)
error_counts <- mosaicList(count_files)
error_files <- list.files("~/NAM-Model-Validation/error_rasters/", pattern = ".grd", full.names = T)
error_sum <- mosaicList(error_files)

# only look at >= 1
geq1 <- error_counts
geq1[geq1 < 1] <- NA
# geq1[geq1 >= 1] <- 1
# plot(error_sum/geq1)

# Convert raster to data frame for ggplot
PWmean1 <- as(error_sum/geq1, "SpatialPixelsDataFrame")
PWM1_df <- as.data.frame(PWmean1)
colnames(PWM1_df) <- c("value", "x", "y")
# write.csv(PWM1_df, "~/NAM-Model-Validation/csv/error_df/PWmean_geq1.csv")
PWMEME <- read.csv("~/NAM-Model-Validation/csv/error_df/PWmean_geq1.csv", row.names = 1)

write.csv(data.frame(x = sprintf("%.20f", PWM1_df$x),
                     y = sprintf("%.20f", PWM1_df$y),
                     value = sprintf("%.20f", PWM1_df$value)),
                     "~/NAM-Model-Validation/csv/error_df/geq1_sprint.csv")

# only look at >= 20
geq20 <- error_counts
geq20[geq20 < 20] <- NA
# geq20[geq20 >= 20] <- 1
# plot(error_sum/geq20)

# Convert raster to data frame for ggplot
PWmean20 <- as(error_sum/geq20, "SpatialPixelsDataFrame")
PWM20_df <- as.data.frame(PWmean20)
colnames(PWM20_df) <- c("value", "x", "y")
# write.csv(PWM20_df, "~/NAM-Model-Validation/csv/error_df/PWmean_geq20.csv")


# Obtain parameter estimates for the portion of PWmean where counts >= 20
# ini.cov.likfit <- cbind(rep(c(0.5,1,1.5,2), each = 4), rep(c(0.05,0.5,1,2), 4)) #sigma2, phi
# fix.nug <- T; nug.val <- 0

# MLtemp <- list()
# likfit_time <- system.time({
#   MLtemp <- likfit(as.geodata(PWM20_df[,c(2,3,1)]), ini.cov.pars = ini.cov.likfit, lik.method = "ML",
#                    cov.model = "matern", fix.kappa = FALSE, fix.nugget = fix.nug, kappa=0.501,
#                    nugget = nug.val,  hessian = T)#, method = "L-BFGS-B", 
#   #lower = c(0,0,0), upper = c(10,10,10))
# })
# # save(MLtemp, file = "~/NAM-Model-Validation/RData/PWmean_geq20_paramsMLE.RData")

# loads MLtemp from above; takes 5-8 hours to run above code
load("~/NAM-Model-Validation/RData/PWmean_geq20_paramsMLE.RData") 

# load MLE estimates, etc from all_storm_res
all_storm_res <- read.csv("~/NAM-Model-Validation/csv/all_storm_res.csv", row.names = 1)

# Read in all error fields (in csv format)
dfs <- list.files("~/NAM-Model-Validation/csv/error_df/subtractPWmeanF", full.names = T)

# Count the number of pixels (rows) in each csv file
n_pixels <- c()
for (i in 1:length(dfs)) { n_pixels[i] <- dim(read.csv(dfs[i]))[1] }
# hist(n_pixels)
dfs[n_pixels < 2000] #11, 17, 18 are smallest areas

# storms 11 and 18 are similar locations and small pixel amts
kat1 <- read.csv(dfs[11], row.names = 1)
alb1 <- read.csv(dfs[16], row.names = 1)
ber1 <- read.csv(dfs[17], row.names = 1)
ern1 <- read.csv(dfs[18], row.names = 1)

# Returns intersection of locns
# plot(rasterFromXYZ(alb1[,c(2,3,1)]) - rasterFromXYZ(kat1[,c(2,3,1)]))

# Look at values/dimensions of raster
# dim(rasterFromXYZ(kat1[c(2,3,1)])) # of unique x, # of unique y, 1
# length(values((rasterFromXYZ(kat1[c(2,3,1)]))))


## Making my own dist to see if it is the same as dist; it is
# deet <- matrix(NA, nrow(s), nrow(s))
# for (i in 1:nrow(s)) {
#   for(j in 1:nrow(s)) {
#     deet[i,j] <- sqrt((s[i,1]-s[j,1])^2 + (s[i,2] - s[j,2])^2)
#   }
# }
# all.equal(deet, D)

# lats <- lons <- c() 
# for (i in 1:length(dfs)) {
#   this_lon <- read.csv(dfs[i], row.names = 1)[,2]
#   this_lat <- read.csv(dfs[i], row.names = 1)[,3]
#   lats <- c(lats, this_lat)
#   lons <- c(lons, this_lon)
# }




##############################
# Matern covariance function #
##############################

library(geoR)

# Without nugget effect

a <- nrow(PWM1_df)
# a <- round(a/2)
all_prc_mtx <- matrix(0, a, a)

# For the prior covariance structure using geq20 params
s <- PWM1_df[,c(2,3)]
D <- as.matrix(dist(s))
sigma2 <- MLtemp$sigmasq
phi <- MLtemp$phi
nu <- MLtemp$kappa
# plot(seq(0,10,0.01),sigma2*matern(seq(0,10,0.01),phi=phi,kappa=nu),type="l") #(1/phi)
priorcovmatrix <- sigma2*matern(D,phi=phi,kappa=nu) # + diag(tau2,nrow(t))
prior.prec.mtx <- solve(priorcovmatrix)



# For i in 1:47 storms; takes about a day when running via sbatch
for (st in storm) { # which(n_pixels < 3000)
  print(st)
  df <- read.csv(dfs[st], row.names = 1) #11 is smallest
  s <- df[,c(2,3)]
  D <- as.matrix(dist(s))
  sigma2 <- all_storm_res[st,"MLEsigma2"]
  phi <- all_storm_res[st,"MLEphi"]
  nu <- all_storm_res[st,"MLEkappa"]
  # plot(seq(0,10,0.01),sigma2*matern(seq(0,10,0.01),phi=phi,kappa=nu),type="l") #(1/phi)
  covmatrix <- sigma2*matern(D,phi=phi,kappa=nu) # + diag(tau2,nrow(t))
  prec.mtx <- solve(covmatrix)
  
  for (i in 1:nrow(prec.mtx)) {
    for (j in 1:i) {

      # which(arl1$x == kat1$x[10] & arl1$y == kat1$y[10])
      k <- which(abs(PWM1_df$x - df$x[i]) < 1e-3 & abs(PWM1_df$y - df$y[i]) < 1e-3)
      l <- which(abs(PWM1_df$x - df$x[j]) < 1e-3 & abs(PWM1_df$y - df$y[j]) < 1e-3)
      
      all_prc_mtx[k,l] <- all_prc_mtx[l,k] <- all_prc_mtx[k,l] + prec.mtx[i,j]
    }
  }
}

save(all_prc_mtx, file=paste0("~/NAM-Model-Validation/RData/all_prc_mtx/",ifelse(storm<10, paste0("0",storm),storm),".RData"))

# a <- dim(all_prc_mtx)[1]
sum_prec_mtx <- prior.prec.mtx
precs <- list.files("~/NAM-Model-Validation/RData/all_prc_mtx/", full.names = T)
for(i in 1:length(precs)){
  print(i)
  load(precs[i])
  sum_prec_mtx <- sum_prec_mtx + all_prc_mtx
}

system.time(sum_cov_mtx <- solve(sum_prec_mtx))
post_cov_mtx <- sum_cov_mtx
save(post_cov_mtx, file="~/NAM-Model-Validation/RData/post_cov_mtx.RData")

times <- c(1140, 1155, 1156, 1208, 1209, 1210, 1224, 1254, 1301, 1313, 1314, 1410, 1420, 1422, 1443, 
           1450, 1502, 1512, 1518, 1522, 1528, 1538, 1541, 1548, 1551, 1601, 1626, 1634, 1645, 1659, 
           1714, 1721, 1735, 1736, 1807, 1808, 1853, 1901, 1915, 1915, 2042, 2116, 2147, 2242, 2502,
           2508, 2724) - 1130
plot(sort(n_pixels), times, main="Size vs Compute Time")
plot(sort(n_pixels), sqrt(times), main="Size vs Sqrt Time")
abline(lm(sqrt(times)~sort(n_pixels)))

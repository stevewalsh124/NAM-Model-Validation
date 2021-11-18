library(knitr)
library(dplyr)
library(kableExtra)

# Model 1
avg_cov95_LM1 <- avg_cov99_LM1 <- npix_LM1 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"subtractpw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  avg_cov95_LM1[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM1[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM1[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM1*npix_LM1)/sum(npix_LM1)
sum(avg_cov99_LM1*npix_LM1)/sum(npix_LM1)

# Model 2
avg_cov95_LM2 <- avg_cov99_LM2 <- npix_LM2 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"_LM2_subtractpw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  avg_cov95_LM2[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM2[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM2[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM2*npix_LM2)/sum(npix_LM2)
sum(avg_cov99_LM2*npix_LM2)/sum(npix_LM2)

# Model 3
avg_cov95_LM3 <- avg_cov99_LM3 <- npix_LM3 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"_LM3_subtractpw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  avg_cov95_LM3[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM3[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM3[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM3*npix_LM3)/sum(npix_LM3)
sum(avg_cov99_LM3*npix_LM3)/sum(npix_LM3)

# Model 4
avg_cov95_LM4 <- avg_cov99_LM4 <- npix_LM4 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"nopw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  ####
  if(nrow(NAM_pred) > nrow(ST4_pred)){
    ST40 <- ST4_r
    values(ST40)[!is.na(values(ST40))] <- 0
    UB_r <- rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, ub_rain)) - ST40
    ub_rain <- as.data.frame(rasterToPoints(UB_r))$layer
    UB_r99 <- rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, ub_rain99)) - ST40
    ub_rain99 <- as.data.frame(rasterToPoints(UB_r99))$layer
    NAM_r <- NAM_r - ST40
    NAM_pred <- as.data.frame(rasterToPoints(NAM_r))
    colnames(NAM_pred) <- c("x","y","value")
  }
  
  ## if ST4 has more rows
  if(nrow(NAM_pred) < nrow(ST4_pred)){
    NAM0 <- NAM_r
    values(NAM0)[!is.na(values(NAM0))] <- 0
    UB_r <- rasterFromXYZ(cbind(ST4_pred$x, ST4_pred$y, ub_rain)) - NAM0
    ub_rain <- as.data.frame(rasterToPoints(UB_r))$layer
    UB_r99 <- rasterFromXYZ(cbind(ST4_pred$x, ST4_pred$y, ub_rain99)) - NAM0
    ub_rain99 <- as.data.frame(rasterToPoints(UB_r99))$layer 
    ST4_r <- ST4_r - NAM0
    ST4_pred <- as.data.frame(rasterToPoints(ST4_r))
    colnames(ST4_pred) <- c("x","y","value")
  } 
  ####
  avg_cov95_LM4[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM4[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM4[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM4*npix_LM4)/sum(npix_LM4)
sum(avg_cov99_LM4*npix_LM4)/sum(npix_LM4)

# Model 5
avg_cov95_LM5 <- avg_cov99_LM5 <- npix_LM5 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"_LM2_nopw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  ####
  if(nrow(NAM_pred) > nrow(ST4_pred)){
    ST40 <- ST4_r
    values(ST40)[!is.na(values(ST40))] <- 0
    UB_r <- rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, ub_rain)) - ST40
    ub_rain <- as.data.frame(rasterToPoints(UB_r))$layer
    UB_r99 <- rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, ub_rain99)) - ST40
    ub_rain99 <- as.data.frame(rasterToPoints(UB_r99))$layer
    NAM_r <- NAM_r - ST40
    NAM_pred <- as.data.frame(rasterToPoints(NAM_r))
    colnames(NAM_pred) <- c("x","y","value")
  }
  
  ## if ST4 has more rows
  if(nrow(NAM_pred) < nrow(ST4_pred)){
    NAM0 <- NAM_r
    values(NAM0)[!is.na(values(NAM0))] <- 0
    UB_r <- rasterFromXYZ(cbind(ST4_pred$x, ST4_pred$y, ub_rain)) - NAM0
    ub_rain <- as.data.frame(rasterToPoints(UB_r))$layer
    UB_r99 <- rasterFromXYZ(cbind(ST4_pred$x, ST4_pred$y, ub_rain99)) - NAM0
    ub_rain99 <- as.data.frame(rasterToPoints(UB_r99))$layer
    ST4_r <- ST4_r - NAM0
    ST4_pred <- as.data.frame(rasterToPoints(ST4_r))
    colnames(ST4_pred) <- c("x","y","value")
    
  } 
  ####
  avg_cov95_LM5[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM5[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM5[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM5*npix_LM5)/sum(npix_LM5)
sum(avg_cov99_LM5*npix_LM5)/sum(npix_LM5)

# Model 6
avg_cov95_LM6 <- avg_cov99_LM6 <- npix_LM6 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"_LM3_nopw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  ####
  if(nrow(NAM_pred) > nrow(ST4_pred)){
    ST40 <- ST4_r
    values(ST40)[!is.na(values(ST40))] <- 0
    UB_r <- rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, ub_rain)) - ST40
    ub_rain <- as.data.frame(rasterToPoints(UB_r))$layer
    UB_r99 <- rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, ub_rain99)) - ST40
    ub_rain99 <- as.data.frame(rasterToPoints(UB_r99))$layer
    NAM_r <- NAM_r - ST40
    NAM_pred <- as.data.frame(rasterToPoints(NAM_r))
    colnames(NAM_pred) <- c("x","y","value")
    
  }
  
  ## if ST4 has more rows
  if(nrow(NAM_pred) < nrow(ST4_pred)){
    NAM0 <- NAM_r
    values(NAM0)[!is.na(values(NAM0))] <- 0
    UB_r <- rasterFromXYZ(cbind(ST4_pred$x, ST4_pred$y, ub_rain)) - NAM0
    ub_rain <- as.data.frame(rasterToPoints(UB_r))$layer
    UB_r99 <- rasterFromXYZ(cbind(ST4_pred$x, ST4_pred$y, ub_rain99)) - NAM0
    ub_rain99 <- as.data.frame(rasterToPoints(UB_r99))$layer
    ST4_r <- ST4_r - NAM0
    ST4_pred <- as.data.frame(rasterToPoints(ST4_r))
    colnames(ST4_pred) <- c("x","y","value")
    
  } 
  ####
  avg_cov95_LM6[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM6[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM6[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM6*npix_LM6)/sum(npix_LM6)
sum(avg_cov99_LM6*npix_LM6)/sum(npix_LM6)

# Model 7
avg_cov95_LM7 <- avg_cov99_LM7 <- npix_LM7 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"nosp_nopw"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  avg_cov95_LM7[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM7[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM7[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM7*npix_LM7)/sum(npix_LM7)
sum(avg_cov99_LM7*npix_LM7)/sum(npix_LM7)

all_avg_covs <- cbind(avg_cov95_LM1, avg_cov99_LM1, 
                      avg_cov95_LM2, avg_cov99_LM2,
                      avg_cov95_LM3, avg_cov99_LM3,
                      avg_cov95_LM4, avg_cov99_LM4,
                      avg_cov95_LM5, avg_cov99_LM5,
                      avg_cov95_LM6, avg_cov99_LM6,
                      avg_cov95_LM7, avg_cov99_LM7)
all_avg_covs

avg_by_type <- c(sum(avg_cov95_LM1*npix_LM1)/sum(npix_LM1),
                 sum(avg_cov99_LM1*npix_LM1)/sum(npix_LM1),
                 sum(avg_cov95_LM2*npix_LM2)/sum(npix_LM2),
                 sum(avg_cov99_LM2*npix_LM2)/sum(npix_LM2),
                 sum(avg_cov95_LM3*npix_LM3)/sum(npix_LM3),
                 sum(avg_cov99_LM3*npix_LM3)/sum(npix_LM3),
                 sum(avg_cov95_LM4*npix_LM4)/sum(npix_LM4),
                 sum(avg_cov99_LM4*npix_LM4)/sum(npix_LM4),
                 sum(avg_cov95_LM5*npix_LM5)/sum(npix_LM5),
                 sum(avg_cov99_LM5*npix_LM5)/sum(npix_LM5),
                 sum(avg_cov95_LM6*npix_LM6)/sum(npix_LM6),
                 sum(avg_cov99_LM6*npix_LM6)/sum(npix_LM6),
                 sum(avg_cov95_LM7*npix_LM7)/sum(npix_LM7),
                 sum(avg_cov99_LM7*npix_LM7)/sum(npix_LM7))

my_table <- rbind(all_avg_covs, avg_by_type)
save(my_table, file="~/NAM-Model-Validation/RData/table/cov_storm_and_type.RData")

my_table %>%
  kbl(caption="Summary Statistics of Financial Well-Being Score by Gender and Education",
      format="latex",
      col.names = c("Gender","Education","Count","Mean","Median","SD"),
      align="r") %>%
  kable_minimal(full_width = F,  html_font = "Source Sans Pro")

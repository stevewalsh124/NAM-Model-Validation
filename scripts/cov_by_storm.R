# Find the weighted averages of coverages for each of the 9 models
# This saves my_table, which can be called into cov_tables.R to make LaTeX table
# and also in prediction_eval to make a supplemented table.

library(knitr)
library(dplyr)
library(kableExtra)

# Model 1
avg_cov95_LM1 <- avg_cov99_LM1 <- npix_LM1 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"nopw.RData"))
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
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"_LM2_nopw.RData"))
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
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"_LM3_nopw.RData"))
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
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"_LM4_nopw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  avg_cov95_LM4[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM4[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM4[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM4*npix_LM4)/sum(npix_LM4)
sum(avg_cov99_LM4*npix_LM4)/sum(npix_LM4)

# Model 5
avg_cov95_LM5 <- avg_cov99_LM5 <- npix_LM5 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"nosp_nopw"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  avg_cov95_LM5[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM5[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM5[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM5*npix_LM5)/sum(npix_LM5)
sum(avg_cov99_LM5*npix_LM5)/sum(npix_LM5)

# Model 6
avg_cov95_LM6 <- avg_cov99_LM6 <- npix_LM6 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"subtractpw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  avg_cov95_LM6[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM6[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM6[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM6*npix_LM6)/sum(npix_LM6)
sum(avg_cov99_LM6*npix_LM6)/sum(npix_LM6)

# Model 7
avg_cov95_LM7 <- avg_cov99_LM7 <- npix_LM7 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"_LM2_subtractpw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  avg_cov95_LM7[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM7[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM7[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM7*npix_LM7)/sum(npix_LM7)
sum(avg_cov99_LM7*npix_LM7)/sum(npix_LM7)

# Model 8
avg_cov95_LM8 <- avg_cov99_LM8 <- npix_LM8 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"_LM3_subtractpw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  avg_cov95_LM8[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM8[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM8[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM8*npix_LM8)/sum(npix_LM8)
sum(avg_cov99_LM8*npix_LM8)/sum(npix_LM8)

# Model 9
avg_cov95_LM9 <- avg_cov99_LM9 <- npix_LM9 <- c()
for (gr in 1:6) {
  load(paste0("~/NAM-Model-Validation/RData/prediction/prediction",gr,"_LM4_subtractpw.RData"))
  print(all.equal(length(NAM_pred$value) , length(ST4_pred$value)))
  avg_cov95_LM9[gr] <- mean(ub_rain + NAM_pred$value > ST4_pred$value)
  avg_cov99_LM9[gr] <- mean(ub_rain99 + NAM_pred$value > ST4_pred$value)
  npix_LM9[gr] <- nrow(NAM_pred)
}
sum(avg_cov95_LM9*npix_LM9)/sum(npix_LM9)
sum(avg_cov99_LM9*npix_LM9)/sum(npix_LM9)

# non-weighted average
all_avg_covs <- cbind(avg_cov95_LM1, avg_cov99_LM1, 
                      avg_cov95_LM2, avg_cov99_LM2,
                      avg_cov95_LM3, avg_cov99_LM3,
                      avg_cov95_LM4, avg_cov99_LM4,
                      avg_cov95_LM5, avg_cov99_LM5,
                      avg_cov95_LM6, avg_cov99_LM6,
                      avg_cov95_LM7, avg_cov99_LM7,
                      avg_cov95_LM8, avg_cov99_LM8,
                      avg_cov95_LM9, avg_cov99_LM9)
all_avg_covs

# weighted average
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
                 sum(avg_cov99_LM7*npix_LM7)/sum(npix_LM7),
                 sum(avg_cov95_LM8*npix_LM8)/sum(npix_LM8),
                 sum(avg_cov99_LM8*npix_LM8)/sum(npix_LM8),
                 sum(avg_cov95_LM9*npix_LM9)/sum(npix_LM9),
                 sum(avg_cov99_LM9*npix_LM9)/sum(npix_LM9))

my_table <- rbind(all_avg_covs, avg_by_type)
save(my_table, file="~/NAM-Model-Validation/RData/table/cov_storm_and_type.RData")

# Use my_table locally (off ARC) to create LaTeX tables (see cov_tables.R)
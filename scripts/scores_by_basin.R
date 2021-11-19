# evaluate basin predictive densities with log score

logS_files <- list.files("~/NAM-Model-Validation/csv/scores/basinwide", full.names = T, pattern = "storm[1-6].csv")
CRPS_files <- list.files("~/NAM-Model-Validation/csv/scores/basinwide", full.names = T, pattern = "CRPS")
logS_marf <- list.files("~/NAM-Model-Validation/csv/scores/basinwide/alt", full.names = T, pattern = ".csv")
logS_sum_marf <- list.files("~/NAM-Model-Validation/csv/scores/basinwide/altsum", full.names = T, pattern = ".csv")

logScores <- do.call(rbind, lapply(logS_files, read.csv, row.names = 1))
logScores <- logScores[complete.cases(logScores),]
CRPScores <- do.call(rbind, lapply(CRPS_files, read.csv, row.names = 1))
CRPScores <- CRPScores[complete.cases(CRPScores),]
logScores_M <- do.call(rbind, lapply(logS_marf, read.csv, row.names = 1))
logScores_M <- logScores_M[complete.cases(logScores_M),]
logScores_SM <- do.call(rbind, lapply(logS_sum_marf, read.csv, row.names = 1))
logScores_SM <- logScores_SM[complete.cases(logScores_SM),]

# find the mean wrt the original (not log) scale
f <- function(X){log(mean(exp(X - max(X)))) + max(X)}
f2 <- function(X){log(median(exp(X - max(X)))) + max(X)}

# mean of logScores
apply(logScores, 2, f) #V5 is best
apply(logScores_M, 2, f) #V4 is best, V5 is 2nd best
apply(CRPScores, 2, f) #V4 is best, V5 is 2nd best
# apply(logScores_SM, 2, f) #V4 is best, V5 is 2nd best

# best logScore overall
colSums(logScores) #V5 is best again
colSums(logScores_M) #V5 is best
colSums(CRPScores)
# colSums(logScores_SM) #V5 is best

# # median of logScores # not as useful 
# apply(logScores, 2, median) #V6 is best; doesn't penalize its shorter tails as much
# apply(logScores_M, 2, median) #V5 is best, V4 is 2nd best
# apply(logScores, 2, median) #V5 is best, V4 is 2nd best

# THESE ARE NOT IDEAL SINCE MEAN, ETC SHOULD BE CALCULATED ON THE ORIGINAL,
# NOT LOG, SCALE. SO, MAKE THIS ADJUSTMENT FIRST WITH F ABOVE
# colMeans(logScores, na.rm = T)
# apply(logScores, 2, median, na.rm=T)
# 
# the_best <- c()
# for (i in 1:length(logS_files)) {
#   logScore <- read.csv(logS_files[i], row.names = 1)
#   the_best[i] <- which.min(colMeans(logScore, na.rm = T))
# }
# the_best
# table(the_best)
# 
# the_best2 <- c()
# for (i in 1:nrow(logScores)) {
#   if(is.na(logScores[i,])) next
#   the_best2[i] <- which.min(logScores[i,])
# }
# the_best2
# table(the_best2)

load("~/NAM-Model-Validation/RData/prediction/prediction1subtractpw.RData")
dim(simvals)
sum(simvals <= 0)
plot(rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, NAM_pred$value)))


load("~/NAM-Model-Validation/RData/prediction/prediction1nopw")
dim(simvals)
sum(simvals <= 0)
plot(rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, NAM_pred$value)))

meanz <- rowMeans(simvals)
plot(rasterFromXYZ(cbind(NAM_pred$x, NAM_pred$y, meanz+NAM_pred$value)))

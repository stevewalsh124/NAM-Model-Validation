# evaluate basin predictive densities with log score (overall)

logS_files <- list.files("~/NAM-Model-Validation/csv/scores/basinwide/", full.names = T, pattern = "storm[1-6].csv")
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

# find the mean wrt the original (NOT log) scale
f <- function(X){log(mean(exp(X - max(X)))) + max(X)}
f2 <- function(X){log(median(exp(X - max(X)))) + max(X)}

# # mean of logScores
# sort(-apply(logScores, 2, f)) #V5 is best
# sort(apply(logScores_M, 2, f)) #V4 is best, V5 is 2nd best
# sort(-apply(CRPScores, 2, f)) #V4 is best, V5 is 2nd best
# # apply(logScores_SM, 2, f) #V4 is best, V5 is 2nd best

# best logScore overall
# sort(-colSums(logScores)) #V5 is best again
sort(colSums(logScores_M)) #V5 is best
# sort(-colSums(CRPScores))
# colSums(logScores_SM) #V5 is best

table(apply(logScores_M,1,which.max))
# table(apply(logScores_M,1,which.min))

# table(apply(CRPScores,1,which.max))
# table(apply(CRPScores,1,which.min))

# # median of logScores # not as useful 
# apply(logScores, 2, median) #V6 is best; doesn't penalize its shorter tails as much
# apply(logScores_M, 2, median) #V5 is best, V4 is 2nd best
# apply(logScores, 2, median) #V5 is best, V4 is 2nd best




# evaluate basin predictive densities with log score (by region)

logS_files <- list.files("~/NAM-Model-Validation/csv/scores/basinwide/by_region/", full.names = T, pattern = "storm[1-6].csv")
CRPS_files <- list.files("~/NAM-Model-Validation/csv/scores/basinwide/by_region/", full.names = T, pattern = "CRPS")
logS_marf <- list.files("~/NAM-Model-Validation/csv/scores/basinwide/by_region/alt", full.names = T, pattern = ".csv")
logS_sum_marf <- list.files("~/NAM-Model-Validation/csv/scores/basinwide/by_region/altsum", full.names = T, pattern = ".csv")

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

# # mean of logScores
# sort(-apply(logScores, 2, f)) #V5 is best
# sort(apply(logScores_M, 2, f)) #V4 is best, V5 is 2nd best
# sort(-apply(CRPScores, 2, f)) #V4 is best, V5 is 2nd best
# # apply(logScores_SM, 2, f) #V4 is best, V5 is 2nd best

# best logScore overall
# sort(-colSums(logScores)) #V5 is best again
sort(colSums(logScores_M[logScores_M[,10]=="A",c(-10,-11)])) #V5 is best
sort(colSums(logScores_M[logScores_M[,10]=="F",c(-10,-11)])) #V5 is best
sort(colSums(logScores_M[logScores_M[,10]=="G",c(-10,-11)])) #V5 is best
# sort(-colSums(CRPScores))
# colSums(logScores_SM) #V5 is best
table(apply(logScores_M,1,which.max))
table(apply(logScores_M[logScores_M[,10]=="A",c(-10,-11)],1,which.max))
table(apply(logScores_M[logScores_M[,10]=="F",c(-10,-11)],1,which.max))
table(apply(logScores_M[logScores_M[,10]=="G",c(-10,-11)],1,which.max))
table(apply(logScores_M[logScores_M[,11]=="N",c(-10,-11)],1,which.max))
table(apply(logScores_M[logScores_M[,11]=="C",c(-10,-11)],1,which.max))
# table(apply(logScores_M,1,which.min))

# best logScore overall
# sort(-colSums(logScores)) #V5 is best again
sort(colSums(logScores_M[logScores_M[,11]=="N",c(-10,-11)])) #V5 is best
sort(colSums(logScores_M[logScores_M[,11]=="C",c(-10,-11)])) #V5 is best
# sort(-colSums(CRPScores))
# colSums(logScores_SM) #V5 is best

table(apply(logScores_M[logScores_M[,11]=="N",c(-10,-11)],1,which.max))
table(apply(logScores_M[logScores_M[,11]=="C",c(-10,-11)],1,which.max))
# table(apply(logScores_M,1,which.min))


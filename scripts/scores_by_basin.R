# evaluate basin predictive densities with log score (overall)

logS_files <- list.files("csv/scores", full.names = T, pattern = "storm[1-6].csv")
CRPS_files <- list.files("csv/scores", full.names = T, pattern = "CRPS")
logS_marf <- list.files("csv/scores/alt", full.names = T, pattern = ".csv")
logS_sum_marf <- list.files("csv/scores/altsum", full.names = T, pattern = ".csv")

logScores <- do.call(rbind, lapply(logS_files, read.csv, row.names = 1))
logScores <- logScores[complete.cases(logScores$X.1),]
CRPScores <- do.call(rbind, lapply(CRPS_files, read.csv, row.names = 1))
CRPScores <- CRPScores[complete.cases(CRPScores$X.1),]
logScores_M <- do.call(rbind, lapply(logS_marf, read.csv, row.names = 1))
logScores_M <- logScores_M[complete.cases(logScores_M$X.1),]
logScores_SM <- do.call(rbind, lapply(logS_sum_marf, read.csv, row.names = 1))
logScores_SM <- logScores_SM[complete.cases(logScores_SM$X.1),]

# find the mean wrt the original (NOT log) scale
f <- function(X){log(mean(exp(X - max(X)))) + max(X)}
f2 <- function(X){log(median(exp(X - max(X)))) + max(X)}

# # mean of logScores
# sort(-apply(logScores, 2, f)) #V5 is best
# sort(apply(logScores_M, 2, f)) #V4 is best, V5 is 2nd best
# sort(-apply(CRPScores, 2, f)) #V4 is best, V5 is 2nd best
# # apply(logScores_SM, 2, f) #V4 is best, V5 is 2nd best

# best logScore overall
sort(colSums(logScores_M[,-(10:11)])) #Model 2 is best overall

table(apply(logScores_M[,-c(10:11)],1,which.max)) #Model 2 best by basin

# evaluate basin predictive densities with log score (by region)

# best logScore by region (ATL, FL, or GULF)
sort(colSums(logScores_M[logScores_M[,10]=="A",c(-10,-11)])) #Models 1&2 are best
sort(colSums(logScores_M[logScores_M[,10]=="F",c(-10,-11)])) #Model 1 is best
sort(colSums(logScores_M[logScores_M[,10]=="G",c(-10,-11)])) #Model 2 is best
# best logScore by region (coastal vs noncoastal)
sort(colSums(logScores_M[logScores_M[,11]=="N",c(-10,-11)])) #Model 2 is best
sort(colSums(logScores_M[logScores_M[,11]=="C",c(-10,-11)])) #Model 2 is best

# counts by basin for which score was best, by region
# Model 2 is most often selected as best, except for FL region (then, Model 1)
table(apply(logScores_M,1,which.max))
table(apply(logScores_M[logScores_M[,10]=="A",c(-10,-11)],1,which.max))
table(apply(logScores_M[logScores_M[,10]=="F",c(-10,-11)],1,which.max))
table(apply(logScores_M[logScores_M[,10]=="G",c(-10,-11)],1,which.max))
table(apply(logScores_M[logScores_M[,11]=="N",c(-10,-11)],1,which.max))
table(apply(logScores_M[logScores_M[,11]=="C",c(-10,-11)],1,which.max))

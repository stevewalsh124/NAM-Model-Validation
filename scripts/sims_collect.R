# Gather the results from the 

seed <- 6 #"1orig"

estvalz <- do.call(rbind, lapply(list.files(path = paste0("~/NAM-Model-Validation/csv/simsMLEout/seed",seed,"/"),
                                            pattern = "sim_", full.names = T), read.csv))[,c(6,7,10)]

# The 36th INLA run when seed=1 crashes, so change to this in that instance
estvals <- as.data.frame(t(do.call(cbind, lapply(list.files(paste0("~/NAM-Model-Validation/csv/simsMLEout/seed",seed,"/"),
                                                            pattern = "simlik", full.names = T), 
                                                 function(x){read.csv(x, row.names = 1)})))[,c(1,2,4)])
colnames(estvals) <- colnames(estvalz)
rownames(estvals) <- 1:47

# inlavals <- do.call(rbind, lapply(list.files(path = paste0("~/NAM-Model-Validation/csv/simsMLEout/seed",seed,"/"),pattern = "sim_",
#                                             full.names = T), read.csv))[,2:4]

# if(seed==1|seed=="1orig"){inlavals <- rbind(inlavals[1:26,],NA,inlavals[27:46,])} #sim 36 fails on INLA when seed=1

# INLA bad
# par(mfrow=c(1,3))
# for (i in 1:3) {
#   plot(exp(theta_i_sim[-c(31,36),i]), if(i==3){sqrt(8)/inlavals[-c(31,36),i]}else{inlavals[-c(31,36),i]}, 
#        main = paste(i, "INLA bad")); abline(0,1)
#   print(cor(exp(theta_i_sim[-c(31,36),i]), inlavals[-c(31,36),i]))
# }

if(nchar(seed)==5){load(paste0("~/NAM-Model-Validation/RData/RDatasim0/theta_i_sim",substr(seed, 0, 1)))
} else {load(paste0("~/NAM-Model-Validation/RData/RDatasim0/theta_i_sim",seed))}

#MLE vs true
par(mfrow=c(1,3))
for (i in 1:3) {
  plot(exp(theta_i_sim[,i]), estvals[,i], main = i); abline(0,1)
  print(cor(exp(theta_i_sim[,i]), estvals[,i]))
}
  
hess_sim <- list()
data_files <- list.files(paste0("~/NAM-Model-Validation/RData/RDatasim0/seed",seed),
                         pattern = "hess", full.names = T)
for (pe in 1:47) {
  load(data_files[pe])
  hess_sim[[pe]] <- this_hess
}

rm(list=setdiff(ls(),c("hess_sim", "estvals", "theta_i_sim", "seed")))

print(paste("seed used was", seed))

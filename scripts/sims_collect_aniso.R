# Collect anisotropic sim results
for (i in 1:9) {
  
set <- i
truths <- read.csv(paste0("~/NAM-Model-Validation/csv/aniso/set",set,"/set",set,"truths.csv"))
box <- 20

truestart <- T

# csvs for aniso sims on [0,box]^2 
if(truestart){
  csvs <- list.files(paste0("~/NAM-Model-Validation/csv/aniso/set",set), 
                     pattern = "truestart", full.names = T)
} else {
  csvs <- list.files(paste0("~/NAM-Model-Validation/csv/aniso/set",set), 
                     pattern = paste0("box",box), full.names = T)
  csvs <- grep(csvs, pattern='truestart', invert=TRUE, value=TRUE)
}

tables <- lapply(csvs, read.csv, header = TRUE)
MLEs <- do.call(rbind , tables)[,-1]
pdf(paste0("~/NAM-Model-Validation/pdf/aniso/collect/aniso_sim_results_box", 
           box,"_set",set,if(truestart){"truestart"},".pdf"))

par(mfrow=c(3,2))
for (i in 1:6) {
  hist(MLEs[,i],  main = paste(truths[i,1],": ", round(truths[i,2],3)))
  abline(v=truths[i,2], lwd=2, col="blue")
}

hist(MLEs[,3]*MLEs[,4], main=paste(expression(sigma2/phi),":", round(truths[3,2]/truths[4,2],3)))
abline(v=truths[3,2]/truths[4,2],col="blue", lwd=2)

# #which(MLEs$sigma2 > 2)
# #which(MLEs$phi > 2)      # same set as above
# #which(MLEs$maj.min < 2)  # disjoint as above two
# 
# 
# ## These have reasonable sigma^2 and phi values
# MLEs_good <- MLEs[-which(MLEs$sigma2 > 2),]
# par(mfrow=c(3,2))
# for (i in 1:6) {
#   hist(MLEs_good[,i], main = paste(truths[i,1],": good sig2, phi"))
#   abline(v=truths[i,2], lwd=2, col="blue")
# }
# 
# 
# MLEs_good <- MLEs[which(MLEs$maj.min > 2),]
# par(mfrow=c(3,2))
# for (i in 1:6) {
#   hist(MLEs_good[,i], main = paste(truths[i,1],": good aniso maj.min, angle"))
#   abline(v=truths[i,2], lwd=2, col="blue")
# }
# 
# # csvs for fixed nugget at 0 (F)
# csvs_F <- list.files("~/NAM-Model-Validation/csv/MLEestimates_sqrt_buffer_ngb_PWmean/",
#                      pattern = "fixnug", full.names = T)
# tables_F <- lapply(csvs_F, read.csv, header = TRUE)
# MLEs_F <- do.call(rbind , tables_F)

dev.off()
}

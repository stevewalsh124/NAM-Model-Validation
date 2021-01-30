# Collect anisotropic sim results

# true values
tau2 <- 0.01
sigma2 <- 0.8
phi <- 0.4
sigvec <- c(1,3)
theta <- pi/4

truths <- c(beta=0, tau2=tau2, sigma2=sigma2, phi=phi, theta=theta, maj.min=sigvec[2]/sigvec[1])

box <- 20

# csvs for aniso sims on [0,box]^2
csvs_a4 <- list.files("~/NAM-Model-Validation/csv/aniso", pattern = paste0("box",box), full.names = T)
tables_a4 <- lapply(csvs_a4, read.csv, header = TRUE)
MLEs_a4 <- do.call(rbind , tables_a4)[,-1]
pdf(paste0("~/NAM-Model-Validation/pdf/aniso_sim_results_box", box, ".pdf"))

par(mfrow=c(3,2))
for (i in 1:6) {
  hist(MLEs_a4[,i],  main = paste(names(truths)[i],": all results"))
  abline(v=truths[i], lwd=2, col="blue")
}

#which(MLEs_a4$sigma2 > 2)
#which(MLEs_a4$phi > 2)      # same set as above
#which(MLEs_a4$maj.min < 2)  # disjoint as above two


## These have reasonable sigma^2 and phi values
MLEs_good <- MLEs_a4[-which(MLEs_a4$sigma2 > 2),]
par(mfrow=c(3,2))
for (i in 1:6) {
  hist(MLEs_good[,i], main = paste(names(truths)[i],": good sig2, phi"))
  abline(v=truths[i], lwd=2, col="blue")
}


MLEs_good <- MLEs_a4[which(MLEs_a4$maj.min > 2),]
par(mfrow=c(3,2))
for (i in 1:6) {
  hist(MLEs_good[,i], main = paste(names(truths)[i],": good aniso maj.min, angle"))
  abline(v=truths[i], lwd=2, col="blue")
}

# csvs for fixed nugget at 0 (F)
csvs_F <- list.files("~/NAM-Model-Validation/csv/MLEestimates_sqrt_buffer_ngb_PWmean/",
                     pattern = "fixnug", full.names = T)
tables_F <- lapply(csvs_F, read.csv, header = TRUE)
MLEs_F <- do.call(rbind , tables_F)

dev.off()
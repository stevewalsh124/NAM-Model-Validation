# look at sim coverages for MLEs vs Gibbs

# Do you want results related to sims where the PWmean is subtracted?
genmean <- T
fixedTruth <- F

if(genmean){
  MLE_files <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/MLEcovers", 
                          full.names = T, pattern = "genmean")
  MLE_files <- MLE_files[!grepl("bad", MLE_files)]
  MLE_files <- MLE_files[1:30]
  # Gibbs_files <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/Gibbscoverages", 
  #                           full.names = T, pattern = "genmean")
  PWtrue_files <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/MLEcovers",
                             full.names = T, pattern = "truePW")
} else {
  MLE_all <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/MLEcovers", 
                        full.names = T)
  MLE_sub <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/MLEcovers", 
                        full.names = T, pattern = "genmean")
  MLE_files <- setdiff(MLE_all, MLE_sub)
  Gibbs_all <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/Gibbscoverages", 
                          full.names = T)
  Gibbs_sub <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/Gibbscoverages", 
                          full.names = T, pattern = "genmean")
  Gibbs_files <- setdiff(Gibbs_all, Gibbs_sub)
}


MLEcoverages <- PWTcoverages <- Gibbscoverages <- c()

for (i in 1:length(MLE_files)) {
  MLEcoverages <- rbind(MLEcoverages, colMeans(read.csv(MLE_files[i])[,-1]))
  PWTcoverages <- rbind(PWTcoverages, colMeans(read.csv(PWtrue_files[i])[,-1]))
  # Gibbscoverages <- rbind(Gibbscoverages, read.csv(Gibbs_files[i])$x)
}

colnames(MLEcoverages) <- c("trutheta1","trutheta2")
colMeans(MLEcoverages, na.rm = T)

colnames(PWTcoverages) <- c("trutheta1","trutheta2")
colMeans(PWTcoverages, na.rm = T)

colnames(Gibbscoverages) <- c("theta1hat","theta2hat","trutheta1","trutheta2","truBcomps","truSTcomps")
colMeans(Gibbscoverages, na.rm = T)


## look at sim results to check approx 95% coverage for theta1 and theta2
# Truth fixed across sims, or the truth is normally distributed by landfall location (A, F, G)
if(fixedTruth){
  tru_sigma2 <- 4
  tru_phi <- 1.5
  
  tru_theta1 <- log(tru_sigma2/tru_phi)
  tru_theta2 <- log(tru_sigma2)
} else {
  trueB <- matrix(c(1.037, 1.001, 0.258, 0.244, 0.013, 0.252), 2, 3)
  trueST<- matrix(c(0.238, 0.054, 0.054, 0.220), 2, 2)
  
  truethetaA <- trueB[1:2]
  truethetaF <- trueB[1:2] + trueB[3:4]
  truethetaG <- trueB[1:2] + trueB[5:6]
}

if(genmean){
  cover_files <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/MLEcovers", 
                            pattern = "genmean", full.names = T)
  theta_files <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/thetas", 
                            full.names = T, pattern = "genmean")
  sd_files <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/thetas_sds", 
                         full.names = T, pattern = "genmean")
} else {
  cover_all <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/MLEcovers", 
                          full.names = T)
  cover_sub <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/MLEcovers", 
                          pattern = "genmean", full.names = T)
  cover_files <- setdiff(cover_all, cover_sub)
  theta_all <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/thetas", 
                          full.names = T)
  theta_sub <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/thetas", 
                          pattern = "genmean", full.names = T)
  theta_files <- setdiff(theta_all, theta_sub)
  sd_all <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/thetas_sds", 
                       full.names = T)
  sd_sub <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/thetas_sds", 
                       pattern = "genmean", full.names = T)
  sd_files <- setdiff(sd_all, sd_sub)
}



# number of theta estimates, sd estimates, and coverage files should be the same
if(!all(length(sd_files) == c(length(theta_files), length(cover_files)))) stop({"file count mismatch"})
m <- length(theta_files)
sim_means <- theta_means <- sd_means <- matrix(NA, m, 3)

for (i in 1:m) {
  sim_means[i,] <- colMeans(read.csv(cover_files[i]))
  theta_means[i,] <- colMeans(read.csv(theta_files[i]))
  sd_means[i,] <- colMeans(read.csv(sd_files[i]))
}

colMeans(sim_means, na.rm = T)
colMeans(theta_means, na.rm = T)
colMeans(sd_means, na.rm = T)



theta1hats <- theta2hats <- theta1sds <- theta2sds <- c()#matrix(NA, 47*m2, 2) 
for (i in 1:m) {
  theta1hats <- c(theta1hats, read.csv(theta_files[i])[,2])
  theta2hats <- c(theta2hats, read.csv(theta_files[i])[,3])
  theta1sds <- c(theta1sds, read.csv(sd_files[i])[,2])
  theta2sds <- c(theta2sds, read.csv(sd_files[i])[,3])
  
  # thetahats[(47*(i-1)+1):(47*i),] <- read.csv(theta_files[i])[,2:3]
}

phihats <- exp(theta2hats - theta1hats)
sigma2hats <- exp(theta2hats)

th_cor <- round(cor(theta1hats, theta2hats), 3)
orig_cor <- round(cor(phihats, sigma2hats), 3)

# Colors ideal for color-blind; full palette in NMV_plotting.R
skyblue <- rgb(86/255,180/255,233/255,1)
yellow <- rgb(240/255,228/255,66/255,1)
verm <- rgb(213/255, 94/255, 0,1)
my_cols <- c(skyblue,yellow,verm)

png(paste0("~/NAM-Model-Validation/png/myMLEsim_results",if(genmean){"_genmean"},".png"), 
    width = 3000, height = 2000, res=350)
par(mfrow=c(2,3))
hist(theta1hats,
     main = expression(theta[1]), 
     xlab = expression(theta[1]),
     prob = T)
if(fixedTruth) {
  abline(v = tru_theta1, col="blue", lwd=2)
} else {
  abline(v = truethetaA[1], col= skyblue, lwd=1.5, lty = 4)
  abline(v = truethetaF[1], col= yellow, lwd=1.5, lty = 2)
  abline(v = truethetaG[1], col= verm, lwd=1.5, lty = 3)
}
hist(theta2hats,
     main = expression(theta[2]), 
     xlab = expression(theta[2]),
     prob = T)
if(fixedTruth){
  abline(v = tru_theta2, col="blue", lwd=2)
} else {
  abline(v = truethetaA[2], col= skyblue, lwd=1.5, lty = 4)
  abline(v = truethetaF[2], col= yellow, lwd=1.5, lty = 2)
  abline(v = truethetaG[2], col= verm, lwd=1.5, lty = 3)
}

plot(theta1hats, theta2hats, xlab = expression(theta[1]), ylab = expression(theta[2]),
     main = bquote("cor("~theta[1]~","~theta[2]~")"== .(th_cor)))
if(fixedTruth){
  points(tru_theta1, tru_theta2, pch=18, col="blue", cex=2)
} else {
  points(truethetaA[1], truethetaA[2], pch=18, col= skyblue, cex=1.5)
  points(truethetaF[1], truethetaF[2], pch=18, col= yellow, cex=1.5)
  points(truethetaG[1], truethetaG[2], pch=18, col= verm, cex=1.5)
}

hist(phihats, main = expression(phi), xlab = expression(phi), prob = T)
if(fixedTruth){
  abline(v = tru_phi, col="blue", lwd=2)
} else {
  abline(v = exp(truethetaA[2]-truethetaA[1]), col= skyblue, lwd = 1.5, lty = 4)
  abline(v = exp(truethetaF[2]-truethetaF[1]), col= yellow, lwd = 1.5, lty = 2)
  abline(v = exp(truethetaG[2]-truethetaG[1]), col= verm, lwd = 1.5, lty = 3)
}

hist(sigma2hats, main = expression(sigma^2), xlab = expression(sigma^2), prob = T)
if(fixedTruth){
  abline(v = tru_sigma2, col="blue", lwd=2)
} else {
  abline(v = exp(truethetaA[2]), col= skyblue, lwd = 1.5, lty = 4)
  abline(v = exp(truethetaF[2]), col= yellow, lwd = 1.5, lty = 2)
  abline(v = exp(truethetaG[2]), col= verm, lwd = 1.5, lty = 3)
}

plot(phihats, sigma2hats, xlab = expression(phi), 
     ylab = expression(sigma^2), main = bquote("cor("~phi~","~sigma^2~")"== .(orig_cor)))
if(fixedTruth){
  points(tru_phi, tru_sigma2, pch=18, col="blue", cex=2)
} else {
  points(exp(truethetaA[2] - truethetaA[1]), exp(truethetaA[2]), pch=18, col= skyblue, cex=1.5)
  points(exp(truethetaF[2] - truethetaF[1]), exp(truethetaF[2]), pch=18, col= yellow, cex=1.5)
  points(exp(truethetaG[2] - truethetaG[1]), exp(truethetaG[2]), pch=18, col= verm, cex=1.5)
} 
dev.off()

hist(theta1sds)
hist(theta2sds)

plot(theta1hats, theta1sds)
cor(theta1hats, theta1sds)

plot(theta2hats, theta2sds)
cor(theta2hats, theta2sds)

hist(theta1sds/theta2sds)
hist(theta2sds/theta1sds)

theta1_evals <- theta2_evals <- matrix(NA, 47, m)
for (i in 1:m) {
  theta1_evals[,i] <- read.csv(cover_files[i])[,2]
  theta2_evals[,i] <- read.csv(cover_files[i])[,3]
}

par(mfrow=c(1,2))
hist(rowMeans(theta1_evals))
hist(rowMeans(theta2_evals))

pixels <- c()
files <- list.files("~/NAM-Model-Validation/csv/error_df_sqrt/subtractPWmeanF", full.names = T)
for (i in 1:length(files)) pixels[i] <- nrow(read.csv(files[i]))

plot(pixels, rowMeans(theta1_evals))
plot(pixels, rowMeans(theta2_evals))

cor(pixels, rowMeans(theta1_evals))	
cor(pixels, rowMeans(theta2_evals, na.rm = T))


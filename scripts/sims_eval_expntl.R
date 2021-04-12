## look at sim results to check approx 95% coverage for theta1 and theta2

# sims 1-50 had these values fixed, with no variation
# sims 51-100 used emp_B and emp_Sigma_theta to vary the true values (still pretty close to true tho)
tru_sigma2 <- 4
tru_phi <- 1.5

tru_theta1 <- log(tru_sigma2/tru_phi)
tru_theta2 <- log(tru_sigma2)

# note the [51-100] at the end of these
cover_files <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers", pattern = "seed", full.names = T)[51:100]
theta_files <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/thetas/", full.names = T)[51:100]
sd_files <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers/thetas_sds/", full.names = T)[51:100]


m <- length(cover_files)
sim_means <- matrix(NA, m, 3)
for (i in 1:m) sim_means[i,] <- colMeans(read.csv(cover_files[i]))

colMeans(sim_means)

m2 <- length(theta_files)
theta_means <- matrix(NA, m2, 3)
for (i in 1:m) theta_means[i,] <- colMeans(read.csv(theta_files[i]))

m3 <- length(sd_files)
sd_means <- matrix(NA, m3, 3)
for (i in 1:m) sd_means[i,] <- colMeans(read.csv(sd_files[i]))


theta1hats <- theta2hats <- theta1sds <- theta2sds <- c()#matrix(NA, 47*m2, 2) 
for (i in 1:m2) {
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

# png("~/NAM-Model-Validation/png/myMLEsim_results.png",width = 3000, height = 2000, res=350)
par(mfrow=c(2,3))
hist(theta1hats,
     main = expression(theta[1]), 
     xlab = expression(theta[1]),
     prob = T); abline(v = tru_theta1, col="blue", lwd=2)
hist(theta2hats,
     main = expression(theta[2]), 
     xlab = expression(theta[2]),
     prob = T); abline(v = tru_theta2, col="blue", lwd=2)

plot(theta1hats, theta2hats, xlab = expression(theta[1]), ylab = expression(theta[2]),
     main = bquote("cor("~theta[1]~","~theta[2]~")"== .(th_cor)))
points(tru_theta1, tru_theta2, pch=18, col="blue", cex=2)

hist(phihats, main = expression(phi), xlab = expression(phi), prob = T)
abline(v = tru_phi, col="blue", lwd=2)
hist(sigma2hats, main = expression(sigma^2), xlab = expression(sigma^2), prob = T)
abline(v = tru_sigma2, col="blue", lwd=2)

plot(phihats, sigma2hats, xlab = expression(phi), 
     ylab = expression(sigma^2), main = bquote("cor("~phi~","~sigma^2~")"== .(orig_cor)))
points(tru_phi, tru_sigma2, pch=18, col="blue", cex=2)
# dev.off()

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
cor(pixels, rowMeans(theta2_evals))


## look at sim results to check approx 95% coverage for theta1 and theta2

cover_files <- list.files("~/NAM-Model-Validation/csv/myMLEsimcovers", pattern = "seed", full.names = T)

m <- length(cover_files)
sim_means <- matrix(NA, m, 3)
for (i in 1:m) sim_means[i,] <- colMeans(read.csv(cover_files[i]))

colMeans(sim_means)

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

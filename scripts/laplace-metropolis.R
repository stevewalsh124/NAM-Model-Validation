# Laplace-Metropolis estimator

library(MCMCpack)

# need theta_hat and hessians from GibbsSamplerHurrSE.R
load("RData/Gibbs_sqrt.RData")

# Sigma_theta and Sigma_theta^{-1}
St <- cov(theta_hat)
sti <- solve(St)

# Gather the L_i's in equation 7 from Lewis/Raftery 1997
Li <- c()
EEs <- matrix(NA, nrow(theta_hat), ncol(theta_hat))
VVs <- matrix(NA, nrow(theta_hat), ncol(theta_hat)*2)

# integrating out the random effects, theta_i, i = 1:47
# still a bivariate normal
# plots show the resulting density of f(theta hat | B, Sigma_theta)
# the point on the image is the theta_hat
pdf("pdf/Laplace-Metropolis.pdf")
par(mfrow = c(2,2))
for (s in 1:nrow(theta_hat)) {
  print(s)
  sii <- hessians[[s]]
  
  VVs[s,] <- VV <- solve(sii - sii%*%solve(sii + sti)%*%sii)
  EEs[s,] <- EE <- VV %*% sii%*%solve(sii + sti)%*%sti%*%matrix(emp_B,2,3)%*%matrix(x[s,])
  
  x1seq <- x2seq <- seq(-0.5,2.5,by=0.1)
  lkhds <- matrix(NA, length(x1seq), length(x2seq))
  
  for (i in 1:length(x1seq)) {
    for (j in 1:length(x2seq)) {
      lkhds[i,j] <- dmvn(x = c(x1seq[i], x2seq[j]), mu = t(EE), Sigma = VV, log=T)
    }
  }
  
  Li[s] <- dmvn(x = theta_hat[s,], mu = t(EE), Sigma = VV, log=T)
  
  image(x1seq, x2seq, lkhds, main = paste("storm #", s),
        xlab = expression(theta[1]), ylab = expression(theta[2]))
  points(theta_hat[s,1], theta_hat[s,2], col="green")
}
dev.off()

# obtain posterior mode and covariance matrix for MCMC output
int_mcmc <- cbind(B_burn, Sigma_burn[,-3])
my_modes <- c()

for(i in 1:ncol(int_mcmc)){
  hist(int_mcmc[,i], main = paste("col #", i), prob = T)
  dens <- density(int_mcmc[,i])
  lines(dens)
  my_modes[i] <- dens$x[which(dens$y == max(dens$y))]
}

my_modes
my_cov <- cov(int_mcmc)

Sigma_tilde <- matrix(my_modes[c(7,8,8,9)],2,2)

# get Laplace-Metropolis estimator
# f(theta^*) is the prior, B is propto 1 so I ignore this and only do IW part
LM_c <- P/2*log(2*pi) + 0.5*log(det(my_cov)) + log(diwish(Sigma_tilde, v = v0, v0*cov(theta_hat)))  + sum(Li)
LM_c

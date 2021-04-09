# Laplace-Metropolis estimator

library(MCMCpack)

# need theta_hat and hessians from GibbsSamplerHurrSE.R
load("RData/Gibbs_sqrt.RData")

plot.it <- F

# Sigma_theta and Sigma_theta^{-1}
St <- cov(theta_hat)

###############
#  Model one  #
###############

# Gather the L_i's in equation 7 from Lewis/Raftery 1997
Li <- c()
EEs <- matrix(NA, nrow(theta_hat), ncol(theta_hat))
VVs <- matrix(NA, nrow(theta_hat), ncol(theta_hat)*2)

# integrating out the random effects, theta_i, i = 1:47
# still a bivariate normal
# plots show the resulting density of f(theta hat | B, Sigma_theta)
# the point on the image is the theta_hat
if(plot.it) {pdf("pdf/Laplace-Metropolis.pdf")}
par(mfrow = c(3,2))
for (s in 1:nrow(theta_hat)) {
  
  VVs[s,] <- VV <- solve(hessians[[s]]) + St
  EEs[s,] <- EE <- matrix(emp_B,2,3)%*%matrix(x[s,])
  
  Li[s] <- dmvn(x = theta_hat[s,], mu = t(EE), Sigma = VV, log=T)
  
  if(plot.it){
  x1seq <- x2seq <- seq(-0.5,2.5,by=0.2)
  lkhds <- matrix(NA, length(x1seq), length(x2seq))
  for (i in 1:length(x1seq)) {
    for (j in 1:length(x2seq)) {
      lkhds[i,j] <- dmvn(x = c(x1seq[i], x2seq[j]), mu = t(EE), Sigma = VV, log=T)
    }
  }
  image(x1seq, x2seq, lkhds, main = paste("storm #", s),
        xlab = expression(theta[1]), ylab = expression(theta[2]))
  points(theta_hat[s,1], theta_hat[s,2], col="green")
  }
}

# obtain posterior mode and covariance matrix for MCMC output
int_mcmc1 <- cbind(B_burn, log(Sigma_burn[,1]),Sigma_burn[,2],log(Sigma_burn[,4]))
my_modes1 <- c()

par(mfrow=c(3,3))
for(i in 1:ncol(int_mcmc1)){
  hist(int_mcmc1[,i], main = paste("model 1, col #", i), prob = T)
  dens <- density(int_mcmc1[,i])
  lines(dens)
  my_modes1[i] <- dens$x[which(dens$y == max(dens$y))]
}
# dev.off()

my_modes1
H_star1 <- cov(int_mcmc1)

# Sigma_tilde <- matrix(my_modes[c(7,8,8,9)],2,2)
Sigma_tilde <- matrix(c(exp(my_modes1[7]),my_modes1[c(8,8)],exp(my_modes1[9])),2,2)

# get Laplace-Metropolis estimator
# f(theta^*) is the prior, B is propto 1 so I ignore this and only do IW part
LM_c <- P/2*log(2*pi) + 0.5*log(det(H_star1)) + 
  log(diwish(W = Sigma_tilde, v = v0, S = v0*cov(theta_hat))) + log(Sigma_tilde[1,1]) + log(Sigma_tilde[2,2])  + 
  sum(Li)
LM_c





###############
#  Model two  #
###############

# Gather the L_i's in equation 7 from Lewis/Raftery 1997
Li2 <- c()
mu_theta <- theta_bar

# integrating out the random effects, theta_i, i = 1:47
# still a bivariate normal
# plots show the resulting density of f(theta hat | B, Sigma_theta)
# the point on the image is the theta_hat
par(mfrow = c(3,2))
for (s in 1:nrow(theta_hat)) {
  
  VV <- solve(hessians[[s]]) + St
  EE <- mu_theta
  
  Li2[s] <- dmvn(x = theta_hat[s,], mu = t(EE), Sigma = VV, log=T)
  
  if(plot.it){
  for (i in 1:length(x1seq)) {
    for (j in 1:length(x2seq)) {
      lkhds[i,j] <- dmvn(x = c(x1seq[i], x2seq[j]), mu = t(EE), Sigma = VV, log=T)
    }
  }

    image(x1seq, x2seq, lkhds, main = paste("storm #", s),
          xlab = expression(theta[1]), ylab = expression(theta[2]))
    points(theta_hat[s,1], theta_hat[s,2], col="green")
  }
}

# obtain posterior mode and covariance matrix for MCMC output
# make common mean mu_theta output based on weighted average from B output
emp_mu_theta <-  (9*(B_burn[,1:2])+21*(B_burn[,1:2]+B_burn[,3:4])+17*(B_burn[,1:2]+B_burn[,5:6]))/47
int_mcmc2 <- cbind(emp_mu_theta, log(Sigma_burn[,1]),Sigma_burn[,2],log(Sigma_burn[,4]))
my_modes2 <- c()

par(mfrow=c(3,2))
for(i in 1:ncol(int_mcmc2)){
  hist(int_mcmc2[,i], main = paste("model 2, col #", i), prob = T)
  dens <- density(int_mcmc2[,i])
  lines(dens)
  my_modes2[i] <- dens$x[which(dens$y == max(dens$y))]
}

my_modes2
H_star2 <- cov(int_mcmc2)

# Same Sigma_tilde as Model 1
Sigma_tilde2 <- matrix(c(exp(my_modes2[3]),my_modes2[c(4,4)],exp(my_modes2[5])),2,2)
all.equal(Sigma_tilde2, Sigma_tilde)

# get Laplace-Metropolis estimator
# f(theta^*) is the prior, B is propto 1 so I ignore this and only do IW part
LM_c2 <- P/2*log(2*pi) + 0.5*log(det(H_star2)) + log(diwish(Sigma_tilde, v = v0, v0*cov(theta_hat))) + 
  log(Sigma_tilde[1,1]) + log(Sigma_tilde[2,2])   + sum(Li2)
LM_c2

# Approximate Bayes Factor; this is about 98, showing evidence for a common mean mu_theta
exp(LM_c2 - LM_c)



#################
#  Model three  #
#################

# Gather the L_i's in equation 7 from Lewis/Raftery 1997
Li3 <- c()

# integrating out the random effects, theta_i, i = 1:47
# still a bivariate normal
# plots show the resulting density of f(theta hat | B, Sigma_theta)
# the point on the image is the theta_hat
par(mfrow = c(3,2))
for (s in 1:nrow(theta_hat)) {

  VV <- solve(hessians[[s]])
  EE <- mu_theta
  
  Li3[s] <- dmvn(x = theta_hat[s,], mu = t(EE), Sigma = VV, log=T)
  
  if(plot.it){
  for (i in 1:length(x1seq)) {
    for (j in 1:length(x2seq)) {
      lkhds[i,j] <- dmvn(x = c(x1seq[i], x2seq[j]), mu = t(EE), Sigma = VV, log=T)
    }
  }
    image(x1seq, x2seq, lkhds, main = paste("storm #", s),
          xlab = expression(theta[1]), ylab = expression(theta[2]))
    points(theta_hat[s,1], theta_hat[s,2], col="green")
  }
}

# obtain posterior mode and covariance matrix for MCMC output
emp_mu_theta <-  (9*(B_burn[,1:2])+21*(B_burn[,1:2]+B_burn[,3:4])+17*(B_burn[,1:2]+B_burn[,5:6]))/47
int_mcmc3 <- cbind(emp_mu_theta)
my_modes3 <- c()

par(mfrow=c(1,2))
for(i in 1:ncol(int_mcmc3)){
  hist(int_mcmc3[,i], main = paste("model 3, col #", i), prob = T)
  dens <- density(int_mcmc3[,i])
  lines(dens)
  my_modes3[i] <- dens$x[which(dens$y == max(dens$y))]
}

my_modes3
H_star3 <- cov(int_mcmc3)

# # Sigma_tilde <- matrix(my_modes[c(7,8,8,9)],2,2)
# Sigma_tilde <- matrix(c(exp(my_modes3[3]),my_modes3[c(4,4)],exp(my_modes3[5])),2,2)

# get Laplace-Metropolis estimator
# f(theta^*) is the prior, B is propto 1 so I ignore this and only do IW part
LM_c3 <- P/2*log(2*pi) + 0.5*log(det(H_star3)) + sum(Li3)
LM_c3

# Approximate Bayes Factor; this is about exp(13000); clearly mu_theta should have variability
(LM_c - LM_c3)
exp(LM_c - LM_c3)

LM_c; LM_c2; LM_c3
exp(LM_c2 - LM_c)

if(plot.it){dev.off()}

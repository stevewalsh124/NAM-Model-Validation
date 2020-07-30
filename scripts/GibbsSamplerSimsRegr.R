# Gibbs sampler for Simulated Hurricane Error Fields
# Steve Walsh Feb 2020

# pdf("~/NAM-Model-Validation/pdf/Gibbs/GibbsSamplerHurrplotsSimsRegr.pdf")
# remove(list=ls())
load("~/NAM-Model-Validation/RData/par_optim_allbut43.RData")
rm(list=setdiff(ls(), c("par_optim")))

set.seed(2)

library(MCMCpack)    #riwish (inverse wishart draws)
library(LaplacesDemon) #rmvn, rmatrixnorm (multivariate normal draws)

as.square <- function(mat){matrix(mat, nrow=sqrt(length(mat)),ncol=sqrt(length(mat)))}
rounder <- 9 #number of decimal places to round matrices so solve doesn't induce asymmetry
nsims <- 46 #number of theta_i, i \in 1:nsims; number of simulated storms

all_avgHess <- T # if FALSE, must have nsims <= 46
if(!all_avgHess & nsims > 46) nsims <- 46

#turns vector of length n^2 into nxn square matrix 
as.square <- function(mat){matrix(mat, nrow=sqrt(length(mat)),ncol=sqrt(length(mat)))}
#testy <- matrix(rnorm(10000),100,100); all.equal(testy, as.square(as.vector(testy)))

# These are the actual hurricane estimates
theta_hat <- matrix(NA, nrow= length(par_optim[[1]]), ncol=length(par_optim[[1]][[1]]))
N <- nrow(theta_hat)
P <- ncol(theta_hat)
for(i in 1:N){ theta_hat[i,] <- par_optim[[1]][[i]] }
colnames(theta_hat) <- c("MLEsigma2","MLEphi","MLEkappa") 

theta_bar <- apply(theta_hat, 2, mean)
hessians <- par_optim[[2]]

true_Sigma_theta <- cov(theta_hat)

# read in location and intensity info
loc_int <- read.csv("~/NAM-Model-Validation/csv/storm_levels_and_locs.csv", row.names = 1)[-43,]
colnames(loc_int) <- c("int","loc")
nA <- table(loc_int$loc)[1]
nF <- table(loc_int$loc)[2]
nG <- table(loc_int$loc)[3]

# Model matrix for locations
x <- with(loc_int, model.matrix(~ loc))
sum_x_xt <- matrix(0,3,3)
for(i in 1:dim(x)[1]){sum_x_xt <- sum_x_xt + x[i,]%*%t(x[i,])}
V <- round(solve(sum_x_xt), rounder)

sum_th_xt_true <- matrix(0, P, P)
for(j in 1:N){sum_th_xt_true <- sum_th_xt_true+ theta_hat[j,]%*%t(x[j,])}

true_M<- sum_th_xt_true %*% round(solve(sum_x_xt), rounder)

# # obtain mean vector by location
theta_A <- lm(theta_hat ~ loc_int$loc)$coefficients[1,]
theta_F <- lm(theta_hat ~ loc_int$loc)$coefficients[2,]
theta_G <- lm(theta_hat ~ loc_int$loc)$coefficients[3,]
# # cbind(theta_A, theta_F, theta_G) and true_M are approximately equal

# obtain covariance matrix by location
# can't use. not enough degrees of freedom
# use common sigma_theta across locations: cov(theta_hat)
# cov_A <- cov(filter(as.data.frame(theta_hat), loc_int$loc=="ATL"))
# cov_F <- cov(filter(as.data.frame(theta_hat), loc_int$loc=="FL"))
# cov_G <- cov(filter(as.data.frame(theta_hat), loc_int$loc=="GULF"))


# Generate simulated data based on the hurricane data:
# Simulate theta_i, i=1:nsims, from a normal

avgHess <- Reduce("+", hessians) / length(hessians)
# is.positive.definite(avgHess) is TRUE, not for solve(avgHess)

theta_i_sim  <- theta_i_hat_sim  <- matrix(NA, nrow = nsims, ncol = 3)
colnames(theta_i_sim) <- colnames(theta_i_hat_sim) <- c("MLEsigma2","MLEphi","MLEkappa") 
for (i in 1:nsims) {
  # theta_i iid \sim N(\bar\mu_\theta, \Sigma_\theta) 
  theta_i_sim[i,] <- rmvn(1, t(true_M%*%x[i,]), true_Sigma_theta)
  if(all_avgHess){
    theta_i_hat_sim[i,] <- rmvn(1, theta_i_sim[i,], round(solve(avgHess), rounder)) #all covmats the same
    hessians[[i]] <- avgHess
  }else{
    theta_i_hat_sim[i,] <- rmvn(1, theta_i_sim[i,], round(solve(hessians[[i]]), rounder))} #or each covmat matches with H_i^-1
}

par(mfrow=c(2,3))
for (i in 1:3) {hist(theta_i_sim[,i], main=paste("theta_",i))}
for (i in 1:3) {hist(theta_i_hat_sim[,i], main=paste("theta_",i,"hat"))}
for (i in 1:length(hessians)) {if(!is.positive.definite(hessians[[i]])){print("one of the hessians is not pos def")}}

# Change the actual data to the simulated data:
# hessians changed to avgHess 9 lines above # hessians[[i]] <- avgHess
theta_hat <- theta_i_hat_sim



# For full conditionals, will need:
#
# N=47
# \bar\theta and \hat\theta_i
# H_i (hessian for theta_i)
#
# B, x_i
# \Sigma_\theta (and its inverse)
# v_0, S_0

# set.seed(Sys.time()) #seed seems to stay the same from the loaded data


#Gibbs sampler for B, theta_i and Sigma_theta
iters <- 5100
burn  <- 100

# mu_theta     <- matrix(NA, nrow = iters, ncol = P)
Sigma_theta  <- matrix(NA, nrow = iters, ncol = P^2)
B            <- matrix(NA, nrow = iters, ncol = P^2)
theta        <- list()

for (j in 1:N) {
  theta[[j]]     <- matrix(NA, nrow = iters, ncol = P)
  theta[[j]][1,] <- c(100,100,100)
}

# Initial values for Markov chains
# mu_theta[1,]    <- c(50,50,50)            #theta_bar#EMPIRICAL BAYES
Sigma_theta[1,] <- diag(rep(100, P))      #cov(theta_hat)#
B[1,]           <- diag(rep(100, P))  
#IW hyperparameters
v0 <- 4                                   #increase for EMP BAYES/debugging
S0 <- diag(c(.01,.01,.01))                #v0*cov(theta_hat)#EMPIRICAL BAYES

for (i in 2:iters) {
  #i <- 2; j <- 1
  if(i %% 500 ==0) print(i)
  
  # Update B
  sum_th_xt <- matrix(0, P, P)
  for(j in 1:N){sum_th_xt <- sum_th_xt+ theta[[j]][i-1,]%*%t(x[j,])}
  M <- sum_th_xt %*% V
  B[i,] <- rmatrixnorm(M = M, U = as.square(Sigma_theta[i-1,]), V = V)
  
  # Update Sigma_theta
  L <- matrix(0, ncol = P, nrow = P)
  for(j in 1:N){
    # Variance of theta_i
    Vtheta_i <- solve(hessians[[j]] + solve(as.square(Sigma_theta[i-1,])))
    Vtheta_i <- round(Vtheta_i, rounder)
    
    # Expectation of theta_i
    Etheta_i <- Vtheta_i %*% ((hessians[[j]]%*%theta_hat[j,]) + (solve(as.square(Sigma_theta[i-1,]))%*%(as.square(B[i,])%*%x[j,])))
    
    theta[[j]][i,]  <- rmvn(1, t(Etheta_i), Vtheta_i)
    
    # L is the component of the scale matrix from the likelihood for the inverse wishart draw
    L <- L + round((theta[[j]][i,]-(as.square(B[i,])%*%x[j,]))%*%t(theta[[j]][i,]-(as.square(B[i,])%*%x[j,])), rounder)
  }
  Sigma_theta[i,] <- riwish(N + v0, L + S0)
  # Sigma_theta[i,] <- LaplacesDemon::rsiw(nu = v0, S = L + diag(3), mu = c(0,0,0), delta = c(1,1,1))
  Sigma_theta[i,] <- round(Sigma_theta[i,], rounder)
}

theta_cols <- c("sigma2","phi","kappa") 

# Histograms and trace plots for B
par(mfrow=c(1,3))
for(i in 1: P^2) {
  hist(B[(burn+1):iters,i], main = paste("B",i))
  abline(v=as.vector(true_M)[i], col="blue")
  abline(v=apply(B[burn:iters,],2,mean)[i],col="green")
}
for(i in 1: P^2) {
  plot(B[(burn+1):iters,i], main = paste("B",i), type="l")
  abline(h=as.vector(true_M)[i], col="blue")
  abline(h=apply(B[burn:iters,],2,mean)[i],col="green")
}

# Trace plots for inidividual storms' parameters
for (j in c(1:5,(nsims-5):nsims)) {
  for(i in 1: P) {
    plot(theta[[j]][(burn+1):iters,i], main = paste("storm",j,loc_int$loc[j],theta_cols[i]), type="l")
         # ylim=c(apply(theta_hat, 2, min)[i] - 0.2, apply(theta_hat, 2, max)[i] + 0.2))
    # abline(h=theta_bar[i], col="red")
    abline(h=apply(theta[[j]][(burn+1):iters,], 2, mean)[i], col="green")
    abline(h=theta_i_sim[j,i], col="blue")
    abline(h=theta_i_hat_sim[j,i], col="orange")
  }
}

# Variance histograms
# true_sigma_theta: generates the sims (from the actual storm data) (blue)
# cov(theta_hat):   the covariance estimated from the sims (orange)
par(mfrow=c(1,3))
# Variance histograms
for (i in c(1,5,9)) {
  hist(as.vector(Sigma_theta[(burn+1):iters,i]), main = paste("Variance of theta_[",i,"]"))
  abline(v=as.vector(cov(theta_hat))[i],col="orange")
  abline(v=as.vector(true_Sigma_theta)[i],col="blue")
  abline(v=apply(Sigma_theta[(burn+1):iters,],2,mean)[i], col="green")
}

# Covariance histograms
for (i in 1:P){ 
  for(j in 1:P){
    # if(i!=j) next
    hist((Sigma_theta[(burn+1):iters, (i-1)*P+j]), main = paste("Covariance of theta_[",i,",",j,"]"))
    abline(v=as.vector(cov(theta_hat)[(i-1)*P+j]),col="orange")
    abline(v=as.vector(true_Sigma_theta[(i-1)*P+j]),col="blue")
    abline(v=apply(Sigma_theta[(burn+1):iters,],2,mean)[(i-1)*P+j], col="green")
  }
}

#Correlations between parts of regression coefficient matrix B
for (i in 1:P^2) { 
  j <- i + 1
  while(j <= P^2) {
    plot(B[(burn+1):iters,i], B[(burn+1):iters,j],
         main=paste("cor of B",i,",",j,":", 
                    round(cor(B[(burn+1):iters,i], B[(burn+1):iters,j]),3)))
    print(round(cor(B[(burn+1):iters,i], B[(burn+1):iters,j]),3))
    j <- j + 1
  }}

var_results <- rbind(as.vector(S0), apply(Sigma_theta[(burn+1):iters,], 2,mean), as.vector(cov(theta_hat)))
rownames(var_results) <- c("S0", "Sigma_theta", "cov(theta_hat)")
var_results

rnorm(3) #to check for seed being the same



###### Evaluate coverage of the Gibbs sampler

B_burn <- B[(burn+1):iters,]
Sigma_burn <- Sigma_theta[(burn+1):iters,]
theta_burn <- list()
for (i in 1:N) {theta_burn[[i]] <- theta[[i]][(burn+1):iters,]}

# "true" data for sim study
# theta_i_sim
# theta_i_hat_sim
# true_M
# true_Sigma_theta

# estimates from MCMC
# B_burn      ; creates emp_B, emp_B_LB, emp_B_UB
# Sigma_burn  ; creates emp_Sigma_theta, emp_Sigma_theta_LB, emp_Sigma_theta_UB
# theta_burn  ; creates emp_thetaMED, emp_thetaLB, emp_thetaUB

emp_thetaLB <- emp_thetaMED <- emp_thetaUB <- matrix(NA, nrow = N, ncol = P)
for (i in 1:N) {
  # emp_theta[i,] <- apply(theta[[i]], 2, mean)
  emp_thetaLB[i,] <- apply(theta_burn[[i]], 2, function(x){quantile(x, 0.025)})
  emp_thetaMED[i,]<- apply(theta_burn[[i]], 2, function(x){quantile(x, 0.5)})
  emp_thetaUB[i,] <- apply(theta_burn[[i]], 2, function(x){quantile(x, 0.975)})
}

apply(theta_i_sim < emp_thetaUB & theta_i_sim > emp_thetaLB, 2, sum)/N

# Collect the "bad thetas": one of the true elements of theta_i 
# isn't contained in the 95% credible interval
bad_thetas <- which(apply(theta_i_sim < emp_thetaUB & theta_i_sim > emp_thetaLB, 1, all)==F)
# cbind(emp_thetaLB[bad_thetas,], theta_i_sim[bad_thetas,], emp_thetaUB[bad_thetas,])

par(mfrow=c(1,3))
for (j in head(bad_thetas,3)) {
  for(i in 1: P) {
    plot(theta[[j]][(burn+1):iters,i], 
         main = paste("storm",j,loc_int$loc[j],theta_cols[i]), type="l")#,
         # ylim=c(apply(theta_hat, 2, min)[i] + 0.2, apply(theta_hat, 2, max)[i] - 0.2))
    abline(h=apply(theta[[j]][(burn+1):iters,], 2, mean)[i], col="green")
    abline(h=theta_i_sim[j,i], col="blue")
    abline(h=theta_i_hat_sim[j,i], col="orange")
    abline(h=theta_bar[i], col="red")
  }
}

# Check if each element of true B and Sigma_theta is contained in the credible intervals
emp_B    <- apply(B_burn, 2, function(x){quantile(x,0.5)})
emp_B_UB <- apply(B_burn, 2, function(x){quantile(x,0.975)})
emp_B_LB <- apply(B_burn, 2, function(x){quantile(x,0.025)})

true_M < emp_B_UB & true_M > emp_B_LB

emp_Sigma_theta    <- apply(Sigma_burn, 2, function(x){quantile(x,0.5)})
emp_Sigma_theta_LB <- apply(Sigma_burn, 2, function(x){quantile(x,0.025)})
emp_Sigma_theta_UB <- apply(Sigma_burn, 2, function(x){quantile(x,0.975)})

true_Sigma_theta < emp_Sigma_theta_UB & true_Sigma_theta > emp_Sigma_theta_LB




# Are the theta_i's close to the estimated \hat\theta^GS_i?
# emp_thetaMED - theta_i_hat_sim
apply(emp_thetaMED - theta_i_sim, 2, function(X)max(abs(X)))
apply(emp_thetaMED - theta_i_sim, 2, mean)

emp_B - true_M
emp_Sigma_theta - true_Sigma_theta

for (i in 1:P) {
  hist(emp_thetaMED[,i] - theta_i_sim[,i], main = paste("theta_hat_GS - theta,", colnames(theta_i_sim)[i]))
  abline(v=0, col="green")
}

rnorm(1)
# dev.off()

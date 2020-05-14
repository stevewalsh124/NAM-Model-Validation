# Gibbs sampler for Simulated Hurricane Error Fields
# Steve Walsh Feb 2020

# run storms_collect.R first

# pdf("~/NAM-Model-Validation/pdf/Gibbs/GibbsSamplerHurrRegr_EMPBAYESIW_NA_GHG.pdf")
# remove(list=ls())
# load("NAM-Model-Validation/RData/par_optim_allbut43.RData")
# rm(list=setdiff(ls(), c("par_optim")))

# set.seed(489)

library(MCMCpack)    #riwish (inverse wishart draws)
library(LaplacesDemon) #rmvn, rmatrixnorm (multivariate normal draws)

rounder <- 9 #number of decimal places to round matrices so solve doesn't induce asymmetry
# nsims <- 46 #number of theta_i, i \in 1:nsims; number of simulated storms

# all_avgHess <- T # if FALSE, must have nsims <= 46
# if(!all_avgHess & nsims > 46) nsims <- 46

#turns vector of length n^2 into nxn square matrix 
as.square <- function(mat){matrix(mat, nrow=sqrt(length(mat)),ncol=sqrt(length(mat)))}
#testy <- matrix(rnorm(10000),100,100); all.equal(testy, as.square(as.vector(testy)))

# These are the actual hurricane estimates
theta_hat <- all_storm_res[,c(5,6,9)]#matrix(NA, nrow= length(par_optim[[1]]), ncol=length(par_optim[[1]][[1]]))
N <- nrow(theta_hat)
P <- ncol(theta_hat)
# for(i in 1:N){ theta_hat[i,] <- par_optim[[1]][[i]] }
# colnames(theta_hat) <- c("MLEsigma2","MLEphi","MLEkappa") 

hessians <- hess_opt#par_optim[[2]]
for (i in 1:length(hessians)) {if(!is.positive.definite(hessians[[i]])){print("one of the hessians is not pos def")}}

###############DELTA METHOD###################
##############################################
for (i in 1:length(hessians)) {
  G <- diag(theta_hat[i,])
  hessians[[i]] <- G%*%hessians[[i]]%*%t(G)
  print(paste("symm",isSymmetric(hessians[[i]])))
  
  hessians[[i]] <- (hessians[[i]]+t(hessians[[i]]))/2
  
  print(paste("pos def",is.positive.definite(hessians[[i]])))
}

theta_hat <- log(theta_hat)

##############################################
##############################################

theta_bar <- apply(theta_hat, 2, mean)

# true_Sigma_theta <- cov(theta_hat)

# read in location and intensity info
loc_int <- read.csv("~/NAM-Model-Validation/csv/storm_levels_and_locs.csv", row.names = 1)[avail,]#[-43,]
colnames(loc_int) <- c("int","loc")

# nA <- table(loc_int$loc)[1]
# nF <- table(loc_int$loc)[2]
# nG <- table(loc_int$loc)[3]

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

# theta_i_sim  <- theta_i_hat_sim  <- matrix(NA, nrow = nsims, ncol = 3)
# colnames(theta_i_sim) <- colnames(theta_i_hat_sim) <- c("MLEsigma2","MLEphi","MLEkappa") 
# for (i in 1:nsims) {
#   # theta_i iid \sim N(\bar\mu_\theta, \Sigma_\theta) 
#   theta_i_sim[i,] <- rmvn(1, t(true_M%*%x[i,]), true_Sigma_theta)
#   if(all_avgHess){
#     theta_i_hat_sim[i,] <- rmvn(1, theta_i_sim[i,], round(solve(avgHess), rounder)) #all covmats the same
#     hessians[[i]] <- avgHess
#   }else{
#     theta_i_hat_sim[i,] <- rmvn(1, theta_i_sim[i,], round(solve(hessians[[i]]), rounder))} #or each covmat H_i^-1
# }

par(mfrow=c(1,3))
for (i in 1:3) {hist(theta_hat[,i], main=paste("theta_",i))}
# for (i in 1:3) {hist(log(theta_hat[,i]), main=paste("log theta_",i))}
# for (i in 1:3) {hist(theta_i_hat_sim[,i], main=paste("theta_",i,"hat"))}
for (i in 1:length(hessians)) {if(!is.positive.definite(hessians[[i]])){print("one of the hessians is not pos def")}}

# Change the actual data to the simulated data:
# hessians changed to avgHess 9 lines above # hessians[[i]] <- avgHess
# theta_hat <- theta_i_hat_sim



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
iters <- 10020
burn  <- 20

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
S0 <- v0*cov(theta_hat)#diag(c(.01,.01,.01))                #EMPIRICAL BAYES

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
    Etheta_i <- Vtheta_i %*% ((hessians[[j]]%*%theta_hat[j,]) + 
                                (solve(as.square(Sigma_theta[i-1,]))%*%(as.square(B[i,])%*%x[j,])))
    
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
  abline(v=as.vector(true_M)[i], lwd=2, col="blue")
  abline(v=apply(B[burn:iters,],2,mean)[i], lwd=2, col="green")
}
for(i in 1: P^2) {
  plot(B[(burn+1):iters,i], main = paste("B",i), type="l")
  abline(h=as.vector(true_M)[i], lwd=2, col="blue")
  abline(h=apply(B[burn:iters,],2,mean)[i], lwd=2, col="green")
}

# Trace plots for inidividual storms' parameters
for (j in 1:N){#c(1:5,(nsims-5):nsims)) {
  for(i in 1: P) {
    plot(theta[[j]][(burn+1):iters,i], main = paste("storm",j,loc_int$loc[j],theta_cols[i]), type="l")
    # ylim=c(apply(theta_hat, 2, min)[i] - 0.2, apply(theta_hat, 2, max)[i] + 0.2))
    # abline(h=theta_bar[i], lwd=2, col="red")
    abline(h=apply(theta[[j]][(burn+1):iters,], 2, mean)[i], lwd=2, col="green")
    # abline(h=theta_i_sim[j,i], lwd=2, col="blue")
    abline(h=theta_hat[j,i], lwd=2, col="blue")
  }
}

# Variance histograms
# true_sigma_theta: generates the sims (from the actual storm data) (blue)
# cov(theta_hat):   the covariance estimated from the sims (orange)
par(mfrow=c(1,3))
# Variance histograms
for (i in c(1,5,9)) {
  hist(as.vector(Sigma_theta[(burn+1):iters,i]), main = paste("Variance of theta_[",i,"]"))
  abline(v=as.vector(cov(theta_hat))[i], lwd=2, col="blue")
  # abline(v=as.vector(true_Sigma_theta)[i], lwd=2, col="blue")
  abline(v=apply(Sigma_theta[(burn+1):iters,],2,mean)[i], lwd=2, col="green")
}

# Covariance histograms
for (i in 1:P){ 
  for(j in 1:P){
    # if(i!=j) next
    hist((Sigma_theta[(burn+1):iters, (i-1)*P+j]), main = paste("Covariance of theta_[",i,",",j,"]"))
    abline(v=as.vector(cov(theta_hat)[(i-1)*P+j]), lwd=2, col="blue")
    # abline(v=as.vector(true_Sigma_theta[(i-1)*P+j]), lwd=2, col="blue")
    abline(v=apply(Sigma_theta[(burn+1):iters,],2,mean)[(i-1)*P+j], lwd=2, col="green")
  }
}

# #Correlations between parts of regression coefficient matrix B
# for (i in 1:P^2) { 
#   j <- i + 1
#   while(j <= P^2) {
#     plot(B[(burn+1):iters,i], B[(burn+1):iters,j],
#          main=paste("cor of B",i,",",j,":", 
#                     round(cor(B[(burn+1):iters,i], B[(burn+1):iters,j]),3)))
#     print(round(cor(B[(burn+1):iters,i], B[(burn+1):iters,j]),3))
#     j <- j + 1
#   }
# }


var_results <- rbind(as.vector(S0), apply(Sigma_theta[(burn+1):iters,], 2,mean), as.vector(cov(theta_hat)))
rownames(var_results) <- c("S0", "Sigma_theta", "cov(theta_hat)")
var_results

rnorm(3) #to check for seed being the same... -1.6760883 -0.8578595  1.2997342



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

apply(theta_hat < emp_thetaUB & theta_hat > emp_thetaLB, 2, sum)/N

# Collect the "bad thetas": one of the true elements of theta_i 
# isn't contained in the 95% credible interval
bad_thetas <- which(apply(theta_hat < emp_thetaUB & theta_hat > emp_thetaLB, 1, all)==F)
# cbind(emp_thetaLB[bad_thetas,], theta_i_sim[bad_thetas,], emp_thetaUB[bad_thetas,])

par(mfrow=c(1,3))
for (j in head(bad_thetas,3)) {
  for(i in 1: P) {
    plot(theta[[j]][(burn+1):iters,i], 
         main = paste("storm",j,loc_int$loc[j],theta_cols[i]), type="l")#,
    # ylim=c(apply(theta_hat, 2, min)[i] + 0.2, apply(theta_hat, 2, max)[i] - 0.2))
    abline(h=apply(theta[[j]][(burn+1):iters,], 2, mean)[i], lwd=2, col="green")
    abline(h=theta_hat[j,i], lwd=2, col="blue")
    # abline(h=theta_i_hat_sim[j,i], lwd=2, col="orange")
    abline(h=theta_bar[i], lwd=2, col="red")
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

cov(theta_hat) < emp_Sigma_theta_UB & cov(theta_hat) > emp_Sigma_theta_LB




# Are the theta_i's close to the estimated \hat\theta^GS_i?
# emp_thetaMED - theta_i_hat_sim
apply(emp_thetaMED - theta_hat, 2, function(X)max(abs(X)))
apply(emp_thetaMED - theta_hat, 2, mean)

emp_B - true_M
emp_Sigma_theta - cov(theta_hat)

for (i in 1:P) {
  hist(emp_thetaMED[,i] - theta_hat[,i], main = paste("theta_hat_GS - theta,", colnames(theta_hat)[i]))
  abline(v=0, lwd=2, col="green")
}


## Make plots comparing theta_i posterior density to 
## theta_i_hat from MLE with corresponding Hessian^-1 variance
# pdf("NAM-Model-Validation/pdf/Gibbs/compare_postthetai_thetaihat_col.pdf")
par(mfrow=c(1,1))
# make it easier to see each of the densities in their corresponding plots
smush <- c(.06,.055,.009)
for (i in 1:P) {
  plot(0, 0, col = "white", xlab = "", ylab = "", 
       xlim=c(min(theta_hat[,i], emp_thetaMED[,i]),
              max(theta_hat[,i], emp_thetaMED[,i])),
       ylim=c(0,2), yaxt='n',
       main = bquote(theta[.(i)]~"bottom and"~hat(theta)[.(i)]~"top"))
  # mult_seg <- data.frame(x0 = c(0.7, 0.2, - 0.9, 0.5, - 0.2),    # Create data frame with line-values
  #                        y0 = c(0.5, 0.5, 0.6, - 0.3, 0.4),
  #                        x1 = c(0, 0.2, 0.4, - 0.3, - 0.6),
  #                        y1 = c(- 0.1, 0.3, - 0.6, - 0.8, 0.9))
  
  segments(x0 = emp_thetaMED[,i],                            # Draw multiple lines
           y0 = 0,
           x1 = theta_hat[,i],
           y1 = 1, col = loc_int$loc)
  points(x=emp_thetaMED[,i], y=rep(0, length(emp_thetaMED[,i])), col= loc_int$loc)
  points(x=theta_hat[,i],    y=rep(1, length(emp_thetaMED[,i])), col= loc_int$loc)
  for (j in 1:N) {
    #top row, theta hats
    xseq <- seq(theta_hat[j,i]-3*sqrt(solve(hessians[[j]])[i,i]),
                theta_hat[j,i]+3*sqrt(solve(hessians[[j]])[i,i]),
                length.out = 1000)
    lines(xseq,smush[i]*dnorm(xseq,
                              theta_hat[j,i],
                              sqrt(solve(hessians[[j]])[i,i]))+1, col=loc_int$loc[j],
          lty = as.integer(loc_int$loc)[j])
    
    #bottom row, thetas from MCMC
    xseq <- seq(min(theta_burn[[j]][,i]),
                max(theta_burn[[j]][,i]), 
                length.out = 1000)
    lines(xseq, smush[i]*dnorm(xseq, 
                               mean(theta_burn[[j]][,i]), 
                               sd(theta_burn[[j]][,i])), col=loc_int$loc[j],
          lty = as.integer(loc_int$loc)[j])
  }
  legend(max(theta_hat[,i], emp_thetaMED[,i])-
           (max(theta_hat[,i], emp_thetaMED[,i])-min(theta_hat[,i], emp_thetaMED[,i]))*.4, 
         2, 
         legend=c(levels(loc_int$loc)[unique(as.integer(loc_int$loc))]),
         col=unique(as.integer(loc_int$loc)), lty=unique(as.integer(loc_int$loc)), cex=0.8)
}


sig2var <- phivar <- kapvar <- matrix(NA,N,2)
for (i in 1:P) {
  for (j in 1:N) {
    if(i==1) sig2var[j,] <- c(sqrt(solve(hessians[[j]])[i,i]),sd(theta_burn[[j]][,i]))
    if(i==2) phivar[j,] <- c(sqrt(solve(hessians[[j]])[i,i]),sd(theta_burn[[j]][,i]))
    if(i==3) kapvar[j,] <- c(sqrt(solve(hessians[[j]])[i,i]),sd(theta_burn[[j]][,i]))
  }
  
}

par(mfrow=c(1,3))
plot(sig2var, main = bquote(paste("Variances for "~ hat(sigma)^2~" from "~ bar(H)^{-1},~" and "~sigma[MCMC]^2)),
     xlab = expression(hat(sigma)^2), ylab = expression(sigma^2), col=loc_int$loc,  pch=as.integer(loc_int$loc)-1,
     xlim = range(sig2var), ylim = range(sig2var),
     cex=2, cex.lab=2, cex.axis=2, cex.main=1.5, cex.sub=2)
abline(0,1)
plot(phivar, main = expression(paste("Variances for ",hat(phi)," from ", bar(H)^{-1}," and ", phi[MCMC])),
     xlab = expression(hat(phi)), ylab = expression(phi), col=loc_int$loc, pch=as.integer(loc_int$loc)-1,
     xlim = range(phivar), ylim = range(phivar),
     cex=2, cex.lab=2, cex.axis=2, cex.main=1.5, cex.sub=2)
abline(0,1)
plot(kapvar, main = expression(paste("Variances for ",hat(kappa)," from ", bar(H)^{-1}," and ", kappa[MCMC])),
     xlab = expression(hat(kappa)), ylab = expression(kappa), col=loc_int$loc,  pch=as.integer(loc_int$loc)-1,
     xlim = range(kapvar), ylim = range(kapvar),
     cex=2, cex.lab=2, cex.axis=2, cex.main=1.5, cex.sub=2)
abline(0,1)

# dev.off()


# Comparing the different variance components
# > sqrt(diag(solve(avgHess)))
# MLEsigma2      MLEphi    MLEkappa 
# 0.016262686 0.022689820 0.006525455 
# > sqrt(diag(cov(theta_hat)))
# MLEsigma2    MLEphi  MLEkappa 
# 0.2544399 0.2556946 0.0822819 
# > apply(theta_burn[[1]], 2, sd)
# [1] 0.04515972 0.07133124 0.01014894
# > apply(theta_var_avg, 2, mean)
# [1] 0.05040922 0.07389148 0.01297306

theta_var_avg <- matrix(NA,N,P)
for(i in 1:N) theta_var_avg[i,] <- apply(theta_burn[[i]], 2, sd)
apply(theta_var_avg, 2, mean)

hess_inv <- matrix(NA, N, P)
for (i in 1:N) {
  hess_inv[i,] <- sqrt(diag(solve(hessians[[i]])))
}
apply(hess_inv, 2, mean)


# Posterior Medians for each element of theta_i
burn_meds <- matrix(NA,N,P)
for (i in 1:N) {
  burn_meds[i,] <- apply(theta_burn[[i]], 2, median)
}

# png("NAM-Model-Validation/png/burn_med_hist.png",width = 1400, height=1000)
par(mfrow=c(1,3))
hist(burn_meds[,1], xlab=expression(sigma^2), main=NULL, cex.axis=2, cex.lab=2)
hist(burn_meds[,2], xlab=expression(phi), main=NULL, cex.axis=2, cex.lab=2)
hist(burn_meds[,3], xlab=expression(kappa), main=NULL, cex.axis=2, cex.lab=2)
# dev.off()

hist(theta_hat[,1], xlab=expression(hat(sigma^2)[MLE]), main=NULL, cex.axis=2, cex.lab=2)
hist(theta_hat[,2], xlab=expression(hat(phi)[MLE]), main=NULL, cex.axis=2, cex.lab=2)
hist(theta_hat[,3], xlab=expression(hat(kappa)[MLE]), main=NULL, cex.axis=2, cex.lab=2)

# dev.off()
rnorm(5)

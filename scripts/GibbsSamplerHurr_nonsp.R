# Gibbs sampler for Hurricane Error Fields
# Nonspatial case
# Steve Walsh Feb 2020

# run storms_collect.R first for geoR MLEs, not myMLEs (see below)

# pick data vs posterior pwmean; change pdf and save.image() names

# pdf("pdf/Gibbs/GibbsSamplerHurrRegr_EMPBAYESIW_NA_GHG_sqrt.pdf")
# remove(list=ls())
# load("NAM-Model-Validation/RData/par_optim_allbut43.RData")
# rm(list=setdiff(ls(), c("par_optim")))

# set.seed(489)

suppressMessages(library(MCMCpack))    #riwish (inverse wishart draws)
suppressMessages(library(LaplacesDemon)) #rmvn, rmatrixnorm (multivariate normal draws)

subtractPWmean <- F
rounder <- 9 #number of decimal places to round matrices so solve doesn't induce asymmetry
# nsims <- 46 #number of theta_i, i \in 1:nsims; number of simulated storms

# all_avgHess <- T # if FALSE, must have nsims <= 46
# if(!all_avgHess & nsims > 46) nsims <- 46

#turns vector of length n^2 into nxn square matrix 
as.square <- function(mat){matrix(mat, nrow=sqrt(length(mat)),ncol=sqrt(length(mat)))}
#testy <- matrix(rnorm(10000),100,100); all.equal(testy, as.square(as.vector(testy)))

# These are the actual hurricane estimates
# lambda_hat <- all_storm_res[,c("MLEsigma2","MLEphi")]
stormMLEfiles <- list.files("csv/myMLEresults/nonsp/MLEs/", full.names = T)
myMLEs   <- do.call(rbind, lapply(stormMLEfiles, read.csv, row.names = 1))$x

lambda_hat <- myMLEs

N <- length(lambda_hat) # number of storms, 47
P <- 1 # number of params, log(sigma^2)
R <- 3 #number of landfall locations (ATL, FL, GULF)

# for(i in 1:N){ lambda_hat[i,] <- par_optim[[1]][[i]] }
# colnames(lambda_hat) <- c("MLEsigma2","MLEphi","MLEkappa") 

theta_hat <- log(lambda_hat)

hessians <- list()
hess_lambda_files <- list.files("csv/myMLEresults/nonsp/hessians/", full.names = T)

if(length(hess_lambda_files) != N){stop("number of MLEs != number of Hessians")}
hessians_lam <- do.call(rbind, lapply(hess_lambda_files, read.csv, row.names = 1))$x

# hessians after delta method
hessians <- -hessians_lam*lambda_hat^2

theta_bar <- mean(theta_hat)
theta_wgt <- sum(theta_hat*hessians)/sum(hessians)
# true_Sigma_theta <- cov(theta_hat)

# read in location and intensity info
loc_int <- read.csv("csv/storm_levels_and_locs.csv", row.names = 1)#[avail,]#[-43,]
colnames(loc_int) <- c("int","loc")

# nA <- table(loc_int$loc)[1]
# nF <- table(loc_int$loc)[2]
# nG <- table(loc_int$loc)[3]

# Model matrix for locations
x <- with(loc_int, model.matrix(~ loc))
sum_x_xt <- matrix(0, R, R)
for(i in 1:dim(x)[1]){sum_x_xt <- sum_x_xt + x[i,]%*%t(x[i,])}
V <- (solve(sum_x_xt)+t(solve(sum_x_xt)))/2

sum_thhat_xt <- matrix(0, P, R)
for(j in 1:N){sum_thhat_xt <- sum_thhat_xt+ matrix(theta_hat[j],1,1)%*%t(x[j,])}

hat_M <- sum_thhat_xt %*% ((solve(sum_x_xt)+t(solve(sum_x_xt)))/2)

# # obtain mean vector by location
theta_A <- lm(theta_hat ~ loc_int$loc)$coefficients[1]
theta_F <- lm(theta_hat ~ loc_int$loc)$coefficients[2]
theta_G <- lm(theta_hat ~ loc_int$loc)$coefficients[3]
# # cbind(theta_A, theta_F, theta_G) and hat_M are approximately equal

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
#   theta_i_sim[i,] <- rmvn(1, t(hat_M%*%x[i,]), true_Sigma_theta)
#   if(all_avgHess){
#     theta_i_hat_sim[i,] <- rmvn(1, theta_i_sim[i,], round(solve(avgHess), rounder)) #all covmats the same
#     hessians[[i]] <- avgHess
#   }else{
#     theta_i_hat_sim[i,] <- rmvn(1, theta_i_sim[i,], round(solve(hessians[[i]]), rounder))} #or each covmat H_i^-1
# }

hist(theta_hat, main=paste("theta_",i))
# for (i in 1:3) {hist(log(theta_hat[,i]), main=paste("log theta_",i))}
# for (i in 1:3) {hist(theta_i_hat_sim[,i], main=paste("theta_",i,"hat"))}

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
iters <- 10100
burn  <- 100

# mu_theta     <- matrix(NA, nrow = iters, ncol = P)
Sigma_theta  <- matrix(NA, nrow = iters, ncol = P^2)
B            <- matrix(NA, nrow = iters, ncol = P*R)
theta        <- list()

for (j in 1:N) {
  theta[[j]]     <- matrix(NA, nrow = iters, ncol = P)
  theta[[j]][1] <- 100
}

# Initial values for Markov chains
# mu_theta[1,] <- rep(100, P)
Sigma_theta[1] <- 100 #diag(rep(100, P))
B[1,]           <- matrix(100, P, R)
#IW hyperparameters
v0 <- P + 1
S0 <- v0*var(theta_hat) #EMPIRICAL BAYES

for (i in 2:iters) {
  #i <- 2; j <- 1
  if(i %% 2000 ==0) print(i)
  
  # Update B
  sum_th_xt <- matrix(0, P, R)
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
    Etheta_i <- Vtheta_i %*% ((hessians[[j]]%*%theta_hat[j]) + 
                                (solve(as.square(Sigma_theta[i-1,]))%*%(matrix(B[i,],P,3)%*%x[j,])))
    
    theta[[j]][i,]  <- rmvn(1, t(Etheta_i), Vtheta_i)
    
    # L is the component of the scale matrix from the likelihood for the inverse wishart draw
    L <- L + round((theta[[j]][i,]-(matrix(B[i,],P,3)%*%x[j,]))%*%t(theta[[j]][i,]-(matrix(B[i,],P,3)%*%x[j,])), rounder)
  }
  Sigma_theta[i,] <- riwish(N + v0, L + S0)
  # Sigma_theta[i,] <- LaplacesDemon::rsiw(nu = v0, S = L + diag(3), mu = c(0,0,0), delta = c(1,1,1))
  Sigma_theta[i,] <- round(Sigma_theta[i,], rounder)
}

theta_cols <- c("log(sigma2)") 

# Histograms and trace plots for B
par(mfrow=c(1,R))
for(i in 1:(P*R)) {
  hist(B[(burn+1):iters,i], main = paste("B",i))
  abline(v=as.vector(hat_M)[i], lwd=2, col="blue")
  abline(v=apply(B[burn:iters,],2,mean)[i], lwd=2, col="green")
}
for(i in 1:(P*R)) {
  plot(B[(burn+1):iters,i], main = paste("B",i), type="l")
  abline(h=as.vector(hat_M)[i], lwd=2, col="blue")
  abline(h=apply(B[burn:iters,],2,mean)[i], lwd=2, col="green")
}

# Trace plots for inidividual storms' parameters
for (j in 1:N){#c(1:5,(nsims-5):nsims)) {
  plot(theta[[j]][(burn+1):iters], main = paste("storm",j,loc_int$loc[j],theta_cols), type="l")
  # ylim=c(apply(theta_hat, 2, min)[i] - 0.2, apply(theta_hat, 2, max)[i] + 0.2))
  # abline(h=theta_bar[i], lwd=2, col="red")
  abline(h=mean(theta[[j]][(burn+1):iters]), lwd=2, col="green")
  # abline(h=theta_i_sim[j,i], lwd=2, col="blue")
  abline(h=theta_hat[j], lwd=2, col="blue")
}

# Variance histograms
# true_sigma_theta: generates the sims (from the actual storm data) (blue)
# cov(theta_hat):   the covariance estimated from the sims (orange)
# par(mfrow=c(1,3))
# Variance histograms
par(mfrow=c(2,2))
for (i in 1) { #1,4 are the diagonals of a vectorized 2x2 matrix
  hist(as.vector(Sigma_theta[(burn+1):iters]), main = paste("Variance of theta_[",i,"]"))
  abline(v=as.vector(var(theta_hat)), lwd=2, col="blue")
  # abline(v=as.vector(true_Sigma_theta)[i], lwd=2, col="blue")
  abline(v=mean(Sigma_theta[(burn+1):iters]), lwd=2, col="green")
}


#Correlations between parts of regression coefficient matrix B
for (i in 1:(R-1)) {
  for(j in (i+1):R) {
    plot(B[(burn+1):iters,i], B[(burn+1):iters,j],
         main=paste("cor of B",i,",",j,":",
                    round(cor(B[(burn+1):iters,i], B[(burn+1):iters,j]),3)))
    print(round(cor(B[(burn+1):iters,i], B[(burn+1):iters,j]),3))
  }
}


var_results <- rbind(as.vector(S0), mean(Sigma_theta[(burn+1):iters,]), as.vector(var(theta_hat)))
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
# hat_M
# true_Sigma_theta

# estimates from MCMC
# B_burn      ; creates emp_B, emp_B_LB, emp_B_UB
# Sigma_burn  ; creates emp_Sigma_theta, emp_Sigma_theta_LB, emp_Sigma_theta_UB
# theta_burn  ; creates emp_thetaMED, emp_thetaLB, emp_thetaUB

emp_thetaLB <- emp_thetaMED <- emp_thetaUB <- matrix(NA, nrow = N, ncol = P)
for (i in 1:N) {
  # emp_theta[i,] <- apply(theta[[i]], 2, mean)
  emp_thetaLB[i] <- quantile(theta_burn[[i]], 0.025)
  emp_thetaMED[i]<- quantile(theta_burn[[i]], 0.5)
  emp_thetaUB[i] <- quantile(theta_burn[[i]], 0.975)
}

# This is expected to be 1; thetahat will be in the 95% credible interval
apply(theta_hat < emp_thetaUB & theta_hat > emp_thetaLB, 2, sum)/N

# Collect the "bad thetas": one of the true elements of theta_i 
# isn't contained in the 95% credible interval
bad_thetas <- which(apply(theta_hat < emp_thetaUB & theta_hat > emp_thetaLB, 1, all)==F)
# cbind(emp_thetaLB[bad_thetas,], theta_i_sim[bad_thetas,], emp_thetaUB[bad_thetas,])

par(mfrow=c(1,P))
for (j in head(bad_thetas,3)) {
  for(i in 1:P) {
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

hat_M < emp_B_UB & hat_M > emp_B_LB

emp_Sigma_theta    <- quantile(Sigma_burn,0.5)
emp_Sigma_theta_LB <- quantile(Sigma_burn,0.025)
emp_Sigma_theta_UB <- quantile(Sigma_burn,0.975)

var(theta_hat) < emp_Sigma_theta_UB & var(theta_hat) > emp_Sigma_theta_LB




# Are the theta_i's close to the estimated \hat\theta^GS_i?
# emp_thetaMED - theta_i_hat_sim
apply(emp_thetaMED - theta_hat, 2, function(X)max(abs(X)))
apply(emp_thetaMED - theta_hat, 2, mean)

emp_B - hat_M
emp_Sigma_theta - var(theta_hat)

hist(emp_thetaMED - theta_hat, main = paste("theta_hat_GS - theta,", colnames(theta_hat)[i]))
abline(v=0, lwd=2, col="green")


## Make plots comparing theta_i posterior density to 
## theta_i_hat from MLE with corresponding Hessian^-1 variance
# pdf("NAM-Model-Validation/pdf/Gibbs/compare_postthetai_thetaihat_col.pdf")
par(mfrow=c(1,P))
# make it easier to see each of the densities in their corresponding plots
smush <- .03
for (i in 1:P) {
  plot(0, 0, col = "white", xlab = "", ylab = "", 
       xlim=c(min(theta_hat, emp_thetaMED),
              max(theta_hat, emp_thetaMED)),
       ylim=c(0,2), yaxt='n',
       main = bquote(theta[.(i)]~"bottom and"~hat(theta)[.(i)]~"top"))
  # mult_seg <- data.frame(x0 = c(0.7, 0.2, - 0.9, 0.5, - 0.2),    # Create data frame with line-values
  #                        y0 = c(0.5, 0.5, 0.6, - 0.3, 0.4),
  #                        x1 = c(0, 0.2, 0.4, - 0.3, - 0.6),
  #                        y1 = c(- 0.1, 0.3, - 0.6, - 0.8, 0.9))
  
  segments(x0 = emp_thetaMED,                            # Draw multiple lines
           y0 = 0,
           x1 = theta_hat,
           y1 = 1, col = as.factor(loc_int$loc))
  points(x=emp_thetaMED, y=rep(0, length(emp_thetaMED)), col= as.factor(loc_int$loc))
  points(x=theta_hat,    y=rep(1, length(emp_thetaMED)), col= as.factor(loc_int$loc))
  for (j in 1:N) {
    #top row, theta hats
    xseq <- seq(c(theta_hat[j]-3*sqrt(solve(hessians[[j]]))),
                c(theta_hat[j]+3*sqrt(solve(hessians[[j]]))),
                length.out = 1000)
    lines(xseq,smush[i]*dnorm(xseq,
                              theta_hat[j],
                              sqrt(solve(hessians[[j]])))+1, col=as.factor(loc_int$loc)[j],
          lty = i)
    
    #bottom row, thetas from MCMC
    xseq <- seq(min(theta_burn[[j]]),
                max(theta_burn[[j]]), 
                length.out = 1000)
    lines(xseq, smush[i]*dnorm(xseq, 
                               mean(theta_burn[[j]]), 
                               sd(theta_burn[[j]])), col=as.factor(loc_int$loc)[j],
          lty = i)
  }
  legend("topleft",
         legend=unique((loc_int$loc)),
         col=unique(as.factor(loc_int$loc)), lty=1:3, cex=0.8)
}


theta1sds <- matrix(NA, N, 2)
for (j in 1:N) { theta1sds[j,] <- c(sqrt(solve(hessians[[j]])),sd(theta_burn[[j]])) }

par(mfrow=c(1,P))
plot(theta1sds, main = bquote(paste("SDs for "~ hat(theta)[1]~" from "~ bar(H)^{-1},~" and "~theta[1,MCMC])),
     xlab = expression(hat(theta)[1]), ylab = expression(theta[1]), col=as.factor(loc_int$loc),  pch=as.numeric(as.factor(loc_int$loc))-1,
     xlim = range(theta1sds), ylim = range(theta1sds),
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
for(i in 1:N) theta_var_avg[i,] <- sd(theta_burn[[i]])
apply(theta_var_avg, 2, mean)

hess_inv <- matrix(NA, N, P)
for (i in 1:N) {
  hess_inv[i,] <- sqrt(diag(solve(hessians[[i]])))
}
apply(hess_inv, 2, mean)


# Posterior Medians for each element of theta_i
burn_meds <- c()
for (i in 1:N) burn_meds[i] <- median(theta_burn[[i]])

# png("NAM-Model-Validation/png/burn_med_hist.png",width = 1400, height=1000)
par(mfrow=c(1,P))
hist(burn_meds, xlab=expression(sigma^2), main=NULL, cex.axis=2, cex.lab=2)
# dev.off()

hist(theta_hat, xlab=expression(hat(sigma^2)[MLE]), main=NULL, cex.axis=2, cex.lab=2)
hist(burn_meds - theta_hat, main = "log(sigma2): posterior medians VS MLE")

# dev.off()
rnorm(5)

matrix(emp_B, P, R)
matrix(emp_Sigma_theta, P, P)

save.image(file = paste0("RData/Gibbs_nosp.RData"))

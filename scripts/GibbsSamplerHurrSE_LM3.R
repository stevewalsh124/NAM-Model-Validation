# Gibbs sampler for Hurricane Error Fields
# LM3: Non-hierarchical version for theta
# Steve Walsh Feb 2020

# pdf("pdf/Gibbs/GibbsSamplerHurrRegr_EMPBAYESIW_NA_GHG_sqrt.pdf")
# remove(list=ls())
# load("NAM-Model-Validation/RData/par_optim_allbut43.RData")
# rm(list=setdiff(ls(), c("par_optim")))

# set.seed(489)

if(!dir.exists("RData/myMLE_precs/")){
  dir.create("RData/myMLE_precs/", recursive = T)
}

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
stormMLEfiles <- grep(list.files("csv/myMLEresults/myMLEs", full.names = T, recursive = F),
                      pattern='sim_vals', invert=TRUE, value=TRUE)
ind <- grepl("subtractPWmean", stormMLEfiles)
if(subtractPWmean){
  myMLEfiles <- stormMLEfiles[ind]
} else {
  myMLEfiles <- stormMLEfiles[!ind]
}
myMLEs   <- do.call(rbind, lapply(myMLEfiles, read.csv))

lambda_hat <- cbind(myMLEs$sigs,myMLEs$phis)

N <- nrow(lambda_hat) # number of storms, 47
P <- ncol(lambda_hat) # number of params, theta1 and theta2
R <- 3 #number of landfall locations (ATL, FL, GULF)

# for(i in 1:N){ lambda_hat[i,] <- par_optim[[1]][[i]] }
# colnames(lambda_hat) <- c("MLEsigma2","MLEphi","MLEkappa") 

theta_hat <- cbind(log(lambda_hat[,1]/lambda_hat[,2]), log(lambda_hat[,1]))

hessians <- list()
all_hess_theta_files <- list.files("csv/myMLEresults/pkgthetahessvecs", full.names = T)	
ind <- grepl("subtractPWmean", all_hess_theta_files)	
if(subtractPWmean){	
  hess_theta_files <- all_hess_theta_files[ind]	
} else {	
  hess_theta_files <- all_hess_theta_files[!ind]
}
if(length(hess_theta_files) != N){stop("number of MLEs != number of Hessians")}
for (i in 1:N) {
  hess <- read.csv(hess_theta_files[i], row.names = 1)
  hess_mtx <- as.square(as.numeric(-1*hess))
  hessians[[i]] <- (hess_mtx + t(hess_mtx))/2
  if(!isSymmetric(hessians[[i]])){stop(paste("storm",i,"is not symm,etric"))}
  if(!is.positive.definite(hessians[[i]])){stop(paste("storm",i,"is not pos def"))}
}

theta_bar <- apply(theta_hat, 2, mean)

sumHess <- Reduce("+", hessians)
sumHessI <- (solve(sumHess) + t(solve(sumHess)))/2
# is.positive.definite(avgHess) is TRUE, not for solve(avgHess)

hess_thetai <- list()
for (i in 1:N) {hess_thetai[[i]] <- hessians[[i]] %*% theta_hat[i,]}
sumHess_thetai <- Reduce("+", hess_thetai)


par(mfrow=c(1,P))
for (i in 1:P) {hist(theta_hat[,i], main=paste("theta_",i))}
# for (i in 1:3) {hist(log(theta_hat[,i]), main=paste("log theta_",i))}
# for (i in 1:3) {hist(theta_i_hat_sim[,i], main=paste("theta_",i,"hat"))}

#Gibbs sampler for mu_theta, theta_i 
iters <- 10100
burn  <- 100

mu_theta     <- matrix(NA, nrow = iters, ncol = P)
theta        <- list()

for (j in 1:N) {
  theta[[j]]     <- matrix(NA, nrow = iters, ncol = P)
  theta[[j]][1,] <- c(100,100)
}

# Initial values for Markov chains
mu_theta[1,]    <- rep(100, P)

#IW hyperparameters
v0 <- P + 1
S0 <- v0*cov(theta_hat) #EMPIRICAL BAYES

for (i in 2:iters) {
  #i <- 2; j <- 1
  if(i %% 2000 ==0) print(i)
  
  # Update mu_theta
  mu_theta[i,] <- rmvn(1, t(sumHessI %*% sumHess_thetai), sumHessI)
  
  # Update Sigma_theta
  L <- matrix(0, ncol = P, nrow = P)
  for(j in 1:N){
    # Variance of theta_i
    Vtheta_i <- (solve(hessians[[j]])+t(solve(hessians[[j]])))/2

    # Expectation of theta_i
    Etheta_i <- theta_hat[j,]
    theta[[j]][i,]  <- rmvn(1, t(Etheta_i), Vtheta_i)
  }

}

theta_cols <- c("theta1","theta2") 

mu_burn <- mu_theta[(burn+1):iters,]
theta_burn <- list()
for (i in 1:N) {theta_burn[[i]] <- theta[[i]][(burn+1):iters,]}

# Histograms and trace plots for B
par(mfrow=c(1,P))
for(i in 1:P) {
  hist(mu_burn[,i], main = paste("mu_theta",i))
  abline(v=theta_bar[i], lwd=2, col="blue")
  abline(v=colMeans(mu_burn)[i], lwd=2, col="green")
}
for(i in 1:(P)) {
  plot(mu_burn[,i], main = paste("mu_theta",i), type="l")
  abline(h=theta_bar[i], lwd=2, col="blue")
  abline(h=colMeans(mu_burn)[i], lwd=2, col="green")
}

# Trace plots for inidividual storms' parameters
for (j in 1:N){#c(1:5,(nsims-5):nsims)) {
  for(i in 1:P) {
    plot(theta_burn[[j]][,i], main = paste("storm",j,theta_cols[i]), type="l")
    # ylim=c(apply(theta_hat, 2, min)[i] - 0.2, apply(theta_hat, 2, max)[i] + 0.2))
    # abline(h=theta_bar[i], lwd=2, col="red")
    abline(h=colMeans(theta_burn[[j]])[i], lwd=2, col="green")
    # abline(h=theta_i_sim[j,i], lwd=2, col="blue")
    abline(h=theta_hat[j,i], lwd=2, col="blue")
  }
}



###### Evaluate coverage of the Gibbs sampler

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
  emp_thetaLB[i,] <- apply(theta_burn[[i]], 2, function(x){quantile(x, 0.025)})
  emp_thetaMED[i,]<- apply(theta_burn[[i]], 2, function(x){quantile(x, 0.5)})
  emp_thetaUB[i,] <- apply(theta_burn[[i]], 2, function(x){quantile(x, 0.975)})
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
         main = paste("storm",j,theta_cols[i]), type="l")#,
    # ylim=c(apply(theta_hat, 2, min)[i] + 0.2, apply(theta_hat, 2, max)[i] - 0.2))
    abline(h=apply(theta[[j]][(burn+1):iters,], 2, mean)[i], lwd=2, col="green")
    abline(h=theta_hat[j,i], lwd=2, col="blue")
    # abline(h=theta_i_hat_sim[j,i], lwd=2, col="orange")
    abline(h=theta_bar[i], lwd=2, col="red")
  }
}

# Check if each element of true B and Sigma_theta is contained in the credible intervals
emp_mu    <- apply(mu_burn, 2, function(x){quantile(x,0.5)})
emp_mu_UB <- apply(mu_burn, 2, function(x){quantile(x,0.975)})
emp_mu_LB <- apply(mu_burn, 2, function(x){quantile(x,0.025)})

theta_bar < emp_mu_UB & theta_bar > emp_mu_LB


# Are the theta_i's close to the estimated \hat\theta^GS_i?
# emp_thetaMED - theta_i_hat_sim
apply(emp_thetaMED - theta_hat, 2, function(X)max(abs(X)))
apply(emp_thetaMED - theta_hat, 2, mean)

emp_mu - theta_bar

for (i in 1:P) {
  hist(emp_thetaMED[,i] - theta_hat[,i], main = paste("theta_hat_GS - theta,", colnames(theta_hat)[i]))
  abline(v=0, lwd=2, col="green")
}


## Make plots comparing theta_i posterior density to 
## theta_i_hat from MLE with corresponding Hessian^-1 variance
# pdf("NAM-Model-Validation/pdf/Gibbs/compare_postthetai_thetaihat_col.pdf")
par(mfrow=c(1,P))
# make it easier to see each of the densities in their corresponding plots
smush <- c(.04,.15)
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
           y1 = 1)
  points(x=emp_thetaMED[,i], y=rep(0, length(emp_thetaMED[,i])))
  points(x=theta_hat[,i],    y=rep(1, length(emp_thetaMED[,i])))
  for (j in 1:N) {
    #top row, theta hats
    xseq <- seq(theta_hat[j,i]-3*sqrt(solve(hessians[[j]])[i,i]),
                theta_hat[j,i]+3*sqrt(solve(hessians[[j]])[i,i]),
                length.out = 1000)
    lines(xseq,smush[i]*dnorm(xseq,
                              theta_hat[j,i],
                              sqrt(solve(hessians[[j]])[i,i]))+1, lty = i)
    
    #bottom row, thetas from MCMC
    xseq <- seq(min(theta_burn[[j]][,i]),
                max(theta_burn[[j]][,i]), 
                length.out = 1000)
    lines(xseq, smush[i]*dnorm(xseq, 
                               mean(theta_burn[[j]][,i]), 
                               sd(theta_burn[[j]][,i])),lty = i)
  }
}


theta1sds <- theta2sds <- matrix(NA,N,P)
for (i in 1:P) {
  for (j in 1:N) {
    if(i==1) theta1sds[j,] <- c(sqrt(solve(hessians[[j]])[i,i]),sd(theta_burn[[j]][,i]))
    if(i==2) theta2sds[j,] <- c(sqrt(solve(hessians[[j]])[i,i]),sd(theta_burn[[j]][,i]))
  }
}

par(mfrow=c(1,P))
plot(theta1sds, main = bquote(paste("SDs for "~ hat(theta)[1]~" from "~ bar(H)^{-1},~" and "~theta[1,MCMC])),
     xlab = expression(hat(theta)[1]), ylab = expression(theta[1]), 
     xlim = range(theta1sds), ylim = range(theta1sds),
     cex=2, cex.lab=2, cex.axis=2, cex.main=1.5, cex.sub=2)
abline(0,1)
plot(theta2sds, main = expression(paste("SDs for ",hat(theta)[2]," from ", bar(H)^{-1}," and ", theta[2,MCMC])),
     xlab = expression(hat(theta)[2]), ylab = expression(theta[2]), 
     xlim = range(theta2sds), ylim = range(theta2sds),
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
par(mfrow=c(1,P))
hist(burn_meds[,1], xlab=expression(sigma^2), main=NULL, cex.axis=2, cex.lab=2)
hist(burn_meds[,2], xlab=expression(phi), main=NULL, cex.axis=2, cex.lab=2)
# dev.off()

hist(theta_hat[,1], xlab=expression(hat(sigma^2)[MLE]), main=NULL, cex.axis=2, cex.lab=2)
hist(theta_hat[,2], xlab=expression(hat(phi)[MLE]), main=NULL, cex.axis=2, cex.lab=2)

# dev.off()
rnorm(5)

emp_mu

save.image(file = paste0("RData/Gibbs_sqrt_LM3",
           if(subtractPWmean){"_subtractPWmean"},".RData"))

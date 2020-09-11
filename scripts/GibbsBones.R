# Basic sim to try and establish problem with other Gibbs sampler
# Steve Walsh August 2020

source("~/NAM-Model-Validation/scripts/sims_collect.R")

#load theta_i_sim's from errorfieldGPsims

library(raster)
library(fields)
library(LaplacesDemon)
library(MCMCpack)

nowtime <- function(){gsub(gsub(gsub(Sys.time(),pattern = " ", replacement = ""),
                                pattern="-", replacement=""),pattern=":", replacement="")}

nreps <- 1

iters <- 5100
burn  <- 100

trueBSig <- F

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

hat_cover <- tru_cover <- matrix(NA, nrow = nreps, ncol = 3)

for (nr in 1:nreps) {
  print(paste("this is rep", nr))
  seed <- 1
  # set.seed(seed)
  
  load("~/NAM-Model-Validation/RData/truthB")
  load("~/NAM-Model-Validation/RData/truthSig")
  # load("~/NAM-Model-Validation/RData/hess_sim")
  
  symmtrz <- function(X){as.matrix((X+t(X))/2)}
  as.square <- function(mat){matrix(mat, nrow=sqrt(length(mat)),ncol=sqrt(length(mat)))}
  
  # read in location and intensity info
  loc_int <- read.csv("~/NAM-Model-Validation/csv/storm_levels_and_locs.csv", row.names = 1)
  colnames(loc_int) <- c("int","loc")
  
  # nA <- table(loc_int$loc)[1]; nF <- table(loc_int$loc)[2]; nG <- table(loc_int$loc)[3]
  
  # Model matrix for locations
  x <- with(loc_int, model.matrix(~ loc))
  sum_x_xt <- matrix(0,3,3)
  for(i in 1:dim(x)[1]){sum_x_xt <- sum_x_xt + x[i,]%*%t(x[i,])}
  V <- symmtrz(solve(sum_x_xt))
  
  # grp_B <- cbind(exp(truthB[,1]), exp(truthB[,1] + truthB[,2]), exp(truthB[,1]+truthB[,3]))
  N <- 47; P <- 3
  
  # hes_vec <- cov_vec <- matrix(NA, N, P^2)
  # for (i in 1:N) {
  #   cov_vec[i,] <- c(solve(hess_sim[[i]]))
  #   hes_vec[i,] <- c(hess_sim[[i]])
  # }
  
  # Psi <- symmtrz(as.square(colMeans(cov_vec)))# symmtrz(solve(as.square(colMeans(hes_vec))))#
  
  # theta_hat_sim <- matrix(NA,N,P)
  # for (i in 1:N) {theta_hat_sim[i,] <- rmvn(n=1, mu=theta_i_sim[i,], Sigma=Psi)}
  
  par(mfrow=c(1,3))
  for (i in 1:P) {
    plot(theta_i_sim[,i], log(estvals)[,i], main=paste("sim true vs sim hat", i))
    abline(0,1)
  }
  
  hessians <- hess_sim
  theta_hat <- log(estvals)#theta_hat_sim
  
  ###############DELTA METHOD###################
  ##############################################
  
  # Note this is modifying the inverse of the asymptotic covariance matrix
  # So, when MVDM says \Sigma ==> G \Sigma G^T,
  # This is doing \Sigma^-1 = H = G^-T H G^-1
  # Thus, the G here is really G^-1, since the diags would otherwise be 1/theta_hat[i,]
  
  GHGt <- list()
  for (i in 1:length(hessians)) {
    G <- diag(estvals[i,]) 
    GHGt[[i]] <-symmtrz(G%*%hessians[[i]]%*%t(G))
    # GHGt[[i]] <- symmtrz(solve(Psi))
    
    if(!isSymmetric.matrix(GHGt[[i]])){
      print("one of the hessians is not symmetric")}
    if(!is.positive.definite(GHGt[[i]])){
      print("one of the hessians is not pos def")}
  }
  
  ##############################################
  ##############################################
  
  theta_bar <- apply(theta_hat, 2, mean)
  
  # true_Sigma_theta <- cov(theta_hat)
  
  # sum_th_xt_true <- matrix(0, P, P)
  # for(j in 1:N){sum_th_xt_true <- sum_th_xt_true+ theta_hat[j,]%*%t(x[j,])}
  
  #truthB #sum_th_xt_true %*% round(solve(sum_x_xt), rounder)
  
  # # obtain mean vector by location
  # theta_A <- lm(theta_hat ~ loc_int$loc)$coefficients[1,]
  # theta_F <- lm(theta_hat ~ loc_int$loc)$coefficients[2,]
  # theta_G <- lm(theta_hat ~ loc_int$loc)$coefficients[3,]
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
  avgGHGt <- Reduce("+", GHGt) / length(GHGt)
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
  v0 <- 4.01                                   #increase for EMP BAYES/debugging
  S0 <- v0*cov(theta_hat)#diag(c(.01,.01,.01))                #EMPIRICAL BAYES
  
  for (i in 2:iters) {
    #i <- 2; j <- 1; i <- 3
    if(i %% 500 ==0) print(i)
    
    # Update B
    sum_th_xt <- matrix(0, P, P)
    for(j in 1:N){sum_th_xt <- sum_th_xt+ theta[[j]][i-1,]%*%t(x[j,])}
    M <- sum_th_xt %*% V
    if(trueBSig){ B[i,] <- truthB 
    } else { B[i,] <- rmatrixnorm(M = M, U = as.square(Sigma_theta[i-1,]), V = V) }
    
    # Update Sigma_theta
    L <- matrix(0, ncol = P, nrow = P)
    for(j in 1:N){
      # Variance of theta_i
      Vtheta_i <- solve(GHGt[[j]] + solve(as.square(Sigma_theta[i-1,])))
      Vtheta_i <- symmtrz(Vtheta_i)
      
      # Expectation of theta_i
      Etheta_i <- Vtheta_i %*% ((GHGt[[j]]%*%t(theta_hat[j,])) + 
                                  (solve(as.square(Sigma_theta[i-1,]))%*%(as.square(B[i,])%*%x[j,])))
      
      theta[[j]][i,]  <- rmvn(1, t(Etheta_i), Vtheta_i)
      
      # L is the component of the scale matrix from the likelihood for the inverse wishart draw
      L <- L + symmtrz((theta[[j]][i,]-(as.square(B[i,])%*%x[j,]))%*%t(theta[[j]][i,]-(as.square(B[i,])%*%x[j,])))
    }
    
    if(trueBSig){ Sigma_theta[i,] <- truthSig 
    } else { Sigma_theta[i,] <- symmtrz(as.square(riwish(N + v0, L + S0))) }
    # Sigma_theta[i,] <- symmtrz(LaplacesDemon::rinvwishart(nu = N + v0, S = L + S0))
    # Sigma_theta[i,] <- symmtrz(LaplacesDemon::rsiw(nu = N + v0, S = L + S0, mu = c(0,0,0), delta = c(1,1,1)))
  }
  
  
  # What do we expect? 95% coverage for each of the params composing theta_i's, coverage for generating B and Sigma_theta
  
  B_burn <- B[(burn+1):iters,]
  Sigma_burn <- Sigma_theta[(burn+1):iters,]
  theta_burn <- list()
  for (i in 1:N) {theta_burn[[i]] <- theta[[i]][(burn+1):iters,]}
  
  emp_thetaLB <- emp_thetaMED <- emp_thetaUB <- matrix(NA, nrow = N, ncol = P)
  for (i in 1:N) {
    # emp_theta[i,] <- apply(theta[[i]], 2, mean)
    emp_thetaLB[i,] <- apply(theta_burn[[i]], 2, function(x){quantile(x, 0.025)})
    emp_thetaMED[i,]<- apply(theta_burn[[i]], 2, function(x){quantile(x, 0.5)})
    emp_thetaUB[i,] <- apply(theta_burn[[i]], 2, function(x){quantile(x, 0.975)})
  }
  
  tru_cover[nr,] <-  colMeans((theta_i_sim)    < emp_thetaUB & (theta_i_sim)    > emp_thetaLB)
  hat_cover[nr,] <- colMeans(log(estvals) < emp_thetaUB & log(estvals) > emp_thetaLB)
  
  # Collect the "bad thetas": one of the true elements of theta_i 
  # isn't contained in the 95% credible interval
  bad_thetas <- which(apply(theta_i_sim < emp_thetaUB & theta_i_sim > emp_thetaLB, 1, all)==F)
  # cbind(emp_thetaLB[bad_thetas,], theta_i_sim[bad_thetas,], emp_thetaUB[bad_thetas,])
  
  theta_cols <- c("sigma2","phi","kappa") 
  
  par(mfrow=c(1,3))
  for (j in bad_thetas) {
    for(i in 1: P) {
      plot(theta[[j]][(burn+1):iters,i], 
           main = paste("storm",j,loc_int$loc[j],theta_cols[i]), type="l")#,
      # ylim=c(apply(theta_hat, 2, min)[i] + 0.2, apply(theta_hat, 2, max)[i] - 0.2))
      abline(h=apply(theta[[j]][(burn+1):iters,], 2, mean)[i], col="green", lty=2)
      abline(h=theta_i_sim[j,i], lwd=2, col="blue")
      # abline(h=theta_i_hat_sim[j,i], lwd=2, col="orange")
      abline(h=log(estvals)[j,i], lwd=2, col="green3")
    }
  }
  
  # Check if each element of true B and Sigma_theta is contained in the credible intervals
  emp_B    <- apply(B_burn, 2, function(x){quantile(x,0.5)})
  emp_B_UB <- apply(B_burn, 2, function(x){quantile(x,0.975)})
  emp_B_LB <- apply(B_burn, 2, function(x){quantile(x,0.025)})
  
  truthB <= emp_B_UB & truthB >= emp_B_LB
  
  emp_Sigma_theta    <- apply(Sigma_burn, 2, function(x){quantile(x,0.5)})
  emp_Sigma_theta_LB <- apply(Sigma_burn, 2, function(x){quantile(x,0.025)})
  emp_Sigma_theta_UB <- apply(Sigma_burn, 2, function(x){quantile(x,0.975)})
  
  truthSig <= emp_Sigma_theta_UB & truthSig >= emp_Sigma_theta_LB
}

# write.csv(hat_cover, file = paste0("~/NAM-Model-Validation/csv/simsMLEout/hat_cover_", 
#                                    if(trueBSig){"trueBSig"}else{"estBSig"}, "_nreps",nreps,
#                                    "_iters", iters,"_",nowtime(),".csv"))
# 
# write.csv(tru_cover, file = paste0("~/NAM-Model-Validation/csv/simsMLEout/tru_cover_", 
#                                    if(trueBSig){"trueBSig"}else{"estBSig"}, "_nreps",nreps,
#                                    "_iters", iters,"_",nowtime(),".csv"))

apply(hat_cover, 2, mean)
apply(tru_cover, 2, mean)

# hat_cover_estBSig_estPsi <- hat_cover
# tru_cover_estBSig_estPsi <- tru_cover

# hat_cover_true_estPsi <- hat_cover
# tru_cover_true_estPsi <- tru_cover

par(mfrow=c(1,3))
xseq <- seq(0.5,1,by=0.005)
for (i in 1:P) {
  hist(tru_cover[,i], main = i, freq=F)
  lines(xseq, dnorm(xseq, mean = .95, sd = sqrt(.95*.05/N)))
}

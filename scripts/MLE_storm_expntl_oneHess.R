## MLE calculation for exponential covariance parameters
## Hessian calculation included to approximate asymptotic covariance matrix
## Steve Walsh, October 2020

## This script will produce a pdf containing images of log lik surface 
## with respect to both parameterizations. It will also only compute
## the numerical approximation to the Hessian (rather than all ways)

library(pracma)
library(raster)

sim_vals <- F
sim_dom <- F
subtractPWmean <- F
if(sim_vals & subtractPWmean) stop("if sim_vals=F, then subtractPWmean should be F")

## lo_sim = length.out for simulated precip error grid
## square grid will have lo_sim^2 pixels on [0, 10]^2
if(sim_dom){lo_sim <- 32}  
if(sim_dom & !sim_vals) stop("if sim_vals=F, then you must have sim_dom=F")

trial <- F
hess_calc <- T
writefiles <- T

get_prec_mtx <- T

plot.it <- F
lo <- 21 #length.out for lkhd grid (do odd so MLE in middle)


if(sim_vals){ 
  seed <- 4
  set.seed(seed) 
}

if(writefiles & plot.it){
  pdf(paste0("pdf/MLE_storm_expntl_oneHess",
             if(sim_vals){"_sim_vals"},if(trial){"_trial"},"_lo",lo,".pdf"))
}

start.time <- Sys.time()

if(subtractPWmean){
  storm_csvs <- list.files("csv/error_df_sqrt/subtractPWmeanT_flat", full.names = T)
} else {
  storm_csvs <- list.files("csv/error_df_sqrt/subtractPWmeanF", full.names = T)
}
small_locs <-  c(11, 18, 17)#,  3, 15, 25,  6, 46, 34, 23, 30, 41, 21, 33, 38, 45, 43, 19, 35, 31) #fastest 20
 #c(3,11,15,17,18,25) #storms less than 3000 pixels
  

if(trial){storms_to_eval <- small_locs} else {storms_to_eval <- 1:length(storm_csvs)}
if(sim_dom){storms_to_eval <- 1:200}

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

## Keep track of MLE estimates for phi and sigma^2, 
## as well as counts of evals of ll fn in optimize, times, and number of pixels
## Hessian calculations, compare three different methods 
## (2 analytical, 1 numerical from pracma package)
phis <- sigs <- cnts <- tims <- pixs <- c()
myhessvecs <- myhessvecsEXP <- pkghessvecs <- 
  mythetahessvecs <- pkgthetahessvecs <- matrix(NA, length(storms_to_eval), 4)
all_grids <- theta1s <- theta2s <- list()

###################################
# Exponential covariance function #
###################################

eps <- sqrt(.Machine$double.eps)
covfunc.exponential <- function(t,phi,sigma2) {sigma2 * exp(-t/phi) + diag(eps, nrow=ncol(t))}

# ## Slow version of negative log likehood calculation
# nl <- function(phi, D, Y)
# {
#   n <- length(Y)
#   # print(n)
#   K <- exp(-D/phi)
#   # print(K[1:6,1:6])
#   Ki <- solve(K)
#   # print(Ki[1:6,1:6])
#   ldetK <- determinant(K, logarithm=TRUE)$modulus
#   # print(ldetK)
#   # ldetR <- 2*sum(log(diag(R)))
#   ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
#   # print(ll)
#   counter <<- counter + 1
#   return(-ll)
# }

## Negative profile log likelihood fn (1D minimization)
## Use Cholesky decomp to speed up computations
## cormatrix = R^T %*% R
## ldetR obtains the log determinant of the cov mtx
## yTCiy obtains the quadratic term
nl_new <- function(phi, D, Y)
{
  n <- length(Y)
  cormatrix <- covfunc.exponential(D, phi = phi, sigma2 = 1)
  R <- chol(cormatrix)
  yTCiy <- sum(backsolve(R, Y, transpose = T)^2)
  sigma2hat <- yTCiy/n
  ldetR <- 2*sum(log(diag(R)))
  ll <- -1/(2*sigma2hat)*yTCiy - (1/2)*ldetR -n/2*log(sigma2hat)
  counter <<- counter + 1
  return(list(obj=-ll, sigma2hat=sigma2hat, cormatrix=cormatrix))
}

## After obtaining MLEs for phi and sigma^2 with profile likelihood,
## Use this log likelihood function to evaluate Hessians w/ pracma pkg
ll_new <- function(vec, D=t, Y=x)
{
  sigma2hat <- vec[1]
  phi       <- vec[2]
  n <- length(Y)
  cormatrix <- covfunc.exponential(D, phi = phi, sigma2 = 1)
  R <- chol(cormatrix)
  yTCiy <- sum(backsolve(R, Y, transpose = T)^2)
  ldetR <- 2*sum(log(diag(R)))
  ll <- -1/(2*sigma2hat)*yTCiy - (1/2)*ldetR -n/2*log(sigma2hat)
  return(ll)
}

## log-likelihood evaluation wrt theta parameterization
ll_theta <- function(thetavec, D=t, Y=x)
{
  sigma2hat <- exp(thetavec[2])
  phi       <- sigma2hat/exp(thetavec[1])
  n <- length(Y)
  cormatrix <- covfunc.exponential(D, phi = phi, sigma2 = 1)
  R <- chol(cormatrix)
  yTCiy <- sum(backsolve(R, Y, transpose = T)^2)
  ldetR <- 2*sum(log(diag(R)))
  ll <- -1/(2*sigma2hat)*yTCiy - (1/2)*ldetR -n/2*log(sigma2hat)
  return(ll)
}

for (i in 1:length(storms_to_eval)) {
  
  print(i) # i <- 1
  if(sim_dom){
    xaxis <- seq(0,10,length.out = lo_sim)
    yaxis <- seq(0,10,length.out = lo_sim)
    s <- expand.grid(xaxis, yaxis)
    # storm <- read.csv(storm_csvs[1], row.names = 1) ## quick hack to do N=50, 500... 
    # s <- storm[,2:3]                                ## of a real location with sim_vals
  } else {
    storm <- read.csv(storm_csvs[storms_to_eval[i]], row.names = 1)
    s <- storm[,2:3]
  }
  
  pixs[i] <- nrow(s)
  
  # plot(s)
  t <- as.matrix(dist(s))
  n <- nrow(s)
  
  if(sim_vals){
    # Without nugget effect
    trueSigma2 <- 4
    truePhi <- 1.5
    # plot(seq(0,1,0.001),covfunc.exponential(seq(0,1,0.001),phi,sigma2),type="l")
    
    true_cov <- covfunc.exponential(t,truePhi,trueSigma2)
    # image(covmatrix)
    x <- t(chol(true_cov)) %*% rnorm(nrow(s))
  } else {
    x <- storm[,1]
  }
  
  # plot(rasterFromXYZ(cbind(s,x)))
  
  counter <- 0
  
  # system.time({ out_new <- nl_new(phi = 2, D = t, Y = x) })
  # print(counter)
  # counter <- 0
  # system.time({ out_old <-     nl(phi = 2, D = t, Y = x) })
  
  tic <- proc.time()[3]
  
  mle_new <- optimize(function(phi) nl_new(phi, D = t, Y = x)[[1]], interval = c(0,1000))
  phihat <- mle_new$minimum
  
  out <- nl_new(phihat, D = t, Y=x)
  sigma2hat <- out$sigma2hat
  cormatrix <- out$cormatrix
  
  if(hess_calc){
    # ## Hessian calculations (phi, sigma2hat, cormatrix all exported from nl_new)
    # yTCiy <- n*sigma2hat
    # d2lds2    <- n/(2*sigma2hat^2)-(yTCiy)/sigma2hat^3 
    # # d2lds2EXP <- -n/(2*sigma2hat^2)
    # Ci <- solve(cormatrix)
    # Cp <- cormatrix * (t/phihat^2)
    # CiCp <- Ci %*% Cp
    # Cpp<- (Cp * t/phihat^2) - (cormatrix * (2*t)/phihat^3)
    # CiCpp <- Ci %*% Cpp
    # d2ldsdp    <- -1/(2*sigma2hat^2)*t(x)%*%CiCp%*%Ci%*%x
    # # d2ldsdpEXP <- -1/(2*sigma2hat)*sum(diag(CiCp))
    # A <- -(CiCp%*%CiCp) + CiCpp
    # B <- CiCpp%*%Ci - 2*CiCp%*%CiCp%*%Ci
    # d2ldp2    <- -0.5*sum(diag(A)) + 1/(2*sigma2hat)*t(x)%*%B%*%x
    # # d2ldp2EXP <- -0.5*sum(diag(CiCp%*%CiCp))
    # 
    # myHess <- c(d2lds2, d2ldsdp, d2ldsdp, d2ldp2)
    # # myHessEXP <- c(d2lds2EXP, d2ldsdpEXP, d2ldsdpEXP, d2ldp2EXP)
    
    pkgHess <- pracma::hessian(f = ll_new, x0=c(sigma2hat, phihat))
    
    # ## multivariate delta method for theta_1, theta_2
    # ## theta_1 = log(sigma^2/phi)
    # ## theta_2 = log(sigma^2)
    # ## M contains the necessary derivatives wrt to these fns & sigma^2, phi
    M <- matrix(c(1/sigma2hat, 1/sigma2hat,-1/phihat, 0), ncol = 2)
    covmtx <- solve(matrix(-pkgHess, 2, 2))
    pkgthetaHess <- -solve(M %*% covmtx %*% t(M))

    # ## evaluate Hessian numerically w/ theta parameterization
    # pkgthetaHess <- pracma::hessian(f = ll_theta, x0=c(log(sigma2hat/phihat),log(sigma2hat)))

    # ## results are the same?
    # all.equal(round(solve(-pkgthetaHess),3),round(mythetaHess,3))
    
  }
  
  toc <- proc.time()[3]
  
  # print(counter)
  print(sigma2hat)
  print(mle_new$minimum)
  if(hess_calc){
    # print(sqrt(-diag(solve(matrix(myHess,ncol = 2)))))
    # print(sqrt(-diag(solve(matrix(myHessEXP,ncol = 2)))))
    print(sqrt(-diag(solve(matrix(pkgHess,ncol = 2)))))
    # print(sqrt(-diag(solve(matrix(mythetaHess,ncol = 2)))))
    print(sqrt(-diag(solve(matrix(pkgthetaHess,ncol = 2)))))
  }
  
  # system.time(mle_old <- optimize(nl, interval = c(0,10), D = t, Y = x))
  # print(counter)
  
  # mle_old$minimum
  phis[i] <- mle_new$minimum
  sigs[i] <- sigma2hat
  cnts[i] <- counter
  tims[i] <- toc-tic
  if(hess_calc){
    # myhessvecs[i,] <- myHess
    # myhessvecsEXP[i,] <- myHessEXP
    pkghessvecs[i,] <- pkgHess
    # mythetahessvecs[i,] <- mythetaHess
    pkgthetahessvecs[i,] <- pkgthetaHess
  }
  
  if(plot.it){
    ## Plot the log-likelihood surface wrt to orig parameterization
    orig_sds <- sqrt(diag(-solve(pkgHess)))
    
    ## convert thetas to sigma^2, phi for calculations
    sigma2vec <- seq(sigma2hat - 1.959964*orig_sds[1],
                     sigma2hat + 1.959964*orig_sds[1], length.out = lo)	#eg: 0.99, 1.11
    phivec <- seq(phihat - 1.959964*orig_sds[2],
                  phihat + 1.959964*orig_sds[2], length.out = lo) #eg: 0.82, 2.20
    
    ## remove any negatives
    bad <- unique(which(phivec < 0), which(sigma2vec < 0))
    if(length(bad)>0){
      sigma2vec <- sigma2vec[-bad]
      phivec <- phivec[-bad]
    }

    
    
    llgrid <- matrix(NA, nrow = length(phivec), ncol = length(sigma2vec))
    
    for (k in 1:length(phivec)) {
      # print(k)
      for(j in 1:length(sigma2vec)){
        llgrid[k,j] <- ll_new(vec = c(sigma2vec[k], phivec[j]))
      }
    }
    
    par(mfrow=c(1,2))
    image(sigma2vec, phivec, llgrid)
    contour(sigma2vec, phivec, llgrid, add=T)
    
    ## Plot the log-likelihood surface wrt to new parameterization
    ## theta_1 = log(sigma^2/phi)
    ## theta_2 = log(sigma^2)
    
    ## MLEs with theta parameterization
    theta1hat <- log(sigma2hat/phihat)
    theta2hat <- log(sigma2hat)
    
    ## corresponding standard deviations
    theta_sds <- sqrt(diag(-solve(pkgthetaHess)))
    
    ## vectors for plotting the 2D log-likelihood fn
    theta1vec <- seq(theta1hat - 1.959964*theta_sds[1],
                     theta1hat + 1.959964*theta_sds[1], length.out = lo)	#eg: 0.99, 1.11
    theta2vec <- seq(theta2hat - 1.959964*theta_sds[2],
                     theta2hat + 1.959964*theta_sds[2], length.out = lo) #eg: 0.82, 2.20
    
    thetagrid <- matrix(NA, nrow = length(theta1vec), ncol = length(theta2vec))
    
    for (k in 1:length(theta1vec)) {
      # print(k)
      for(j in 1:length(theta2vec)){
        thetagrid[k,j] <- ll_theta(thetavec = c(theta1vec[k], theta2vec[j]))
      }
    }
    
    all_grids[[i]] <- thetagrid
    theta1s[[i]] <- theta1vec
    theta2s[[i]] <- theta2vec
    
    image(theta1vec, theta2vec, thetagrid)
    contour(theta1vec, theta2vec, thetagrid, add=T)
    
    image(sigma2vec, phivec, exp(llgrid-min(llgrid)))
    contour(sigma2vec, phivec, exp(llgrid-min(llgrid)), add=T)
    image(theta1vec, theta2vec, exp(thetagrid-min(thetagrid)))
    contour(theta1vec, theta2vec, exp(thetagrid-min(thetagrid)), add=T)
    
    par(mfrow=c(1,2))
    ind <- which(theta1hat == theta1vec)
    if(!(all(exp(thetagrid[,ind]) == Inf) | all(exp(thetagrid[ind,]) == Inf))){
    plot(theta1vec, exp(thetagrid[,ind]-min(thetagrid[,ind])), type = "l", main="slice of likhd at theta_2_hat")
    plot(theta2vec, exp(thetagrid[ind,]-min(thetagrid[ind,])), type = "l", main="slice of likhd at theta_1_hat")
    
    plot(theta1vec, thetagrid[,ind], type = "l", main="slice of loglik at theta_2_hat")
    plot(theta2vec, thetagrid[ind,], type = "l", main="slice of loglik at theta_1_hat")
    }
  }
  
}

# nowtime <- function(){gsub(gsub(gsub(Sys.time(),pattern = " ", replacement = ""),
#                                 pattern="-", replacement=""),pattern=":", replacement="")}

colMeans(cbind(sigs, phis, cnts, tims), na.rm = T)

if(plot.it){
  par(mfrow=c(2,2))
  hist(phis, main = expression(phi), xlab = expression(phi))
  hist(sigs, main = expression(sigma^2), xlab = expression(sigma^2))
  hist(log(phis), main = "log phi")#, xlab = expression(phi))
  hist(log(sigs), main = "log sigma2")#, xlab = expression(sigma^2))
  
  par(mfrow=c(1,1))
  plot(phis,sigs, main = "Correlation of Phi and Sigma^2")
  if(length(storms_to_eval) > 3){
    abline(lm(sigs~phis))
    summary(lm(sigs~phis))
  }
  
  par(mfrow=c(1,3))
  hist(phis, main = expression(phi), xlab = expression(phi))
  if(sim_vals){ abline(v=truePhi, col="blue") }
  hist(sigs, main = expression(sigma^2), xlab = expression(sigma^2))
  if(sim_vals){ abline(v=trueSigma2, col="blue") }
  hist(sigs/phis, main = expression(sigma^2/phi), xlab = expression(sigma^2/phi))
  if(sim_vals){ abline(v=trueSigma2/truePhi, col="blue") }
  
}

sig_phi_sds <- thetas_sds <- matrix(NA, length(storms_to_eval), 2)
if(length(storms_to_eval)==1){
    sig_phi_sds <- sqrt(diag(-solve(matrix(pkghessvecs, 2, 2))))
    thetas_sds  <- sqrt(diag(-solve(matrix(pkgthetahessvecs, 2, 2))))
} else {
  for (k in 1:length(storms_to_eval)) {
    sig_phi_sds[k,] <- sqrt(diag(-solve(matrix(pkghessvecs[k,], 2, 2))))
    thetas_sds[k,]  <- sqrt(diag(-solve(matrix(pkgthetahessvecs[k,], 2, 2))))
  }
}


if(plot.it & length(storms_to_eval) > 1){
  ## hists of theta1's and theta2's
  par(mfrow=c(1,2))
  hist(log(sigs/phis), main = "theta1 hats")
  hist(log(sigs), main = "theta2 hats")
  
  ## hists of asymptotic sd's for both params
  hist(sig_phi_sds[,1], main = "asymp SD for sigma^2")
  hist(sig_phi_sds[,2], main = "asymp SD for phi")
  
  hist(thetas_sds[,1], main = "asymp SD for log(sigma^2/phi)")
  hist(thetas_sds[,2], main = "asymp SD for log(sigma^2)")
}

## Write files
if(!dir.exists("csv/myMLEresults/myMLEs/")){
  dir.create("csv/myMLEresults/myMLEs/", recursive = T)
}
if(!dir.exists("csv/myMLEresults/pkgthetahessvecs/")){
  dir.create("csv/myMLEresults/pkgthetahessvecs/", recursive = T)
}
if(!dir.exists("csv/myMLEresults/pkghessvecs/")){
  dir.create("csv/myMLEresults/pkghessvecs/", recursive = T)
}
if(!dir.exists("RData/myMLE_precs/")){
  dir.create("RData/myMLE_precs/", recursive = T)
}

if(writefiles){
  write.csv(cbind(phis, sigs, cnts, tims), 
            file=paste0("csv/myMLEresults/myMLEs/",
                        if(storms_to_eval[1] < 10){"0"}, storms_to_eval[1],
                        if(subtractPWmean){"subtractPWmean"},
                        if(sim_vals){paste0("seed",seed,"_")},
                        if(trial){"trial_"},if(sim_vals){"sim_vals"},
                        if(sim_dom){paste0("_sim_dom", lo_sim)},".csv"))
  if(hess_calc){
    # write.csv(myhessvecs, paste0("csv/myMLEresults/myhessvecs_",
    #                              if(storms_to_eval[1] < 10){"0"}, storms_to_eval[1],
    #                              if(sim_vals){paste0("seed",seed,"_")},
    #                              if(trial){"trial_"},if(sim_vals){"sim_vals"},".csv"))
    # write.csv(myhessvecsEXP, paste0("csv/myMLEresults/myhessvecsEXP_",
    #                                 if(storms_to_eval[1] < 10){"0"}, storms_to_eval[1],
    #                                 if(sim_vals){paste0("seed",seed,"_")},
    #                                 if(trial){"trial_"},if(sim_vals){"sim_vals"},".csv"))
    write.csv(pkghessvecs, paste0("csv/myMLEresults/pkghessvecs/",
                                  if(storms_to_eval[1] < 10){"0"}, storms_to_eval[1],
                                  if(subtractPWmean){"subtractPWmean"},
                                  if(sim_vals){paste0("seed",seed,"_")},
                                  if(trial){"trial_"},if(sim_vals){"sim_vals"},
                                  if(sim_dom){paste0("_sim_dom", lo_sim)},".csv"))
    # write.csv(mythetahessvecs, paste0("csv/myMLEresults/mythetahessvecs_",
    #                                   if(storms_to_eval[1] < 10){"0"}, storms_to_eval[1],
    #                                   if(sim_vals){paste0("seed",seed,"_")},
    #                                   if(trial){"trial_"},if(sim_vals){"sim_vals"},".csv"))
    write.csv(pkgthetahessvecs, paste0("csv/myMLEresults/pkgthetahessvecs/",
                                       if(storms_to_eval[1] < 10){"0"}, storms_to_eval[1],
                                       if(subtractPWmean){"subtractPWmean"},
                                       if(sim_vals){paste0("seed",seed,"_")},
                                       if(trial){"trial_"},if(sim_vals){"sim_vals"},
                                       if(sim_dom){paste0("_sim_dom", lo_sim)},".csv"))
  }
}

if(get_prec_mtx){
  prec_mtx <- solve(sigma2hat * cormatrix)
  save(prec_mtx, file=paste0("RData/myMLE_precs/", 
                             if(storms_to_eval[1] < 10){"0"},
                             storms_to_eval[1],
                             if(subtractPWmean){"subtractPWmean"},
                             ".RData"))
}

end.time <- Sys.time()
end.time - start.time

if(writefiles & plot.it){ dev.off() }

if(sim_vals){
  trueTheta1 <- log(trueSigma2/truePhi)
  trueTheta2 <- log(trueSigma2)
  
  phix <- phis[complete.cases(phis)]
  sigx <- sigs[complete.cases(sigs)]
  thet1x <- log(sigx/phix)
  thet2x <- log(sigx)
  
  mean(thet1x + 1.96*thetas_sds[,1] > trueTheta1 & thet1x - 1.96*thetas_sds[,1] < trueTheta1)
  mean(thet2x + 1.96*thetas_sds[,2] > trueTheta2 & thet2x - 1.96*thetas_sds[,2] < trueTheta2)
}

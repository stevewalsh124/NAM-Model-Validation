# By M.A.R. Ferreira, 2014. Updated by M.A.R. Ferreira, 2016.
library(raster)
library(fields)
library(LaplacesDemon)

load("~/NAM-Model-Validation/RData/all_storm_res.RData") #load all_storm_res
load("~/NAM-Model-Validation/RData/emp_B")
load("~/NAM-Model-Validation/RData/emp_Sigma_theta")
loc_int <- read.csv("~/NAM-Model-Validation/csv/storm_levels_and_locs.csv", row.names = 1)

n <- 1
if(n != 1){
  d <- loc_int
  loc_int <-  do.call("rbind", replicate(n, d, simplify = FALSE))
}



as.square <- function(mat){matrix(mat, nrow=sqrt(length(mat)),ncol=sqrt(length(mat)))}
emp_B <- as.square(emp_B)
emp_Sigma_theta <- as.square(emp_Sigma_theta)

grp_B <- cbind(exp(emp_B[,1]), exp(emp_B[,1] + emp_B[,2]), exp(emp_B[,1]+emp_B[,3]))
N <- 47 * n; P <- 3

sim_B <- as.square(c(0,0,0,3,3,3,9,9,9))
xi <- with(loc_int, model.matrix(~ storm_locs))

set.seed(1)
theta_i_sim <- matrix(nrow=N, ncol=P)
# for (i in 1:47) {theta_i_sim[i,] <- rmvn(n = 1, mu = t(emp_B %*% xi[i, ]), Sigma = emp_Sigma_theta)}

for (i in 1:length(levels(loc_int$storm_locs))) {
  if(i==1){
    theta_i_sim[which(loc_int$storm_locs == levels(loc_int$storm_locs)[i]),] <- 
      rmvn(n = length(which(loc_int$storm_locs == levels(loc_int$storm_locs)[i])), 
           mu = exp(emp_B[,i]), 
           Sigma = emp_Sigma_theta)
  } else {
    theta_i_sim[which(loc_int$storm_locs == levels(loc_int$storm_locs)[i]),] <- 
      rmvn(n = length(which(loc_int$storm_locs == levels(loc_int$storm_locs)[i])), 
           mu = exp(emp_B[,1]+ emp_B[,i]), 
           Sigma = emp_Sigma_theta)
  }
  
  print(apply(theta_i_sim[which(loc_int$storm_locs == levels(loc_int$storm_locs)[i]),], 2, mean))
}

t(grp_B)
save(theta_i_sim, file="~/NAM-Model-Validation/RData/theta_i_sim")


# s <- 11
# 
# args <- commandArgs(TRUE)
# if(length(args) > 0)
#   for(i in 1:length(args))
#     eval(parse(text=args[[i]]))

par(mfrow=c(3,2))
for (st in 1:47) {
  print(paste("this is storm", st))
  storm <- read.csv(list.files("~/NAM-Model-Validation/csv/error_df/subtractPWmeanF/", full.names = T)[st], row.names = 1)
  
  # Simulated locations from unit square
  # xaxis <- seq(0,1,0.01)#storm$x#seq(0,1,0.01)
  # yaxis <- seq(0,1,0.01)#storm$y#seq(0,1,0.01)
  # s <- matrix(NA,nrow=length(xaxis)*length(yaxis),ncol=2)
  # for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) s[i+(j-1)*length(xaxis),] <- c(xaxis[i],yaxis[j])
  
  # Simulated locations based on actual locations
  s <- storm[,c(2,3)]
  # plot(s)
  
  t <- as.matrix(dist(s))
  # image(t)
  
  ##############################
  # Matern covariance function #
  ##############################
  
  library(geoR)
  
  # With nugget effect
  tau2 <- 0
  sigma2 <- theta_i_sim[st,1]
  phi <- theta_i_sim[st,2]
  nu <- theta_i_sim[st,3]
  # plot(seq(0,1,0.001),c(tau2,rep(0,length(seq(0,1,0.001))-1))+sigma2*matern(seq(0,1,0.001),phi=(1/phi),kappa=nu))
  covmatrix <- diag(tau2,nrow(t)) + sigma2*matern(t,phi=(1/phi),kappa=nu)
  # image(covmatrix)
  x <- t(chol(covmatrix)) %*% rnorm(nrow(s))
  # z <- matrix(NA,nrow=length(xaxis),ncol=length(yaxis))
  # for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) z[i,j] <- x[i+(j-1)*length(xaxis)]
  # persp(x=storm$x, y=storm$y, z=z)
  # image(z)
  write.csv(cbind(storm[,c(2,3)], z=x), file = paste0("~/NAM-Model-Validation/csv/sim_df/sim",st,".csv"), row.names = F)
  
  plot(rasterFromXYZ(cbind(storm[,c(2,3)], x)), main = st); US(add=T)

  plot(rasterFromXYZ(storm[,c(2,3,1)]), main = loc_int$storm_locs[st]); US(add=T)
}








# 
# likfit_time <- system.time({
#   MLtemp <- likfit(as.geodata(cbind(s,x)), ini.cov.pars = c(1,1), lik.method = "ML",
#                    cov.model = "matern", fix.kappa = FALSE, fix.nugget = FALSE, kappa=.501, hessian=T,
#                    method = "L-BFGS-B", lower = c(0,0,0), upper = c(10,10,10))
# })
# 
# loglikhd.spatial <- function (geodata, theta_GRF)
# {
#   cov.pars = theta_GRF[1:2]
#   kappa = theta_GRF[3]
#   obj.model = NULL
#   cov.model = "matern"
#   coords = geodata$coords
#   data = geodata$data
#   lambda = 1
#   psiR = 1
#   psiA = 0
#   trend = "cte"
#   nugget = MLtemp$nugget
#   method.lik = "ML"
#   compute.dists = TRUE
#   realisations = NULL
# 
#   sigmasq <- cov.pars[1]
#   phi <- cov.pars[2]
# 
#   if (is.null(realisations))
#     realisations <- as.factor(rep(1, length(data)))
#   else realisations <- as.factor(realisations)
#   nrep <- length(levels(realisations))
#   if (kappa < 1e-04)
#     return(-(.Machine$double.xmax^0.5))
#   if ((nugget + sigmasq) < 1e-16)
#     return(-(.Machine$double.xmax^0.5))
#   if (missing(geodata))
#     xmat <- unclass(trend.spatial(trend = trend, geodata = list(coords = coords,
#                                                                 data = data)))
#   else xmat <- unclass(trend.spatial(trend = trend, geodata = geodata))
#   if (nrow(xmat) != nrow(coords))
#     stop("coords and trend have incompatible sizes")
#   beta.size <- ncol(xmat)
#   xmat <- split(as.data.frame(unclass(xmat)), realisations)
#   xmat <- lapply(xmat, as.matrix)
#   vecdist <- function(x) {
#     as.vector(dist(x))
#   }
#   if (psiR != 1 | psiA != 0) {
#     coords.c <- coords.aniso(coords, aniso.pars = c(psiA,
#                                                     psiR))
#     .likGRF.dists.vec <- lapply(split(as.data.frame(coords.c),
#                                       as.factor(realisations)), vecdist)
#   }
#   else if (compute.dists)
#     .likGRF.dists.vec <- lapply(split(as.data.frame(coords),
#                                       as.factor(realisations)), vecdist)
#   z <- data
#   if (abs(lambda - 1) < 1e-04)
#     log.jacobian <- 0
#   else {
#     if (any(z <= 0))
#       stop("Transformation not allowed for zero or negative data")
#     data <- z^(lambda - 1)
#     if (any(data <= 0))
#       log.jacobian <- log(prod(data))
#     else log.jacobian <- sum(log(data))
#     data <- NULL
#     if (abs(lambda) < 1e-04)
#       data <- log(z)
#     else data <- ((z^lambda) - 1)/lambda
#   }
#   data <- split(data, as.factor(realisations))
#   sumnegloglik <- 0
#   for (i in 1:nrep) {
#     n <- length(data[[1]])
#     if ((phi < 1e-16) | (sigmasq < 1e-16)) {
#       V <- list(varcov = diag(x = (nugget + sigmasq), n),
#                 log.det.to.half = (n/2) * log(nugget + sigmasq))
#     }
#     else {
#       V <- varcov.spatial(dists.lowertri = .likGRF.dists.vec[[i]],
#                           cov.model = cov.model, kappa = kappa, nugget = nugget,
#                           cov.pars = c(sigmasq, phi), det = TRUE)
#     }
#     if (!is.null(V$crash.parms)) {
#       cat("varcov.spatial: improper matrix for following the given parameters:")
#       print(V$crash.parms)
#       stop()
#     }
#     ivx <- solve(V$varcov, xmat[[i]])
#     xivx <- crossprod(ivx, xmat[[i]])
#     betahat <- .solve.geoR(xivx, crossprod(ivx, data[[i]]))
#     res <- data[[i]] - xmat[[i]] %*% betahat
#     ssres <- drop(crossprod(res, solve(V$varcov, res)))
#     if (method.lik == "ML") {
#       negloglik <- (n/2) * (log(2 * pi)) + V$log.det.to.half +
#         0.5 * ssres
#     }
#     if (method.lik == "RML") {
#       choldet <- sum(log(diag(chol(xivx))))
#       negloglik <- V$log.det.to.half + 0.5 * ssres + choldet
#       xx.eigen <- eigen(crossprod(xmat[[i]]), symmetric = TRUE,
#                         only.values = TRUE)
#       negloglik <- negloglik + ((n - beta.size)/2) * (log(2 *
#                                                             pi)) - 0.5 * sum(log(xx.eigen$values))
#     }
#     sumnegloglik <- sumnegloglik + negloglik
#   }
#   sumnegloglik <- sumnegloglik - log.jacobian
#   if (sumnegloglik > (.Machine$double.xmax^0.5))
#     sumnegloglik <- .Machine$double.xmax^0.5
#   return(as.vector(sumnegloglik)) ### I removed a minus sign so that I minimize the -LL, aka maximize the LL
# }
# 
# 
# hur_opt <- list()
# hess_time <- system.time({hur_opt <- optim(par = c(MLtemp$cov.pars, MLtemp$kappa),
#                                            fn = loglikhd.spatial, geodata=as.geodata(cbind(s,x)),
#                                            hessian = T, method = "L-BFGS-B",
#                                            lower = c(0,0,0), upper = c(10,10,10))})
# 
# print(hur_opt$par);print(hur_opt$hessian);print(paste("The above was from storm #", co))
# 

# By M.A.R. Ferreira, 2014. Updated by M.A.R. Ferreira, 2016.
library(raster)
load("~/NAM-Model-Validation/RData/all_storm_res.RData")
all_storm_res
loc_int <- read.csv("~/NAM-Model-Validation/csv/storm_levels_and_locs.csv", row.names = 1)
storm <- read.csv(list.files("~/NAM-Model-Validation/csv/errordf", full.names = T)[3], row.names = 1)
plot(rasterFromXYZ(storm[,c(2,3,1)]))
xaxis <- seq(0,1,0.01)#storm$x#seq(0,1,0.01)
yaxis <- seq(0,1,0.01)#storm$y#seq(0,1,0.01)

s <- matrix(NA,nrow=length(xaxis)*length(yaxis),ncol=2)
for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) s[i+(j-1)*length(xaxis),] <- c(xaxis[i],yaxis[j])
s <- storm[,c(2,3)]
# plot(s)

t <- as.matrix(dist(s))
# image(t)

##############################
# Matern covariance function #
##############################

library(geoR)

# With nugget effect

tau2 <- 0.05
sigma2 <- 1
phi <- 1.4
nu <- 0.5
# plot(seq(0,1,0.001),c(tau2,rep(0,length(seq(0,1,0.001))-1))+sigma2*matern(seq(0,1,0.001),phi=(1/phi),kappa=nu))

covmatrix <- diag(tau2,nrow(t)) + sigma2*matern(t,phi=(1/phi),kappa=nu)
# image(covmatrix)
x <- t(chol(covmatrix)) %*% rnorm(nrow(s))
z <- matrix(NA,nrow=length(xaxis),ncol=length(yaxis))
for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) z[i,j] <- x[i+(j-1)*length(xaxis)]
persp(x=storm$x, y=storm$y, z=z)
image(z)

plot(rasterFromXYZ(cbind(storm[,c(2,3)], x)))





likfit_time <- system.time({
  MLtemp <- likfit(as.geodata(cbind(s,x)), ini.cov.pars = c(1,1), lik.method = "ML",
                   cov.model = "matern", fix.kappa = FALSE, fix.nugget = FALSE, kappa=.501, hessian=T,
                   method = "L-BFGS-B", lower = c(0,0,0), upper = c(10,10,10))
})

loglikhd.spatial <- function (geodata, theta_GRF) 
{
  cov.pars = theta_GRF[1:2]
  kappa = theta_GRF[3]
  obj.model = NULL
  cov.model = "matern"
  coords = geodata$coords
  data = geodata$data
  lambda = 1
  psiR = 1
  psiA = 0
  trend = "cte"
  nugget = MLtemp$nugget
  method.lik = "ML"
  compute.dists = TRUE
  realisations = NULL
  
  sigmasq <- cov.pars[1]
  phi <- cov.pars[2]
  
  if (is.null(realisations)) 
    realisations <- as.factor(rep(1, length(data)))
  else realisations <- as.factor(realisations)
  nrep <- length(levels(realisations))
  if (kappa < 1e-04) 
    return(-(.Machine$double.xmax^0.5))
  if ((nugget + sigmasq) < 1e-16) 
    return(-(.Machine$double.xmax^0.5))
  if (missing(geodata)) 
    xmat <- unclass(trend.spatial(trend = trend, geodata = list(coords = coords, 
                                                                data = data)))
  else xmat <- unclass(trend.spatial(trend = trend, geodata = geodata))
  if (nrow(xmat) != nrow(coords)) 
    stop("coords and trend have incompatible sizes")
  beta.size <- ncol(xmat)
  xmat <- split(as.data.frame(unclass(xmat)), realisations)
  xmat <- lapply(xmat, as.matrix)
  vecdist <- function(x) {
    as.vector(dist(x))
  }
  if (psiR != 1 | psiA != 0) {
    coords.c <- coords.aniso(coords, aniso.pars = c(psiA, 
                                                    psiR))
    .likGRF.dists.vec <- lapply(split(as.data.frame(coords.c), 
                                      as.factor(realisations)), vecdist)
  }
  else if (compute.dists) 
    .likGRF.dists.vec <- lapply(split(as.data.frame(coords), 
                                      as.factor(realisations)), vecdist)
  z <- data
  if (abs(lambda - 1) < 1e-04) 
    log.jacobian <- 0
  else {
    if (any(z <= 0)) 
      stop("Transformation not allowed for zero or negative data")
    data <- z^(lambda - 1)
    if (any(data <= 0)) 
      log.jacobian <- log(prod(data))
    else log.jacobian <- sum(log(data))
    data <- NULL
    if (abs(lambda) < 1e-04) 
      data <- log(z)
    else data <- ((z^lambda) - 1)/lambda
  }
  data <- split(data, as.factor(realisations))
  sumnegloglik <- 0
  for (i in 1:nrep) {
    n <- length(data[[1]])
    if ((phi < 1e-16) | (sigmasq < 1e-16)) {
      V <- list(varcov = diag(x = (nugget + sigmasq), n), 
                log.det.to.half = (n/2) * log(nugget + sigmasq))
    }
    else {
      V <- varcov.spatial(dists.lowertri = .likGRF.dists.vec[[i]], 
                          cov.model = cov.model, kappa = kappa, nugget = nugget, 
                          cov.pars = c(sigmasq, phi), det = TRUE)
    }
    if (!is.null(V$crash.parms)) {
      cat("varcov.spatial: improper matrix for following the given parameters:")
      print(V$crash.parms)
      stop()
    }
    ivx <- solve(V$varcov, xmat[[i]])
    xivx <- crossprod(ivx, xmat[[i]])
    betahat <- .solve.geoR(xivx, crossprod(ivx, data[[i]]))
    res <- data[[i]] - xmat[[i]] %*% betahat
    ssres <- drop(crossprod(res, solve(V$varcov, res)))
    if (method.lik == "ML") {
      negloglik <- (n/2) * (log(2 * pi)) + V$log.det.to.half + 
        0.5 * ssres
    }
    if (method.lik == "RML") {
      choldet <- sum(log(diag(chol(xivx))))
      negloglik <- V$log.det.to.half + 0.5 * ssres + choldet
      xx.eigen <- eigen(crossprod(xmat[[i]]), symmetric = TRUE, 
                        only.values = TRUE)
      negloglik <- negloglik + ((n - beta.size)/2) * (log(2 * 
                                                            pi)) - 0.5 * sum(log(xx.eigen$values))
    }
    sumnegloglik <- sumnegloglik + negloglik
  }
  sumnegloglik <- sumnegloglik - log.jacobian
  if (sumnegloglik > (.Machine$double.xmax^0.5)) 
    sumnegloglik <- .Machine$double.xmax^0.5
  return(as.vector(sumnegloglik)) ### I removed a minus sign so that I minimize the -LL, aka maximize the LL
}


hur_opt <- list()
hess_time <- system.time({hur_opt <- optim(par = c(MLtemp$cov.pars, MLtemp$kappa), 
                                           fn = loglikhd.spatial, geodata=as.geodata(cbind(s,x)),
                                           hessian = T, method = "L-BFGS-B", 
                                           lower = c(0,0,0), upper = c(10,10,10))})

print(hur_opt$par);print(hur_opt$hessian);print(paste("The above was from storm #", co))


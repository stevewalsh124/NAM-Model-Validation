# Obtain MLE estimates and corresponding Hessian matrices
# for the spatial covariance parameters phi, kappa, sigma2
# Based on the csvs from LogPrecipVariograms.R

require(geoR)
# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
require(gstat)
require(fields)
library(stringr)
library(foreach)
# library(snow)
# library(doParallel)
# library(doSNOW)

#compare results with WLS approach, here alpha = 2, so nu = 1 in these WLS results:
# WLS_results_fixnu1 <- read.csv("csv/INLAvsWLS/svg.param.ests.error_KM_LCC_MaternSVGestimates.csv")

before <- Sys.time()

nowtime <- function(){gsub(gsub(gsub(Sys.time(),pattern = " ", replacement = ""),
                                pattern="-", replacement=""),pattern=":", replacement="")}

# pdf("~/NAM-Model-Validation/pdf/posterior_plots_allstorms_1to24_NAmask.pdf")

include.elevation = F #include/exclude elevation 
small_sample <- F
# I for INLA results
sig2om_Ivec <- c()
sig2ep_Ivec <- c()
range_Ivec <- c()
sig2om_Ivec_MED <- c()
sig2ep_Ivec_MED <- c()
range_Ivec_MED <- c()
name_Ivec <- c()

storms_to_eval <- 11
seed <- 1
fix.nug = T
nug.val <- 0#.01

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

print(paste("seed is", seed))
nug.prt <- gsub("\\.", "", nug.val)

storms_char <- paste0(if(min(storms_to_eval) < 10){"0"},min(storms_to_eval),
                      if(max(storms_to_eval) < 10){"0"},max(storms_to_eval))
cores <- 1
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

cl <- parallel::makeCluster(cores, outfile="")
doParallel::registerDoParallel(cl)

par_calcy <- foreach(co=1:cores, .combine= 'comb',
                     .init=list(list(), list(), list(), list(),
                                list(), list(), list(), list(),
                                list(), list(), list(), list(),
                                list(), list(), list())) %dopar% {
                                  require(geoR)
                                  library(INLA)
                                  require(gstat)
                                  require(fields)
                                  library(stringr)
                                  # i <- 17;co <- 1
                                  if(!exists("co")) co <- 1
                                  
                                  name_file <- list.files("~/NAM-Model-Validation/csv/error_df/subtractPWmeanT_flat",
                                                          full.names = T)[storms_to_eval[co]] #"csv/INLAvsWLS/error"
                                  name_start_INLA <- str_locate_all(pattern ='_', name_file)[[1]][,1]
                                  name_Ivec <- substr(name_file,
                                                      name_start_INLA[length(name_start_INLA)-1]+1,
                                                      name_start_INLA[length(name_start_INLA)]-1) #name_Ivec[i]
                                  print(name_Ivec)#name_Ivec[i]
                                  
                                  storm_file <- list.files(paste0("~/NAM-Model-Validation/csv/sim_df/seed",seed),
                                                           full.names = T)[storms_to_eval[co]] #"csv/INLAvsWLS/error"
                                  hurricane <- read.csv(storm_file)
                                  if(small_sample){hurricane <- hurricane[sample(1:dim(hurricane)[1], 1000),]}
                                  est.coord <- hurricane[,1:2]
                                  center.coord <- est.coord
                                  center.coord$x <- center.coord$x - min(center.coord$x)
                                  center.coord$y <- center.coord$y - min(center.coord$y)
                                  est.coord <- center.coord
                                  scale.coord <- center.coord
                                  scale_factor <- max(abs(center.coord$x),abs(center.coord$y))
                                  scale.coord$x <- scale.coord$x/scale_factor
                                  scale.coord$y <- scale.coord$y/scale_factor
                                  est.coord <- scale.coord
                                  est.data  <- hurricane[,3]
                                  
                                  ini.cov.likfit <- cbind(rep(c(0.5,1,1.5,2), each = 4), rep(c(0.05,0.5,1,2), 4)) #sigma2, phi
                                  
                                  likfit_res1 <- matrix(NA, nrow=1, ncol=5)
                                  likfitreml_res1 <- matrix(NA, nrow=1, ncol=5)
                                  colnames(likfitreml_res1) <- colnames(likfit_res1) <- c("sigma2","phi","tau2","kappa","beta")
                                  
                                  # geoR::likfit
                                  # MLtemp <- list()
                                  # likfit_time <- system.time({
                                  #   MLtemp <- likfit(as.geodata(hurricane), ini.cov.pars = ini.cov.likfit, lik.method = "ML",
                                  #                     cov.model = "matern", fix.kappa = FALSE, fix.nugget = FALSE, kappa=0.501)
                                  # })
                                  # 
                                  # likfit_res1 <- c(MLtemp$cov.pars, MLtemp$nugget, MLtemp$kappa, MLtemp$beta)
                                  # print(likfit_res1)
                                  
                                  MLtemp <- list()
                                  likfit_time <- system.time({
                                    MLtemp <- likfit(as.geodata(hurricane), ini.cov.pars = ini.cov.likfit, lik.method = "REML",
                                                     cov.model = "matern", fix.kappa = FALSE, fix.nugget = fix.nug, kappa=0.501,
                                                     nugget = nug.val,  hessian = T)#, method = "L-BFGS-B", 
                                    #lower = c(0,0,0), upper = c(10,10,10))
                                  })
                                  
                                  likfit_res1 <- c(MLtemp$cov.pars, MLtemp$nugget, MLtemp$kappa, MLtemp$beta)
                                  write.csv(likfit_res1, file = paste0("~/NAM-Model-Validation/csv/simsMLEout/seed",seed,"/simlik",
                                                                       if(storms_to_eval<10){"0"}, storms_to_eval, ".csv"))
                                  
                                  print(likfit_res1)
                                  print(MLtemp$info.minimisation.function$hessian)
                                  print(likfit_time)
                                  # REMLtemp <- likfit(as.geodata(cbind(s,as.vector(z))), ini.cov.pars = ini.cov.likfit[i,], lik.method = "REML", 
                                  # cov.model = "matern", fix.kappa = T, kappa = 1)
                                  # likfitreml_res1[i,] <- c(REMLtemp$cov.pars, REMLtemp$nugget, REMLtemp$kappa)
                                  # print(likfitreml_res1[i,])
                                  
                                  LL_REML <- function (geodata, theta_GRF) # based on geoR::loglik.GRF
                                  {
                                    coords = geodata$coords
                                    data = geodata$data 
                                    obj.model = NULL
                                    cov.model = "matern" 
                                    cov.pars = theta_GRF[1:2]
                                    nugget = 0 
                                    kappa = theta_GRF[3]
                                    lambda = 1 
                                    psiR = 1 
                                    psiA = 0
                                    trend = "cte"
                                    method.lik = "REML" 
                                    compute.dists = TRUE
                                    realisations = NULL
                                    
                                    if (!is.null(obj.model)) {
                                      if (!is.null(obj.model$cov.model)) 
                                        cov.model <- obj.model$cov.model
                                      if (!is.null(obj.model$cov.pars)) 
                                        cov.pars <- obj.model$cov.pars
                                      if (!is.null(obj.model$nugget)) 
                                        nugget <- obj.model$nugget
                                      if (!is.null(obj.model$kappa)) 
                                        kappa <- obj.model$kappa
                                      if (!is.null(obj.model$lambda)) 
                                        lambda <- obj.model$lambda
                                      if (!is.null(obj.model$psiR)) 
                                        psiR <- obj.model$psiR
                                      if (!is.null(obj.model$psiA)) 
                                        psiA <- obj.model$psiA
                                      if (!is.null(obj.model$trend)) 
                                        trend <- eval(obj.model$trend)
                                    }
                                    sigmasq <- cov.pars[1]
                                    phi <- cov.pars[2]
                                    if (method.lik == "REML" | method.lik == "reml" | method.lik == 
                                        "rml") 
                                      method.lik <- "RML"
                                    if (method.lik == "ML" | method.lik == "ml") 
                                      method.lik <- "ML"
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
                                    return(as.vector(sumnegloglik)) ## Removed a minus sign to min the -LL
                                  }
                                  

                                  hur_opt <- NULL
                                  hess_time <- system.time({hur_opt <- optim(par = c(MLtemp$cov.pars, MLtemp$kappa), 
                                                                             fn = LL_REML, geodata=as.geodata(hurricane),
                                                                             hessian = T, method = "L-BFGS-B", 
                                                                             lower = c(0,0,0), upper = c(10,10,10))})
                                  
                                  print(hur_opt$par);print(hur_opt$hessian);print(paste("The above was from storm #", storms_to_eval[co]))
                                  print(hess_time)
                                  this_hess <- hur_opt$hessian
                                  save(this_hess, file = paste0("~/NAM-Model-Validation/RData/RDatasim0/seed",seed,"/hess",
                                                                if(storms_to_eval<10){"0"}, storms_to_eval))
                                  
                                  ###########################################################################
                                  #SPDE/INLA 
                                  ###########################################################################
                                  #-- Create the mesh --#
                                  loc.domain = cbind(c(-1,-1,1,1),c(-1,1,-1,1))
                                  mesh =inla.mesh.2d(loc.domain = loc.domain - 0.5, max.edge = c(0.1,0.3), cutoff = .04, offset = .2)  
                                  mesh$n
                                  mesh$loc[,1:2] = mesh$loc[,1:2]+0.5   #to avoid infinite search for best mesh
                                  
                                  
                                  #-- Plot the triangulation --#
                                  # par(mar=c(0,2,2,0))
                                  # plot(mesh)
                                  # lines(sic.borders,lwd=3)
                                  # points(est.coord,pch=19, col="red", cex=1.2)
                                  
                                  #-- Construct the SPDE object --#
                                  spde = inla.spde2.matern(mesh=mesh, alpha = 1.5, theta.prior.prec = c(1e-4, 1e-4)) #alpha=2 by default
                                  # ?inla.spde2.matern
                                  
                                  #-- Create A matrices (mapping between the mesh nodes and station locations) --# 
                                  # observation matrix A
                                  A.est =
                                    inla.spde.make.A(mesh, loc=as.matrix(est.coord))
                                  # A.val =
                                  # inla.spde.make.A(mesh, loc=val.coord)
                                  # A.pred =
                                  # inla.spde.make.A(mesh)
                                  
                                  #-- Create index matrix --#
                                  field.indices =
                                    inla.spde.make.index("field", n.spde=mesh$n) 
                                  #Warning message: In inla.spde.make.index("field", n.mesh = mesh$n) :  
                                  #'n.mesh' is deprecated, please use 'n.spde' instead.
                                  
                                  
                                  #-- Create data stacks --#
                                  if(include.elevation){
                                    stack.est =
                                      inla.stack(data=list(rain=est.data),
                                                 A=list(A.est, 1),
                                                 effects=
                                                   list(c(field.indices,
                                                          list(Intercept=1)),
                                                        list(Elevation=est.elevation)),
                                                 tag="est")
                                    
                                    stack.val =
                                      inla.stack(data=list(rain=NA),
                                                 A=list(A.val,1),
                                                 effects=
                                                   list(c(field.indices,
                                                          list(Intercept=1)),
                                                        list(Elevation=val.elevation)),
                                                 tag="val")
                                    
                                    stack.pred =
                                      inla.stack(data=list(rain=NA),
                                                 A=list(A.pred),
                                                 effects=
                                                   list(c(field.indices,
                                                          list(Intercept=1))),
                                                 tag="pred")
                                    
                                    stack = inla.stack(stack.est, stack.val, stack.pred)
                                    
                                    
                                    #-- Define the formula --#
                                    formula <- rain ~ -1 + Intercept + Elevation + f(field, model=spde)
                                    
                                  }else{
                                    #no elevation
                                    stack.est =
                                      inla.stack(data=list(rain=est.data),
                                                 A=list(A.est),
                                                 effects=
                                                   list(c(field.indices,
                                                          list(Intercept=1))),
                                                 tag="est")
                                    # stack.val =
                                    #   inla.stack(data=list(rain=NA),
                                    #              A=list(A.val),
                                    #              effects=
                                    #                list(c(field.indices,
                                    #                       list(Intercept=1))),
                                    #              tag="val")
                                    # 
                                    # stack.pred =
                                    #   inla.stack(data=list(rain=NA),
                                    #              A=list(A.pred),
                                    #              effects=
                                    #                list(c(field.indices,
                                    #                       list(Intercept=1))),
                                    #              tag="pred")
                                    
                                    stack = inla.stack(stack.est)#, stack.val, stack.pred)
                                    
                                    #-- Define the formula --#
                                    formula <- rain ~ -1 + Intercept +  f(field, model=spde)
                                    
                                  }
                                  
                                  #-- Call INLA and get results --#
                                  INLA:::inla.dynload.workaround()
                                  mod =   inla(formula,
                                               data=inla.stack.data(stack, spde=spde),
                                               family="gaussian",
                                               control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
                                               control.compute=list(cpo=TRUE, dic=TRUE),
                                               # control.inla = list(h = .001),
                                               keep=FALSE, verbose=F)
                                  
                                  mod$dic$dic
                                  print(summary(mod))
                                  
                                  write.csv(c(co,storms_to_eval[co]),
                                            paste0("/home/walsh124/NAM-Model-Validation/csvtest/MLE_INLA/storm",storms_to_eval[co],".csv"))
                                  
                                  
                                  ###########################################################################
                                  #EXTRACT POSTERIOR SUMMARY STATISTICS
                                  ###########################################################################
                                  
                                  #-- Extract results for fixed effects - covariate coeffs --#
                                  beta = mod$summary.fixed[,"mean"]
                                  beta_sd = mod$summary.fixed[,"sd"]
                                  
                                  #-- Extract results for sigma2eps (1/precision)
                                  sigma2e_marg =
                                    inla.tmarginal(function(x) 1/x,
                                                   mod$marginals.hyperpar$"Precision for the Gaussian observations")
                                  sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
                                  sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
                                  sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
                                  sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
                                  cat("-----Results for sigma2eps-----\n",
                                      c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),
                                      "\n-----------------------------")
                                  
                                  #-- Extract results for the covariance function parameters --#
                                  mod.field = inla.spde2.result(mod, name="field", spde)
                                  
                                  var.nom.marg = mod.field$marginals.variance.nominal[[1]]
                                  var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
                                  var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
                                  var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
                                  var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
                                  cat("-----Results for sigma2omega-----\n",
                                      c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
                                      "\n-----------------------------")
                                  
                                  range.nom.marg = mod.field$marginals.range.nominal[[1]]
                                  range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
                                  range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
                                  range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
                                  range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
                                  cat("-----Results for the range-----\n",
                                      c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
                                      "\n-----------------------------")
                                  
                                  range_Ivec <- range.nom.m1 #_Ivec[i] for all 6 of these lines
                                  range_Ivec_MED <- range.nom.quantiles[2]
                                  sig2ep_Ivec <- sigma2e_m1
                                  sig2ep_Ivec_MED <-  sigma2e_quantiles[2]
                                  sig2om_Ivec <- var.nom.m1
                                  sig2om_Ivec_MED <-  var.nom.quantiles[2]
                                  
                                  theta1 <- mod$summary.hyperpar$`0.5quant`[2]
                                  theta1mode <- mod$summary.hyperpar$mode[2]
                                  theta2 <- mod$summary.hyperpar$`0.5quant`[3]
                                  theta2mode <- mod$summary.hyperpar$mode[3]
                                  tau_ <- exp(theta1)
                                  taumode_ <- exp(theta1mode)
                                  kappa_ <- exp(theta2)
                                  kappamode_ <- exp(theta2mode)
                                  (range_ <- sqrt(8)/kappa_)
                                  rangemode_ <- sqrt(8)/kappamode_
                                  (sig2_  <- 1/(4*pi*kappa_^2*tau_^2))
                                  (sig2mode_  <- 1/(4*pi*kappamode_^2*taumode_^2))
                                  (sig2eps_ <- 1/mod$summary.hyperpar$`0.5quant`[1])
                                  (sig2epsmode_ <- 1/mod$summary.hyperpar$mode[1])
                                  
                                  c(sig2eps_, sig2eps_-sig2ep_Ivec)
                                  c(sig2_, sig2_-sig2om_Ivec)
                                  c(range_, range_-range_Ivec)
                                  
                                  
                                  #Plotting posterior densities for range, variance, 
                                  par(mfrow=c(3,2))
                                  result = inla.spde2.result(mod, "field", spde)
                                  plot(result[["marginals.range.nominal"]][[1]], type = "l",
                                       main = paste(name_Ivec, "Nominal range, posterior density")) #name_Ivec[i]
                                  abline(v=range_, col="blue")
                                  abline(v=rangemode_, col="green")
                                  
                                  plot(result$marginals.variance.nominal[[1]], type="l",
                                       main="Nominal Variance, posterior density")
                                  abline(v=sig2_, col="blue")
                                  abline(v=sig2mode_, col="green")
                                  
                                  plot(mod$marginals.hy[[1]], type = "l", ylab = "Density", xlab = names(mod$marginals.hy)[1],
                                       main="Raw: Precision results")
                                  abline(v=1/sig2eps_, col="blue")
                                  abline(v=1/sig2epsmode_, col="green")
                                  
                                  plot(mod$marginals.hyperpar[[1]][mod$marginals.hy[[1]][,2] > 1e-10,], type="l",
                                       main="Raw: Precision (zoomed)")
                                  abline(v=1/sig2eps_, col="blue")
                                  abline(v=1/sig2epsmode_, col="green")
                                  
                                  plot(mod$marginals.hy[[2]], type = "l", ylab = "Density", xlab = names(mod$marginals.hy)[2],
                                       main="Raw: Theta1 results")
                                  abline(v=theta1, col="blue")
                                  abline(v=theta1mode, col="green")
                                  
                                  plot(mod$marginals.hy[[3]], type = "l", ylab = "Density", xlab = names(mod$marginals.hy)[3],
                                       main="Raw: Theta2 results")
                                  abline(v=theta2, col="blue")
                                  abline(v=theta2mode, col="green")
                                  #}
                                  list(name_Ivec, sig2om_Ivec, range_Ivec,sig2ep_Ivec, beta, 
                                       MLtemp$cov.pars[1], MLtemp$cov.pars[2], MLtemp$nugget, MLtemp$beta, MLtemp$kappa,
                                       MLtemp, as.geodata(hurricane), hur_opt$par, hur_opt$hessian, hur_opt)
                                }

stopCluster(cl)

if(exists("par_calcy")){save.image(file=paste0("~/NAM-Model-Validation/RData/RDatasim0/seed",seed,"/",
                                               if(storms_to_eval < 10){"0"},
                                               storms_to_eval,".RData"))}

par_calcy_vec <- cbind(unlist(par_calcy[[1]]), unlist(par_calcy[[2]]), unlist(par_calcy[[3]]),
                       unlist(par_calcy[[4]]), unlist(par_calcy[[5]]), unlist(par_calcy[[6]]),
                       unlist(par_calcy[[7]]), unlist(par_calcy[[8]]), unlist(par_calcy[[9]]),
                       unlist(par_calcy[[10]]))
par_calcy_vec2 <- par_calcy_vec
if(is.null(dim(par_calcy_vec))) par_calcy_vec <- t(par_calcy_vec)

colnames(par_calcy_vec) <- c("name_Ivec", "sig2om_Ivec", "range_Ivec","sig2ep_Ivec", "beta",
                             "MLEsigma2","MLEphi","MLEnugget","MLEbeta","MLEkappa")
# colnames(par_calcy) <- c("name_Ivec", "range_Ivec","sig2ep_Ivec","sig2om_Ivec",
#                          "MLEsigma2","MLEphi","MLEnugget","MLEkappa")
name_Ivec <- par_calcy_vec[,"name_Ivec"]
par_calcy_vec <- par_calcy_vec[,-which(colnames(par_calcy_vec) %in% c("name_Ivec"))]
ifelse(is.null(dim(par_calcy_vec)), 
       {par_calcy_vec <- as.numeric(par_calcy_vec); 
       names(par_calcy_vec) <- c("sig2om_Ivec", "range_Ivec","sig2ep_Ivec", "beta",
                                 "MLEsigma2","MLEphi","MLEnugget","MLEbeta","MLEkappa")
       par_calcy_vec <- t(par_calcy_vec)},
       par_calcy_vec <- apply(par_calcy_vec,2,as.numeric))
rownames(par_calcy_vec) <- name_Ivec

# if(all.equal(storms_to_eval,1:24)){rownames(par_calcy) <- as.character(read.csv(
#   "~/NAM-Model-Validation/csv/MLEestimates_estKappa/MLEvsINLAstorms124_20200123134241.csv")[,1])
# }
# if(all.equal(storms_to_eval,25:47)){rownames(par_calcy) <- as.character(read.csv(
#   "~/NAM-Model-Validation/csv/MLEestimates_estKappa/MLEvsINLAstorms2547_20200123134552.csv")[,1])
# }

svg.param.ests.errorINLA <- par_calcy_vec[,c("range_Ivec", "sig2ep_Ivec", "sig2om_Ivec")]
# svg.param.ests.errorINLA_MED <- cbind(name_Ivec, range_Ivec_MED,sig2ep_Ivec_MED,sig2om_Ivec_MED)
print(svg.param.ests.errorINLA)

# write.csv(svg.param.ests.errorINLA,
# "csv/INLAvsWLS/svg.param.ests.error_deg_MaternSVGestimatesINLA_center_scale2ndhalf.csv")
# write.csv(svg.param.ests.errorINLA_MED,
# "csv/INLAvsWLS/svg.param.ests.error_deg_MaternSVGestimatesINLA_MED_center_scale2ndhalf.csv")


# dev.off()


after <- Sys.time()
before; after


#par(mfrow=c(2,3))
#for (i in 1:6) {hist(par_calcy_vec[,i], main = colnames(par_calcy_vec)[i])}
#for (i in 1:6) {hist(log(par_calcy_vec[,i]), main = paste("log",colnames(par_calcy_vec)[i]))}
#par_calcy_vec[,"MLEnugget"]
#
#apply(par_calcy_vec, 2, mean)
#
## pdf(paste0("~/NAM-Model-Validation/pdf/MLEvsINLA/buffer_ngb_PWmean/MLEvsINLAstorms",storms_char,"_est_smooth_05_", nowtime(),".pdf"))
#par(mfrow=c(2,3))
#MLEnug      <- par_calcy_vec[,"MLEnugget"]
#MLEsig2     <- par_calcy_vec[,"MLEsigma2"]
#MLEphi      <- par_calcy_vec[,"MLEphi"]
#sig2ep_Ivec <- par_calcy_vec[,"sig2ep_Ivec"]
#sig2om_Ivec <- par_calcy_vec[,"sig2om_Ivec"]
#range_Ivec  <- par_calcy_vec[,"range_Ivec"]
#nug_bounds  <- c(min(MLEnug, sig2ep_Ivec), max(MLEnug, sig2ep_Ivec))
#phi_bounds  <- c(min(MLEphi, range_Ivec),max(MLEphi,range_Ivec))
#sig2_bounds <- c(min(MLEsig2, sig2om_Ivec),max(MLEsig2,sig2om_Ivec))

#hist(MLEnug, main="nug MLE")
#hist(MLEphi, main="phi MLE")
#hist(MLEsig2,main="sigma2 MLE")
#hist(sig2ep_Ivec, main="nug INLA")
#hist(range_Ivec*sqrt(8), main="phi INLA")
#hist(sig2om_Ivec,main="sigma2 INLA")
#
#hist(MLEnug, main="nug MLE", xlim = nug_bounds)
#hist(MLEphi, main="phi MLE",xlim = phi_bounds)
#hist(MLEsig2,main="sigma2 MLE",xlim = sig2_bounds)
#hist(sig2ep_Ivec, main="nug INLA", xlim = nug_bounds)
#hist(range_Ivec*sqrt(8), main="phi INLA")
#hist(sig2om_Ivec,main="sigma2 INLA")
#
#if(all(log(MLEnug)==-Inf)){MLEnug <- 1e-300}
#hist(log(MLEnug), main="log nug MLE")
#hist(log(MLEphi), main="log phi MLE")
#hist(log(MLEsig2),main="log sigma2 MLE")
#hist(log(sig2ep_Ivec), main="log nug INLA")
#hist(log(range_Ivec*sqrt(8)), main="log phi INLA")
#hist(log(sig2om_Ivec),main="log sigma2 INLA")
#
#hist(par_calcy_vec[,"MLEkappa"], main = "MLE kappa"); hist(log(par_calcy_vec[,"MLEkappa"]), main = "log MLE kappa")
# dev.off()


write.csv(par_calcy_vec,paste0("~/NAM-Model-Validation/csv/simsMLEout/seed",seed,"/sim_", if(storms_to_eval<10){"0"},storms_to_eval))

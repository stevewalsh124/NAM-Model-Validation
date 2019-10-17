# Based on code from Blangiardo
remove(list=ls())
require(geoR)
require(INLA)
require(gstat)
require(fields)
library(stringr)

#compare results with WLS approach, here alpha = 2, so nu = 1 in these WLS results:
# WLS_results_fixnu1 <- read.csv("csv/INLAvsWLS/svg.param.ests.error_KM_LCC_MaternSVGestimates.csv")

before <- Sys.time()

storms_to_eval <- 47
pdf("posterior_plots_storm47.pdf")

include.elevation = F #include/exclude elevation 
# I for INLA results
sig2om_Ivec <- c()
sig2ep_Ivec <- c()
range_Ivec <- c()
sig2om_Ivec_MED <- c()
sig2ep_Ivec_MED <- c()
range_Ivec_MED <- c()
name_Ivec <- c()

for(i in 1:length(list.files("csv/INLAvsWLSdeg/",full.names = T)[storms_to_eval])){
  # i <- 1
  storm_file <- list.files("csv/INLAvsWLSdeg/",full.names = T)[storms_to_eval[i]] #"csv/INLAvsWLS/error"
  name_start_INLA <- str_locate_all(pattern ='_', storm_file)[[1]][,1]
  name_Ivec[i] <- substr(storm_file,name_start_INLA[1]+1,name_start_INLA[2]-1)
  print(name_Ivec[i])
  
  hurricane <- read.csv(storm_file)
  est.coord <- hurricane[,3:4]
  center.coord <- est.coord
  center.coord$x <- center.coord$x - mean(center.coord$x)
  center.coord$y <- center.coord$y - mean(center.coord$y)
  est.coord <- center.coord
  scale.coord <- center.coord
  scale.coord$x <- scale.coord$x/max(abs(scale.coord$x))
  scale.coord$y <- scale.coord$y/max(abs(scale.coord$y))
  est.coord <- scale.coord
  est.data  <- hurricane[,2]
  
  ###########################################################################
  #SPDE/INLA 
  ###########################################################################
  #-- Create the mesh --#
  mesh =
    inla.mesh.create.helper(points=est.coord,
                            max.edge=c(40,100), min.angle=c(21,21))
  mesh$n
  
  #-- Plot the triangulation --#
  # par(mar=c(0,2,2,0))
  # plot(mesh)
  # lines(sic.borders,lwd=3)
  # points(est.coord,pch=19, col="red", cex=1.2)
  
  #-- Construct the SPDE object --#
  spde = inla.spde2.matern(mesh=mesh, alpha = 2) #alpha=2 by default
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
  mod =   inla(formula,
               data=inla.stack.data(stack, spde=spde),
               family="gaussian",
               control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
               control.compute=list(cpo=TRUE, dic=TRUE),
               keep=FALSE, verbose=TRUE)
  
  mod$dic$dic
  print(summary(mod))
  
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
  
  range_Ivec[i] <- range.nom.m1
  range_Ivec_MED[i] <- range.nom.quantiles[2]
  sig2ep_Ivec[i] <- sigma2e_m1
  sig2ep_Ivec_MED[i] <-  sigma2e_quantiles[2]
  sig2om_Ivec[i] <- var.nom.m1
  sig2om_Ivec_MED[i] <-  var.nom.quantiles[2]
  
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
       main = paste(name_Ivec[i], "Nominal range, posterior density"))
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
}

svg.param.ests.errorINLA <- cbind(name_Ivec, range_Ivec,sig2ep_Ivec,sig2om_Ivec)
svg.param.ests.errorINLA_MED <- cbind(name_Ivec, range_Ivec_MED,sig2ep_Ivec_MED,sig2om_Ivec_MED)
print(svg.param.ests.errorINLA)

# write.csv(svg.param.ests.errorINLA,
# "csv/INLAvsWLS/svg.param.ests.error_deg_MaternSVGestimatesINLA_center_scale2ndhalf.csv")
# write.csv(svg.param.ests.errorINLA_MED,
# "csv/INLAvsWLS/svg.param.ests.error_deg_MaternSVGestimatesINLA_MED_center_scale2ndhalf.csv")


dev.off()


after <- Sys.time()
before; after

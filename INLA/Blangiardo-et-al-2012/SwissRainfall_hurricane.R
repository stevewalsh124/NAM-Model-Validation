remove(list=ls())
require(geoR)
require(INLA)
require(gstat)
require(fields)
library(stringr)

#compare results with WLS approach, here alpha = 2, so nu = 1 in these WLS results:
WLS_results_fixnu1 <- read.csv("csv/INLAvsWLS/svg.param.ests.error_KM_LCC_MaternSVGestimates.csv")


include.elevation = F #include/exclude elevation 
# I for INLA results
sig2om_Ivec <- c()
sig2ep_Ivec <- c()
range_Ivec <- c()
sig2om_Ivec_MED <- c()
sig2ep_Ivec_MED <- c()
range_Ivec_MED <- c()
name_Ivec <- c()

for(i in 1:length(list.files("csv/INLAvsWLS/error",full.names = T))){
  # i <- 1
  storm_file <- list.files("csv/INLAvsWLS/error",full.names = T)[i]
  name_start_INLA <- str_locate_all(pattern ='_', storm_file)[[1]][,1]
  name_Ivec[i] <- substr(storm_file,name_start_INLA[1]+1,name_start_INLA[2]-1)
  print(name_Ivec[i])
  
  hurricane <- read.csv(storm_file)
  est.coord <- hurricane[,3:4]
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
  par(mar=c(0,2,2,0))
  plot(mesh)
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
}

svg.param.ests.errorINLA <- cbind(name_Ivec, range_Ivec,sig2ep_Ivec,sig2om_Ivec)
svg.param.ests.errorINLA_MED <- cbind(name_Ivec, range_Ivec_MED,sig2ep_Ivec_MED,sig2om_Ivec_MED)
print(svg.param.ests.errorINLA)

# write.csv(svg.param.ests.errorINLA, 
          # "csv/INLAvsWLS/svg.param.ests.error_KM_LCC_MaternSVGestimatesINLA.csv")
# write.csv(svg.param.ests.errorINLA_MED, 
          "csv/INLAvsWLS/svg.param.ests.error_KM_LCC_MaternSVGestimatesINLA_MED.csv")

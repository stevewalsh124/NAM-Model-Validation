# remove(list=ls()) 
require(geoR)
require(INLA)
require(gstat)
require(fields)
library(stringr)

#compare results with WLS approach, here alpha = 2, so nu = 1 in these WLS results:
WLS_results_fixnu1 <- read.csv("csv/INLAvsWLS/svg.param.ests.error_KM_LCC_MaternSVGestimates.csv")
sim_params <- read.csv("sim_params.csv")

iters <- 5

INLAmedians <- list()
INLAmeans <- list()

for (c in 1:length(sim_params)) {
  print(c)
  xaxis <- seq(sim_params$xmin[c], sim_params$xmax[c], sim_params$res[c])#seq(140,1500, by=35)#seq(-105,-90,.5)#seq(0,2,0.05)
  yaxis <- seq(sim_params$ymin[c], sim_params$ymax[c], sim_params$res[c])#seq(445,1930, by=35)#seq(24,35,.5)#seq(0,2,0.05)
  # xaxis <- (xaxis - mean(xaxis))/(max(xaxis)-min(xaxis))
  # yaxis <- (yaxis - mean(yaxis))/(max(yaxis)-min(yaxis))
    
  units <- "9sims"#0to1"#"KM_LCC"#"deg" #
  
  s <- matrix(NA,nrow=length(xaxis)*length(yaxis),ncol=2)
  for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) s[i+(j-1)*length(xaxis),] <- c(xaxis[i],yaxis[j])
  t <- as.matrix(dist(s))
  
  ##############################
  # Matern covariance function #
  ##############################
  
  # With and without nugget effect (see covmatrix below)
  
  sigma2 <- sim_params$sig2[c]#1#0.2#.2#0.5# 0.2 #sig2om
  phi <-    sim_params$phi[c]#98#1.75#2#4#1.75 #raange
  nu <- 1.0 #this is fixed in inla (alpha=2 is max, alpha = nu + d/2)
  tau2 <-   sim_params$tau2[c]#.05#.04#.1#.04 #sig2ep
  # plot(seq(0,1,0.001),sigma2*matern(seq(0,1,0.001),phi=(1/phi),kappa=nu),type="l")
  
  covmatrix <- diag(tau2,nrow(t)) + sigma2*matern(t,phi=(1/phi),kappa=nu) #nugget
  
  include.elevation = F #include/exclude elevation 
  # I for INLA results
  sig2om_Ivec <- c()
  sig2ep_Ivec <- c()
  range_Ivec <- c()
  sig2om_Ivec_MED <- c()
  sig2ep_Ivec_MED <- c()
  range_Ivec_MED <- c()
  name_Ivec <- c()
  
  for(i in 1:iters){
    print(c(c, i))
    # i <- 1
    # storm_file <- list.files("csv/INLAvsWLS/error",full.names = T)[i]
    # name_start_INLA <- str_locate_all(pattern ='_', storm_file)[[1]][,1]
    # name_Ivec[i] <- substr(storm_file,name_start_INLA[1]+1,name_start_INLA[2]-1)
    # print(name_Ivec[i])
    
    #covmatrix.sqrt already loaded from simulate_error_field.R for the WLS bootstrap
    x <- t(chol(covmatrix)) %*% rnorm(nrow(s))
    # z <- matrix(x,nrow=length(xaxis),ncol=length(yaxis),byrow=FALSE)
    # # persp(x=xaxis, y=yaxis, z=z)
    # image(z)
    
    # hurricane <- read.csv(storm_file)
    est.coord <- s#hurricane[,3:4]
    est.data  <- x#hurricane[,2]
    
    
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
    # plot(mod.field$marginals.range.nominal[[1]][,1], mod.field$marginals.range.nominal[[1]][,2],type="l") #post distbn
    range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
    range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
    range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
    range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
    cat("-----Results for the range-----\n",
        c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
        "\n-----------------------------")
    
    range_Ivec[i] <- sqrt(8)/range.nom.m1
    range_Ivec_MED[i] <- sqrt(8)/range.nom.quantiles[2]
    sig2ep_Ivec[i] <- sigma2e_m1
    sig2ep_Ivec_MED[i] <-  sigma2e_quantiles[2]
    sig2om_Ivec[i] <- var.nom.m1
    sig2om_Ivec_MED[i] <-  var.nom.quantiles[2]
  }
  
  svg.param.ests.errorINLA_SIM <- cbind(range_Ivec,sig2ep_Ivec,sig2om_Ivec) #name_Ivec,   
  svg.param.ests.errorINLA_MED_SIM <- cbind(range_Ivec_MED,sig2ep_Ivec_MED,sig2om_Ivec_MED) #name_Ivec, 
  print(svg.param.ests.errorINLA_SIM)
  INLAmedians[[c]] <- svg.param.ests.errorINLA_SIM
  INLAmeans[[c]] <-svg.param.ests.errorINLA_MED_SIM
  
  par(mar=c(5,4,4,1))
  par(mfrow=c(2,2))
  # Note the necessary transformation above for range: range = sqrt(8)/kappa => kappa=phi=sqrt(8)/range
  hist(svg.param.ests.errorINLA_SIM[,1], main=paste(c,"phi",phi,iters,"sims"), xlab = expression(phi)); abline(v=phi, col="green")
  hist(svg.param.ests.errorINLA_SIM[,2], main=paste(c,"tau2",tau2,iters,"sims"), xlab=expression(tau^2)); abline(v= tau2, col="red")
  hist(svg.param.ests.errorINLA_SIM[,3], main=paste(c,"sig2",sigma2,iters,"sims"), xlab=expression(sigma^2)); abline(v= sigma2, col="blue")
  hist(svg.param.ests.errorINLA_SIM[,3], main=paste(c,"sig2",sigma2,iters,"sims"), breaks=30, xlab=expression(sigma^2)); abline(v= sigma2, col="blue")
  
  cbind(mean(svg.param.ests.errorINLA_SIM[,1]), phi) #phi from range
  cbind(mean(svg.param.ests.errorINLA_SIM[,2]),tau2)         #sig2ep aka tau2 nugget
  cbind(mean(svg.param.ests.errorINLA_SIM[,3]),sigma2)         #sig2om aka sigma2 partial sill
  
  #plotting posterior distributions
  par(mfrow = c(2, 2))
  plot(mod$marginals.fix[[1]], type = "l", xlab = "Intercept", ylab = "Density", 
       main= "1 sim: intercept distbn"); abline(v=0, col="blue")
  plot(mod$marginals.hyp[[1]], type = "l", ylab = "Density", xlab = "range", main="1 sim range distbn"); abline(v=phi, col="green") #check what phi is...
  plot(mod$marginals.hyperpar[[2]], type="l")
  plot(mod$marginals.hyperpar[[3]], type="l")
  

  # write.csv(svg.param.ests.errorINLA_SIM,
  #           paste0("csv/INLAvsWLS/svg.param.ests.error_",
  #                  units, "_MaternSVGestimatesINLA_SIM_",
  #                  "nug",tau2,"phi",phi,"sigsq",sigma2,"sims",iters,".csv"))
  # write.csv(svg.param.ests.errorINLA_MED_SIM,
  #           paste0("csv/INLAvsWLS/svg.param.ests.error_",
  #                  units, "_MaternSVGestimatesINLA_MED_SIM_",
  #                  "nug",tau2,"phi",phi,"sigsq",sigma2,"sims",iters,".csv"))

  
  # # Trial with data simulated from a GAUSSIAN cov fn with true params shown below in vertical lines
  # tempy <- read.csv("csv/INLAvsWLS/svg.param.ests.error_KM_LCC_MaternSVGestimatesINLA_MED_SIM_harvey.csv")
  # hist(sqrt(8)/tempy$range_Ivec_MED); abline(v=1.763) #these were estimates from harvey
  # hist(tempy$sig2ep_Ivec_MED); abline(v=.036)
  # hist(tempy$sig2om_Ivec_MED); abline(v=.201)
  
  Sys.time()
}  

#INLAmedians
#INLAmeans
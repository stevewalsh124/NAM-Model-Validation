#WLS results when using KM (not degrees) and nu was fixed to be 1 (equivalent to alpha = nu + d/2 = 2 for inla) 
WLS_results_fixnu1 <- read.csv("csv/INLAvsWLS/svg.param.ests.error_KM_LCC_MaternSVGestimates.csv")

#INLA results when alpha fixed at 2
avg1 <- read.csv("svg.param.ests.error_KM_LCC_MaternSVGestimatesINLA.csv") #1-6
avg2 <- read.csv("svg.param.ests.error_KM_LCC_MaternSVGestimatesINLA8to40.csv") #8-40
med1 <- read.csv("svg.param.ests.error_KM_LCC_MaternSVGestimatesINLA_MED.csv")
med2 <- read.csv("svg.param.ests.error_KM_LCC_MaternSVGestimatesINLA_MED8to40.csv")
mix  <- read.csv("svg.param.ests.error_KM_LCC_MaternSVGestimatesINLAmix.csv")
mix2  <- read.csv("svg.param.ests.error_KM_LCC_MaternSVGestimatesINLAmix2.csv")
mixmed  <- read.csv("svg.param.ests.error_KM_LCC_MaternSVGestimatesINLA_MEDmix.csv")
mixmed2  <- read.csv("svg.param.ests.error_KM_LCC_MaternSVGestimatesINLA_MEDmix2.csv")

avg2 <- avg2[-c(1:7),] # remove NA rows
med2 <- med2[-c(1:7),] 

inla_avg <- rbind(avg1,avg2) #just missing no 7
inla_med <- rbind(med1,med2)


#WLS for degrees (gaussian and nu set to 1)
gauss <- read.csv("svg.param.ests.error_deg_ARC_GAUSS_deg.csv")
matern <- read.csv("svg.param.ests.error_deg_ARC_MATERN_deg.csv")


#INLA for degrees: orig, center, scale
inla_deg_orig <- read.csv("csv/INLAvsWLS/svg.param.ests.error_deg_MaternSVGestimatesINLA.csv")
inla_deg_center1 <- read.csv("csv/INLAvsWLS/svg.param.ests.error_deg_MaternSVGestimatesINLA_center.csv")
inla_deg_center2 <- read.csv("csv/INLAvsWLS/svg.param.ests.error_deg_MaternSVGestimatesINLA_center2ndhalf.csv")
inla_deg_center <- rbind(inla_deg_center1, inla_deg_center2[30:47,])
inla_deg_scale1 <- read.csv("csv/INLAvsWLS/svg.param.ests.error_deg_MaternSVGestimatesINLA_center_scale.csv")
inla_deg_scale2 <- read.csv("csv/INLAvsWLS/svg.param.ests.error_deg_MaternSVGestimatesINLA_center_scale2ndhalf.csv")
inla_deg_scale <- rbind(inla_deg_scale1, inla_deg_scale2[42:47,])

par(mar=c(5,4,4,1))
plot(inla_deg_orig$range_Ivec, inla_deg_scale$range_Ivec, main="INLA range\n Orig Vs Scale")
lm(inla_deg_orig$range_Ivec~ inla_deg_scale$range_Ivec)

plot(inla_deg_orig$sig2ep_Ivec, inla_deg_scale$sig2ep_Ivec, main="INLA tau^2\n Orig Vs Scale")
lm(inla_deg_orig$sig2ep_Ivec~ inla_deg_scale$sig2ep_Ivec)

plot(inla_deg_orig$sig2om_Ivec, inla_deg_scale$sig2om_Ivec, main="INLA sig^2\n Orig Vs Scale")
lm(inla_deg_orig$sig2om_Ivec~ inla_deg_scale$sig2om_Ivec)



plot(inla_deg_orig$range_Ivec, inla_deg_center$range_Ivec, main="INLA range\n Orig Vs Center")
lm(inla_deg_orig$range_Ivec~ inla_deg_center$range_Ivec)

plot(inla_deg_orig$sig2ep_Ivec, inla_deg_center$sig2ep_Ivec, main="INLA tau^2\n Orig Vs Center")
lm(inla_deg_orig$sig2ep_Ivec~ inla_deg_center$sig2ep_Ivec)

plot(inla_deg_orig$sig2om_Ivec, inla_deg_center$sig2om_Ivec, main="INLA sig^2\n Orig Vs Center")
lm(inla_deg_orig$sig2om_Ivec~ inla_deg_center$sig2om_Ivec)




#all WLS comparisons
par(mfrow=c(2,2))
par(mar=c(5,4,4,1))
plot(gauss$phivec,matern$phivec , main="WLS: Gauss Vs Matern")
plot(gauss$prRangevec, matern$prRangevec)
plot(gauss$tau2vec, matern$tau2vec)
plot(gauss$sig2vec, matern$sig2vec)

cor(gauss$phivec,matern$phivec )
cor(gauss$prRangevec, matern$prRangevec)
cor(gauss$tau2vec, matern$tau2vec)
cor(gauss$sig2vec, matern$sig2vec)

#inla avg vs inla med
par(mfrow=c(1,3))
plot(inla_avg$range_Ivec, inla_med$range_Ivec_MED , main="INLA avg vs med")
plot(inla_avg$sig2ep_Ivec, inla_med$sig2ep_Ivec_MED)
plot(inla_avg$sig2om_Ivec, inla_med$sig2om_Ivec_MED)

cor(inla_avg$range_Ivec, inla_med$range_Ivec_MED )
cor(inla_avg$sig2ep_Ivec, inla_med$sig2ep_Ivec_MED)
cor(inla_avg$sig2om_Ivec, inla_med$sig2om_Ivec_MED)

WLS_results_fixnu1[-7,"phivec"]
par(mfrow=c(3,1))
par(mar=c(5,4,4,1))
plot(WLS_results_fixnu1[-c(7,22),"phivec"], inla_avg$range_Ivec[-22], main="INLA vs WLS")
plot(WLS_results_fixnu1[-c(7,22),"tau2vec"], inla_avg$sig2ep_Ivec[-22])
plot(WLS_results_fixnu1[-c(7,22),"sig2vec"], inla_avg$sig2om_Ivec[-22])


cor(WLS_results_fixnu1[-c(7,22),"phivec"], inla_avg$range_Ivec[-22])
cor(WLS_results_fixnu1[-c(7,22),"tau2vec"], inla_avg$sig2ep_Ivec[-22])
cor(WLS_results_fixnu1[-c(7,22),"sig2vec"], inla_avg$sig2om_Ivec[-22])

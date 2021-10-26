# log scoring rule

# let y = observations = ST4
# D = NAM for that storm, and all ST4, NAM pairs of previous storms
# theta is spatial parameters

# consider P(y|D) = \int p(y|D, theta) * pi(theta | D) dtheta
# we have MCMC draws for pi(theta | D), so approximate with
# 1/G * sum_{g=1}^G p(y | D, theta_g)

# by equation 2 in our paper, p(y | D, theta_g) implies
# p(ST4_new | NAM_new, all old NAM&ST4, theta_g) = MVN(mu + NAM_new, Sigma(theta_g)))

library(plgp)
library(raster)
library(mvtnorm)
library(LaplacesDemon)
library(scoringRules)

#informative prior yields informed posterior (based on geq20 params)
PW_post <- raster("~/NAM-Model-Validation/error_rasters_summary_sqrt/PW_post_informPrior.grd")

# storm to evaluate; 1-6for 2018 and 2019 storms
ste <- 5

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

# pdf(paste0("~/NAM-Model-Validation/pdf/logScoring_storm",ste,".pdf"))
load(paste0("~/NAM-Model-Validation/RData/prediction",ste))

dim(theta_pred) # 1000 x 2
G <- 100#nrow(theta_pred)

# error field = ST4 - NAM - mu has mean 0, so ST4 has mean NAM + mu
# NAM_r already has PWmean subtracted, so add it twice to have NAM + mu
NAM_plus_PW <- rasterToPoints(NAM_r + PW_mean + PW_mean)[,3]
NAM_plus_PWpost <- rasterToPoints(NAM_r + PW_mean + PW_post)[,3]
NAM_no_PW <- rasterToPoints(NAM_r + PW_mean )[,3]

my_dist_mtx <- sqrt(plgp::distance(X1 = coords))

# bootstrap samples
B <- 1000
theta_boot <- matrix(NA, B, 2)
for (b in 1:B) { theta_boot[b,] <- theta_hat[sample(1:nrow(theta_hat),1),] }

# weighted average of thetas based on Hessian info
sum_wi_ti <- c(0,0)
for (i in 1:nrow(theta_hat)) sum_wi_ti <- sum_wi_ti + hessians[[i]] %*% theta_hat[i,]
sum_prec <- Reduce("+", hessians)
theta_wgt <- solve(sum_prec) %*% sum_wi_ti
  
# option 1: theta = theta bar, avg of 47 MLEs
# option 2: theta ~ MVN(thetabar, solve(sum(H_i))) # sum of precision matrices
# option 3: theta ~ MVN(Bx_i, Sigma_theta) # the hierarchical model we employ
# option 4: theta ~ generate from bootstrap of 47 MLEs
theta1 <- theta_wgt
theta2 <- rmvn(1000, theta_wgt, solve(sum_prec))
theta3 <- theta_pred
theta4 <- theta_boot

log_score <- function(thet, dist_mtx, obsvn, mean_fn, G=(length(thet)/2)){
  thet <- matrix(thet, length(thet)/2, 2)
  log_pred_dens <- CRPSs <- c()
  for(g in 1:G){
    if(g %% 500 ==0) print(g)
    sig2 = exp(thet[g,2])
    phi = exp(thet[g,2]-thet[g,1])
    my_cov_mtx <- sig2 * exp(-0.5 * my_dist_mtx / phi)
    log_pred_dens[g] <-  dmvnorm(x = obsvn, mean = mean_fn, 
                                 sigma = my_cov_mtx, checkSymmetry = F, log = T)
    # CRPSs[g] <- crps_norm(y = obsvn, mean = mean_fn, sd = chol(my_cov_mtx))
    if(length(thet)==2) break
  }
  return(log_pred_dens)
}

# Obtain log score from each of the 4 sampling schemes for theta
out1 <- log_score(theta1, my_dist_mtx, ST4_pred$value, NAM_plus_PW)
out2 <- log_score(theta2, my_dist_mtx, ST4_pred$value, NAM_plus_PW)
out3 <- log_score(theta3, my_dist_mtx, ST4_pred$value, NAM_plus_PW)
out4 <- log_score(theta4, my_dist_mtx, ST4_pred$value, NAM_plus_PW)

out1n <- log_score(theta1, my_dist_mtx, ST4_pred$value, NAM_no_PW)
out2n <- log_score(theta2, my_dist_mtx, ST4_pred$value, NAM_no_PW)
out3n <- log_score(theta3, my_dist_mtx, ST4_pred$value, NAM_no_PW)
out4n <- log_score(theta4, my_dist_mtx, ST4_pred$value, NAM_no_PW)

out1p <- log_score(theta1, my_dist_mtx, ST4_pred$value, NAM_plus_PWpost)
out2p <- log_score(theta2, my_dist_mtx, ST4_pred$value, NAM_plus_PWpost)
out3p <- log_score(theta3, my_dist_mtx, ST4_pred$value, NAM_plus_PWpost)
out4p <- log_score(theta4, my_dist_mtx, ST4_pred$value, NAM_plus_PWpost)

log(mean(exp(out1 - max(out1)))) + max(out1)
log(mean(exp(out2 - max(out2)))) + max(out2)
log(mean(exp(out3 - max(out3)))) + max(out3)
log(mean(exp(out4 - max(out4)))) + max(out4)

log(mean(exp(out1n - max(out1n)))) + max(out1n)
log(mean(exp(out2n - max(out2n)))) + max(out2n)
log(mean(exp(out3n - max(out3n)))) + max(out3n)
log(mean(exp(out4n - max(out4n)))) + max(out4n)

log(mean(exp(out1p - max(out1p)))) + max(out1p)
log(mean(exp(out2p - max(out2p)))) + max(out2p)
log(mean(exp(out3p - max(out3p)))) + max(out3p)
log(mean(exp(out4p - max(out4p)))) + max(out4p)

# write log scores (ie, log predictive densities) to csv's
write.csv(out1, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta1_storm",ste,".csv"))
write.csv(out2, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta2_storm",ste,".csv"))
write.csv(out3, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta3_storm",ste,".csv"))
write.csv(out4, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta4_storm",ste,".csv"))

write.csv(out1n, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta1_storm",ste,"_noPW.csv"))
write.csv(out2n, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta2_storm",ste,"_noPW.csv"))
write.csv(out3n, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta3_storm",ste,"_noPW.csv"))
write.csv(out4n, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta4_storm",ste,"_noPW.csv"))

write.csv(out1p, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta1_storm",ste,"_PWpost.csv"))
write.csv(out2p, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta2_storm",ste,"_PWpost.csv"))
write.csv(out3p, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta3_storm",ste,"_PWpost.csv"))
write.csv(out4p, file=paste0("~/NAM-Model-Validation/csv/scores/logS_theta4_storm",ste,"_PWpost.csv"))

# # plots
# par(mfrow=c(2,2))
#   p1 <- hist(out1, prob=T, main = paste0("storm ", ste, ": logS for theta1\n mean: ",
#                                          round(log(mean(exp(out1 - max(out1)))) + max(out1),4)))
#   p2 <- hist(out2, prob=T, main = paste0("storm ", ste, ": logS for theta2\n mean: ",
#                                          round(log(mean(exp(out2 - max(out2)))) + max(out2),4)))
#   p3 <- hist(out3, prob=T, main = paste0("storm ", ste, ": logS for theta3\n mean: ",
#                                          round(log(mean(exp(out3 - max(out3)))) + max(out3),4)))
#   p4 <- hist(out4, prob=T, main = paste0("storm ", ste, ": logS for theta4\n mean: ",
#                                          round(log(mean(exp(out4 - max(out4)))) + max(out4),4)))
# 
# par(mfrow=c(1,1))
# plot( p4, col=rgb(1,0,0,1/4), xlim=range(out1,out2,out3,out4),
#       main = paste0("Storm ", ste, ", Some of HM tail cut off")) # first histogram
# plot( p2, col=rgb(0,1,0,1/4), xlim=range(out1,out2,out3,out4), add=T)  # second
# plot( p1, col=rgb(0,0,1,1/4), xlim=range(out1,out2,out3,out4), add=T)  # second
# plot( p3, col=rgb(0,0,1,1/4), xlim=range(out1,out2,out3,out4), add=T)  # second
# legend("topleft",legend = c("thetabar (2 schemes)","HM","Boot"),
#        col=c("black", "blue","purple"), lty=1)
# abline(v=log(mean(exp(out1 - max(out1)))) + max(out1), col="black", lty=1, lwd=2)
# abline(v=log(mean(exp(out2 - max(out2)))) + max(out2), col="black", lty=2, lwd=2)
# abline(v=log(mean(exp(out3 - max(out3)))) + max(out3), col="blue", lty=1, lwd=5)
# abline(v=log(mean(exp(out4 - max(out4)))) + max(out4), col="purple", lty=4, lwd=5)
# 
# dev.off()

# # combine all of the storms into one pdf below
# boot_better <- c()
# pdf("~/NAM-Model-Validation/pdf/logScoring_allstorms.pdf")
# for (s in 1:6) {
#   out1 <- read.csv(paste0("~/NAM-Model-Validation/csv/scores/logS_theta1_storm",s,".csv"), row.names = 1)$x
#   out2 <- read.csv(paste0("~/NAM-Model-Validation/csv/scores/logS_theta2_storm",s,".csv"), row.names = 1)$x
#   out3 <- read.csv(paste0("~/NAM-Model-Validation/csv/scores/logS_theta3_storm",s,".csv"), row.names = 1)$x
#   out4 <- read.csv(paste0("~/NAM-Model-Validation/csv/scores/logS_theta4_storm",s,".csv"), row.names = 1)$x
# 
#   par(mfrow=c(2,2))
#   p1 <- hist(out1, prob=T, main = paste0("storm ", s, ": logS for theta1\n mean: ",
#                                          round(log(mean(exp(out1 - max(out1)))) + max(out1),4)))
#   p2 <- hist(out2, prob=T, main = paste0("storm ", s, ": logS for theta2\n mean: ",
#                                          round(log(mean(exp(out2 - max(out2)))) + max(out2),4)))
#   p3 <- hist(out3, prob=T, main = paste0("storm ", s, ": logS for theta3\n mean: ",
#                                          round(log(mean(exp(out3 - max(out3)))) + max(out3),4)))
#   p4 <- hist(out4, prob=T, main = paste0("storm ", s, ": logS for theta4\n mean: ",
#                                          round(log(mean(exp(out4 - max(out4)))) + max(out4),4)))
# 
#   par(mfrow=c(1,1))
#   plot( p4, col=rgb(1,0,0,1/4), xlim=range(out1,out2,out3,out4),
#         main = paste0("Storm ", s)) # first histogram
#   plot( p2, col=rgb(0,1,0,1/4), xlim=range(out1,out2,out3,out4), add=T)  # second
#   plot( p1, col=rgb(0,0,1,1/4), xlim=range(out1,out2,out3,out4), add=T)  # second
#   plot( p3, col=rgb(0,0,1,1/4), xlim=range(out1,out2,out3,out4), add=T)  # second
#   legend("topleft",legend = c("thetabar (2 schemes)","HM","Boot"),
#          col=c("black", "blue","purple"), lty=1)
#   abline(v=log(mean(exp(out1 - max(out1)))) + max(out1), col="black", lty=1, lwd=2)
#   abline(v=log(mean(exp(out2 - max(out2)))) + max(out2), col="black", lty=2, lwd=2)
#   abline(v=log(mean(exp(out3 - max(out3)))) + max(out3), col="blue", lty=1, lwd=5)
#   abline(v=log(mean(exp(out4 - max(out4)))) + max(out4), col="purple", lty=4, lwd=5)
#   boot_better[s] <- log(mean(exp(out4 - max(out4)))) + max(out4) > log(mean(exp(out3 - max(out3)))) + max(out3)
# }
# dev.off()

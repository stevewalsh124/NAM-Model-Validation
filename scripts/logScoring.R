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

# storm to evaluate; 1-6for 2018 and 2019 storms
ste <- 6

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

load(paste0("~/NAM-Model-Validation/RData/prediction",ste))

dim(theta_pred) # 1000 x 2
G <- nrow(theta_pred)

# error field = ST4 - NAM - mu has mean 0, so ST4 has mean NAM + mu
# NAM_r already has PWmean subtracted, so add it twice to have NAM + mu
NAM_plus_PW <- rasterToPoints(NAM_r + PW_mean + PW_mean)[,3]

my_dist_mtx <- sqrt(plgp::distance(X1 = coords))

log_pred_dens <- log_pred_straw <- c()

for(g in 1:G){
  if(g %% 500 ==0) print(g)
  sig2 = exp(theta_pred[g,2])
  phi = exp(theta_pred[g,2]-theta_pred[g,1])
  my_cov_mtx <- sig2 * exp(-0.5 * my_dist_mtx / phi)
  log_pred_dens[g] <-  dmvnorm(x = ST4_pred$value, mean = NAM_plus_PW, 
                               sigma = my_cov_mtx, checkSymmetry = F, log = T)
}

sig2 <- exp(theta_bar[2])
phi <- exp(theta_bar[2]-theta_bar[1])
my_cov_mtx <- sig2 * exp(-0.5 * my_dist_mtx / phi)
log_pred_straw <-  dmvnorm(x = ST4_pred$value, mean = NAM_plus_PW, 
                           sigma = my_cov_mtx, checkSymmetry = F, log = T)

mean(log_pred_dens)
log_pred_straw

# hist(log_pred_dens)
# abline(v=log_pred_straw)
# abline(v=mean(log_pred_dens), col="blue")

write.csv(log_pred_dens, file=paste0("~/NAM-Model-Validation/csv/log_score_",s))
write.csv(log_pred_straw, file=paste0("~/NAM-Model-Validation/csv/log_straw_",s))

# pdf("~/NAM-Model-Validation/pdf/logScoring_hierarchical_vs_thetabar.pdf")
# for (ste in 1:6) {
#   log_pred <- read.csv(paste0("~/NAM-Model-Validation/csv/log_score_",ste))[,2]
#   log_straw <- as.numeric(read.csv(paste0("~/NAM-Model-Validation/csv/log_straw_",ste))[2])
#   
#   hist(log_pred, main = paste(ste,"mean of hierarch. is blue"))
#   abline(v=log_straw)
#   abline(v=mean(log_pred), col="blue")
#   
# }
# dev.off()

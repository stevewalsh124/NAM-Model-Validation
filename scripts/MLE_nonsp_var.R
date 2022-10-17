# get MLEs of variance in the non-spatial setting

library(pracma) # ::hessian
library(LaplacesDemon) #:: dmvn

if(!dir.exists("csv/myMLEresults/nonsp/MLEs/")) {dir.create("csv/myMLEresults/nonsp/MLEs/", recursive = T)}
if(!dir.exists("csv/myMLEresults/nonsp/hessians/")) {dir.create("csv/myMLEresults/nonsp/hessians/", recursive = T)}

# Error fields for each of the 47 training storms
EF_files <- list.files("csv/error_df_sqrt/subtractPWmeanF/", full.names = T)
N <- length(EF_files)

# Store the MLEs for variance when mean is assumed to be 0 (and not assuming 0)
my_nonsp_vars <- matrix(NA, N, 2)
means <- hessies <- n_pix <- c()

ste <- 11

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

# Find each MLE, assuming 0 mean and not assuming 0 mean
for (i in ste) {
  print(i)
  EF <- read.csv(EF_files[i], row.names = 1)$value
  n_pix[i] <- length(EF)
  means[i] <- mean(EF)
  my_nonsp_vars[i,1] <- sum(EF^2)/length(EF)
  # my_nonsp_vars[i,2] <- sum((EF-mean(EF))^2)/length(EF)
  
  hessies[i] <- pracma::hessian(f = function(sig2){dmvn(EF, 0, diag(sig2, nrow = length(EF)),
                                                        log = T)}, x0 = sum(EF^2)/length(EF))
  
  write.csv(my_nonsp_vars[i,1], file=paste0("csv/myMLEresults/nonsp/MLEs/",if(i<10){"0"},i,".csv"))
  # write.csv(my_nonsp_vars[i,2], file=paste0("csv/myMLEs_nonsp/MLEs_nonzeromean/",if(i<10){"0"},i,".csv"))
  write.csv(hessies[i], file=paste0("csv/myMLEresults/nonsp/hessians/",if(i<10){"0"},i,".csv"))
}

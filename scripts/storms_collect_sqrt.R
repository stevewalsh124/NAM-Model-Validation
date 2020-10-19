# After running storm_MLE_Hessian.R for all 47 storms, use this
# code to collect all 47 results into a single data frame
# This will be input into the Gibbs sampler

path <- "/home/walsh124/NAM-Model-Validation/RData/RDatafixnug0/"
data_files <- list.files(path, pattern = "sqrt_2020101", full.names = T)

all_storm_res <- matrix(NA, length(data_files), 5)
rownamez <- lik_vals <- c()
hess_opt <- hess_lik <- list()
for (pe in 1:length(data_files)) {
  # pe <- 1
  load(data_files[pe])
  # rm(list=setdiff(ls(), c("par_calcy", "all_storm_res", "pe", "data_files","rownamez", "hess_opt", "hess_lik")))
  
  hess_lik[[pe]] <- par_calcy[[11]][[1]]$info.minimisation.function$hessian#^-1
  hess_opt[[pe]] <- par_calcy[[14]][[1]]#^-1 
  
  lik_vals[pe] <- par_calcy[[11]][[1]]$info.minimisation.function$value
  
  par_calcy_vec <- cbind(unlist(par_calcy[[1]]), unlist(par_calcy[[2]]), unlist(par_calcy[[3]]),
                         unlist(par_calcy[[4]]), unlist(par_calcy[[5]]), unlist(par_calcy[[6]]),
                         unlist(par_calcy[[7]]), unlist(par_calcy[[8]]), unlist(par_calcy[[9]]),
                         unlist(par_calcy[[10]]))
  par_calcy_vec2 <- par_calcy_vec
  if(is.null(dim(par_calcy_vec))) par_calcy_vec <- t(par_calcy_vec)
  
  colnames(par_calcy_vec) <- c("name_Ivec","MLEsigma2","MLEphi","MLEnugget","MLEbeta","MLEkappa")
  
  name_Ivec <- par_calcy_vec[,"name_Ivec"]
  par_calcy_vec <- par_calcy_vec[,-which(colnames(par_calcy_vec) %in% c("name_Ivec"))]
  ifelse(is.null(dim(par_calcy_vec)), 
         {par_calcy_vec <- as.numeric(par_calcy_vec); 
         names(par_calcy_vec) <- c("MLEsigma2","MLEphi","MLEnugget","MLEbeta","MLEkappa")
         par_calcy_vec <- t(par_calcy_vec)},
         par_calcy_vec <- apply(par_calcy_vec, 2, as.numeric))
  rownamez[pe] <- name_Ivec
  all_storm_res[pe,] <- par_calcy_vec
}

rownames(all_storm_res) <- rownamez
colnames(all_storm_res) <- c("MLEsigma2","MLEphi","MLEnugget","MLEbeta","MLEkappa")

# all_storm_res <- all_storm_res[order(rownames(all_storm_res)),]
# write.csv(all_storm_res, "~/NAM-Model-Validation/csv/all_storm_res.csv")

old_fixnug <- read.csv("~/NAM-Model-Validation/csv/MLEestimates_buffer_ngb_PWmean/MLEvsINLAstorms_allbut43_nugfix0.csv", row.names = 1)
allstormnames <- rownames(old_fixnug)
old_fixnug <- old_fixnug[rownames(all_storm_res),]

par(mfrow=c(2,2))
for (q in c("MLEsigma2","MLEphi","MLEbeta","MLEkappa")) {
  plot(all_storm_res[,q], old_fixnug[rownames(all_storm_res),q], main=q,
       xlab = NA, ylab=NA)
  print(lm(all_storm_res[,q]~ old_fixnug[rownames(all_storm_res),q]))
  abline(0,1)
  abline(lm(old_fixnug[rownames(all_storm_res),q]~all_storm_res[,q]), col="blue")
  
}

par(mfrow=c(2,3))
for (q in c(1:2,5)) {hist(all_storm_res[,q],main=colnames(all_storm_res)[q])}
for (q in c(1:2,5)) {hist(log(all_storm_res[,q]),main=paste("log",colnames(all_storm_res)[q]))}

for (k in 1:length(data_files)) {
  print(c(rownames(all_storm_res)[k],sqrt(diag(hess_opt[[k]]))))
}

for (k in 1:length(data_files)) {
  print(c(rownames(all_storm_res)[k],sqrt(diag(hess_lik[[k]]))))
}

hess_opt_diagsqrt <- matrix(NA, nrow=length(data_files),3)
for(i in 1:length(data_files)) hess_opt_diagsqrt[i,] <- (sqrt(diag(hess_opt[[i]])))


hess_lik_diagsqrt <- matrix(NA, nrow=length(data_files), ncol = 2)
for(i in 1:length(data_files)) hess_lik_diagsqrt[i,] <- (sqrt(diag(hess_lik[[i]])))



fileza <- list.files("~/NAM-Model-Validation/INLAvsWLSdeg_minusPWmeanNA/", full.names = T)
for (j in 1:length(fileza)) {
  print(c(j,dim(read.csv(fileza[j]))[1]))
}

hist(hess_opt_diagsqrt[,1])
hist(hess_opt_diagsqrt[,2])
hist(hess_opt_diagsqrt[,3])

all_storm_res
old_fixnug

hess_opt

avail<-c(); for(i in 1:47) avail[i] <- (i %in% as.numeric(substr(data_files,nchar(path) + 8, nchar(path) + 9)))
numvec <- 1:47
numvec[!avail]

locs <- read.csv("~/NAM-Model-Validation/csv/storm_levels_and_locs.csv")
rownames(locs) <- sort(c(allstormnames, "2016matthew"))
dim(locs[avail,])

# save(all_storm_res, file="~/NAM-Model-Validation/RData/all_storm_res.RData")
# save(hess_opt, file="~/NAM-Model-Validation/RData/hess_opt.RData")

# est0 <- all_storm_res
# fix0 <- all_storm_res
# fix001 <- all_storm_res
# fix0t <- fix0[rownames(fix001),]
# all_fix <- rbind(fix0t[,c(5,6,7,9)], fix001[,c(5,6,7,9)])

# plot(fix0[,"MLEsigma2"], est0[,"MLEsigma2"])
# merger0 <- merge(fix0,est0,by="X")
# mergerfix <- merge(fix0, fix001, by="X")
# 
# plot(merger0$MLEsigma2.x, merger0$MLEsigma2.y, main="nugs0 sigma2", xlab = "fix0", ylab = "est0"); abline(0,1)
# plot(merger0$MLEphi.x,    merger0$MLEphi.y,    main="nugs0 phi",    xlab = "fix0", ylab = "est0"); abline(0,1)
# plot(merger0$MLEkappa.x,  merger0$MLEkappa.y,  main="nugs0 kappa",  xlab = "fix0", ylab = "est0"); abline(0,1)
# 
# plot(mergerfix$MLEsigma2.x, mergerfix$MLEsigma2.y, main="fixnugs sigma2", xlab = "fix0", ylab = "fix001"); abline(0,1)
# plot(mergerfix$MLEphi.x,    mergerfix$MLEphi.y,    main="fixnugs phi",    xlab = "fix0", ylab = "fix001"); abline(0,1)
# plot(mergerfix$MLEkappa.x,  mergerfix$MLEkappa.y,  main="fixnugs kappa",  xlab = "fix0", ylab = "fix001"); abline(0,1)


## Compare fixnug0 with fixnug001
## run each separately and get all_fix above to run below

# # pdf("~/NAM-Model-Validation/pdf/nug0vs001.pdf")
# hist(all_storm_res[,"MLEnugget"], breaks=20)
# summary(all_storm_res[,"MLEnugget"])
# 
# pairs(all_fix, col = ifelse(all_fix[,"MLEnugget"]==0,"green","red"), cex=1.5)
# par(mfrow=c(3,2))
# for (i in 1:4) {
#   for (j in 1:4) {
#     if (j > i) {
#       plot(all_fix[,c(i,j)], col=ifelse(all_fix[,3]==0,"green","red"))
#       segments(all_fix[1:39,i], all_fix[1:39,j],all_fix[40:78,i], all_fix[40:78,j],)
#     }
#   }
# }
# 
# par(mfrow=c(2,2))
# for (i in 1:ncol(all_fix)) {
#   hist(all_fix[1:39,i]- all_fix[40:78,i], main=paste(colnames(all_fix)[i],"differences"),
#        xlab = paste(colnames(all_fix)[i],"differences"))
# }
# 
# summary(all_fix[1:39,1]- all_fix[40:78,1])
# # dev.off()



# Compare MLEs before and after the change in PWmean map (all24 vs 12/24/12 and the Maine bug)
# newer_res <- all_storm_res

if(exists("flat_storm_res")){
  both <- merge(flat_storm_res, all_storm_res, by="row.names", all=T)
  
  order(-1*(abs(both$MLEsigma2.x - both$MLEsigma2.y)))
  order(-1*(abs(both$MLEphi.x - both$MLEphi.y)))
  order(-1*(abs(both$MLEkappa.x - both$MLEkappa.y)))
  
  order(-1*(abs(both$MLEkappa.x - both$MLEkappa.y) + 
              abs(both$MLEphi.x - both$MLEphi.y) + 
              abs(both$MLEsigma2.x - both$MLEsigma2.y)))
  rownames(all_storm_res)[order(-1*(abs(both$MLEkappa.x - both$MLEkappa.y) + 
                                      abs(both$MLEphi.x - both$MLEphi.y) + 
                                      abs(both$MLEsigma2.x - both$MLEsigma2.y)))]
  
  library(plotly)
  plot_ly(data=both,x=~MLEkappa.x, y=~MLEkappa.y) %>% 
    add_annotations(x = both$MLEkappa.x, y = both$MLEkappa.y, text = rownames(both), showarrow = FALSE)
  
  plot_ly(data=both,x=~MLEphi.x, y=~MLEphi.y) %>% 
    add_annotations(x = both$MLEphi.x, y = both$MLEphi.y, text = rownames(both), showarrow = FALSE)
  
  plot_ly(data=both,x=~MLEsigma2.x, y=~MLEsigma2.y) %>% 
    add_annotations(x = both$MLEsigma2.x, y = both$MLEsigma2.y, text = rownames(both), showarrow = FALSE)
}

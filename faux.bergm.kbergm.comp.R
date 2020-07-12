#########################################
### This script is written for the comparison between bergm and kbergm given faux.mesa.high network
#Load required libraries and such
wants <- c("MASS", "parallel", "data.table", "parallel", 
           "dplyr", "stringr", "statnet", 
           "ergm", "foreach",
           "np", "mvtnorm", "network",
           "doParallel", "foreach","Bergm", "Matrix",
           "sna", "moments", "coda",  "ggplot2") 
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, library, character.only=TRUE)
############################################
data("faux.mesa.high")
############################################
faux_mod <- readRDS("faux_mod3.rds")
faux_mod_mple <- ergm(faux_mod$formula, estimate = "MPLE")
faux_kbergm_list <- readRDS("faux_kbergm_list_4442_24000_96000_sqrt1.rds")
faux_bergm_list_truth <- readRDS("faux_bergm_list_20505_sqrt5_super_high_quality_long.rds")
faux_bergm_list_standard <- readRDS("faux_bergm_list_standard_1000_4000_50000.rds")
##################################################################
faux_bergm_list_truth_post.mean <-  apply(faux_bergm_list_truth$Theta,2,mean) 
faux_bergm_post.mean=do.call(rbind, lapply(faux_bergm_list_standard, function(x) apply(x$Theta, 2, mean) ))
faux_kbergm_post.mean=do.call(rbind, lapply(faux_kbergm_list, function(x) x$post.mean))
####################################################
############# RMSE and MAE
####################################################

bergm_true_coef <- faux_bergm_list_truth_post.mean
# RMSE 
RMSE_bergm_standard <- (sweep(faux_bergm_post.mean,2,bergm_true_coef)^2) %>% apply(2,mean) %>% sqrt()
RMSE_kbergm_standard <- (sweep(faux_kbergm_post.mean,2,bergm_true_coef)^2) %>% apply(2,mean) %>% sqrt()
RMSE_kbergm_standard
RMSE_bergm_standard
# MAE
MAE_bergm_standard <- abs(sweep(faux_bergm_post.mean,2,bergm_true_coef)) %>% apply(2,mean) 
MAE_kbergm_standard <- abs(sweep(faux_kbergm_post.mean,2,bergm_true_coef)) %>% apply(2,mean) 
MAE_bergm_standard
MAE_kbergm_standard

do.call(c,lapply(faux_kbergm_list, function(x) x$runtime)) %>% mean()
do.call(c,lapply(faux_bergm_list_standard, function(x) x$Time)) %>% mean()
########################################################
resample.IS <- function(w,m,seed=0, replacement=TRUE){
  # This function implements the improved SIR with replacement and return the index
  # it returns the index, assuming that w has equal length of the original samples
  # make sure w is normalized
  # improved SIR with replacement
  # if replacement is true, by default, usiung the Skare(2003)
  # Improved SIR with replacement
  if(replacement){
    nw <- w/sum(w)
    S <- 1-nw
    # m means the samples we want
    q <- nw/S # corrected probability
    set.seed(seed)
    sample_id <- sample(1:length(w), size = m, replace = TRUE, prob = q)
  } else{
    # Andrew Gelman, BDA, without replacement
    set.seed(seed)
    sample_id <- sample(1:length(w), size = m, replace = FALSE, prob=w)
  }
  return(sample_id) # return the vector of index
}
###################################################################################
resample_id <- resample.IS(faux_kbergm_list[[7]]$nw, m=400, replacement = TRUE)
SIR_samples <- faux_kbergm_list[[7]]$param.samp[resample_id,]
par(mfrow=c(1,3))
# density of resampled data
for(i in 1:length(faux_mod$coef)){
  # density(SIR_samples[,i]
  plot(density(faux_bergm_list_truth$Theta[,i]), main=names(faux_mod$coef)[i],xlab="",ylab="", col="black",lwd=3,cex.main=2,cex.lab=2)
  lines(density(faux_bergm_list_standard[[8]]$Theta[,i]), type="l", col="red",lwd=3)
  lines(density(SIR_samples[,i]), type="l", col="blue",lwd=3)
  # density(SIR_samples[,i])
  #temp_kde <- kde( SIR_samples[,i])
  #lines(x=temp_kde$eval.points, y=temp_kde$estimate, type="l", col="red")
  abline(v=faux_mod$coef[i],col="grey",lwd=2)
  abline(v=faux_mod_mple$coef[i],lty=2,col="grey",lwd=2)
  #legend("topleft", col = c("blue","red","black"), lty = 1,  legend = c("K-ABC","AEA","Ground Truth"), lwd=1)
  if(i==1){
    legend("topleft", legend = c("Truth", "K-ABC-AIS", "AEA"),
           col=c("black", "blue", "red"), lwd=3, lty = 1, cex = 1)
  }
}
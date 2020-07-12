########################################
### This script is written for the comparison between bergm and kbergm given karate network
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
data("karate")
karate_network <- intergraph::asNetwork(karate)
# teenage_network_list <- readRDS("E:/UCI/Research/ABC for ERGMs/teenage_network_list.rds")
###################################
######## posterior mean
# faux_mod <- readRDS("E:/UCI/Research/ABC for ERGMs/faux_mod3.rds")
# faux_mod_mple <- ergm(faux_mod$formula, estimate = "MPLE")
karate_mod <- readRDS("karate_mod1.rds")
karate_mod_mple <- ergm(karate_mod$formula, estimate = "MPLE")
##################################
karate_bergm_list_truth <- readRDS("karate_bergm_list_high_quality.rds")
karate_bergm_list_standard <- readRDS("karate_bergm_list_500_1500_10000.rds")
karate_kbergm_list <- readRDS("karate_kbergm_list_4442_8000_24000.rds")
karate_kbergm_list_IS <- readRDS("karate_kbergm_list_44_32000.rds")
#######################################
karate_bergm_list_truth_post.mean <- apply(karate_bergm_list_truth$Theta,2,mean)
karate_bergm_post.mean=do.call(rbind, lapply(karate_bergm_list_standard, function(x) apply(x$Theta, 2, mean) ))
karate_kbergm_post.mean=do.call(rbind, lapply(karate_kbergm_list, function(x) x$post.mean))
karate_kbergm_IS_post.mean=do.call(rbind, lapply(karate_kbergm_list_IS, function(x) x$post.mean))
########################################
karate_post.mean <- data.frame(rbind(karate_bergm_post.mean,
                                     karate_kbergm_post.mean), factor(c(rep("AEA",nrow(karate_bergm_post.mean)),
                                                                        rep("K-ABC", nrow(karate_kbergm_post.mean))),
                                                                      levels = c("AEA", "K-ABC"),
                                                                      labels = c("AEA", "K-ABC")) )  
colnames(karate_post.mean) <- c(colnames(karate_kbergm_post.mean), "method")
########################################
################################
#####
bergm_true_coef <- karate_bergm_list_truth_post.mean

# RMSE 
RMSE_bergm_standard <- (sweep(karate_bergm_post.mean,2,bergm_true_coef)^2) %>% apply(2,mean) %>% sqrt()
RMSE_kbergm_standard <- (sweep(karate_kbergm_post.mean,2,bergm_true_coef)^2) %>% apply(2,mean) %>% sqrt()
RMSE_kbergm_IS_standard <- (sweep(karate_kbergm_IS_post.mean,2,bergm_true_coef)^2) %>% apply(2,mean) %>% sqrt()

RMSE_kbergm_IS_standard
RMSE_kbergm_standard
RMSE_bergm_standard
# MAE for posterior mean estimate
MAE_bergm_standard <- abs(sweep(karate_bergm_post.mean,2,bergm_true_coef)) %>% apply(2,mean) 
MAE_kbergm_standard <- abs(sweep(karate_kbergm_post.mean,2,bergm_true_coef)) %>% apply(2,mean) 
MAE_kbergm_IS_standard <- abs(sweep(karate_kbergm_IS_post.mean,2,bergm_true_coef)) %>% apply(2,mean) 

MAE_bergm_standard
MAE_kbergm_standard
MAE_kbergm_IS_standard


# computational time
do.call(c,lapply(karate_kbergm_list, function(x) x$runtime)) %>% mean()
do.call(c,lapply(karate_bergm_list_standard, function(x) x$Time)) %>% mean()
##################
# sampling resampling for the recovery of marginal distribution
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
####################### Karate bergm 3, kbergm 11
resample_id <- resample.IS(karate_kbergm_list[[1]]$nw, m=300, replacement = FALSE)
SIR_samples <- karate_kbergm_list[[1]]$param.samp[resample_id,]
par(mfrow=c(1,2))
# density of resampled data
for(i in 1:length(karate_mod$coef)){
  # density(SIR_samples[,i]
  plot(density(SIR_samples[,i]), main=names(karate_mod$coef)[i],xlab="",ylab="", col="blue",lwd=2,cex.main=2,cex.lab=2)
  lines(density(karate_bergm_list_standard[[19]]$Theta[,i]), type="l", col="red",lwd=3)
  lines(density(karate_bergm_list_truth$Theta[,i]), type="l", col="black",lwd=3)
  #temp_kde <- kde( SIR_samples[,i])
  #lines(x=temp_kde$eval.points, y=temp_kde$estimate, type="l", col="red")
  abline(v=karate_mod$coef[i],col="grey",lwd=3)
  abline(v=karate_mod_mple$coef[i],lty=2,col="grey",lwd=2)
  # 
  if(i==2){
     legend("topright", legend = c("Truth", "K-ABC-AIS", "AEA"),
            col=c("black", "blue", "red"), lwd=3, lty = 1, cex = 0.9)
   }
  
}
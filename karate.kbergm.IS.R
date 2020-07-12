###############
wants <- c("stringr", "parallel", "ergm",
           "mvtnorm", "network", "NetData", "sand", 
           "intergraph", "doParallel", "foreach","Bergm") 
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, library, character.only=TRUE)
#####################################################################################################
data("karate")
karate_network <- intergraph::asNetwork(karate)
####################################
### check the replicability of the optimal setting for karate network 
####################################
target_mod <- readRDS("../karate_mod1.rds")
source("./KABC_submission_codes/krergm_1.R")
######################################
dist_type = "mahalanobis"
#readline(prompt = "Enter the distance type, euclidean, mahalanobis or scaled_euclidean : ")
mc_cores = 30
prior_sd = 10
est_type = "N-D"
imp_rounds <- 1
samp_size <- 32000
samp_df <- 4
samp_scale_fact <- 4
# for(i in 1:imp_rounds){
#  samp_size[i] <- readline(prompt = paste("Enter the sample size for the ", i, "-th importance sampling step : ", sep="") ) %>% as.numeric() 
#  samp_df[i] <- readline(prompt = paste("Enter the df for the ", i, "-th importance sampling step : ", sep="") ) %>% as.numeric()
#  samp_scale_fact[i] <- readline(prompt = paste("Enter the scale factor for the ", i, "-th importance sampling step : ", sep="") ) %>% as.numeric()
# }
#samp_size <- readline(prompt = "Enter the sample size for importance sampling step : ") %>% as.numeric()
MPLE_type <- "glm"
bw_type <- "single"
mcmc.burnin <- 10000
#bw_type <- "single"
# if(MPLE_type == "NULL"){ MPLE_type = NULL }
# if(bw_type == "elementwise"){ k=readline(prompt = "Enter the k for KNN CV : ") %>% as.numeric() }
# if(bw_type == "percentage"){ percentage=readline(prompt = "Enter the percentage for bandwidth: ") %>% as.numeric() }
transf <- "none"
kernel_fun = "dnorm"
kernel_type = "gaussian"
bw_fun = "bw.nrd0"
splinefun_method = "monoH.FC"
M <- 20
percentage = 0.01
##################################################
kbergm_list <- vector(mode="list", length = M)
for(i in 1:M){
 start_time <- Sys.time()
 kbergm_list[[i]] <- kbergm(form=target_mod$formula,prior.mean=0, prior.sd=prior_sd, samp.size=samp_size, 
                            samp.df = samp_df, samp.scale.fact = samp_scale_fact , est.type = est_type, dist.type = dist_type,
                            transf=transf, bw.fun = get(bw_fun),kernel.type=kernel_type, kernel.fun = get(kernel_fun),
                            bw.type = bw_type, percentage = percentage,
                            mc.cores = mc_cores, mcmc.burnin=mcmc.burnin,splinefun.method = splinefun_method,resample = TRUE,
                            MPLE.type = MPLE_type,threshold=NULL, imp.rounds = imp_rounds,k=k,seed.base=201812+i)
 end_time <- Sys.time()
 kbergm_list[[i]]$time = end_time - start_time
 print( paste("This run takes ", end_time-start_time, " to fit", sep="") )
 print( paste("The estimated posterior mean is", kbergm_list[[i]]$post.mean) )
 cat("\n")
}
############################################################################################
saveRDS(kbergm_list,
        "karate_kbergm_list_44_32000.rds")
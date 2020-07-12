###############
wants <- c("stringr", "parallel", "ergm",
           "mvtnorm", "network", "NetData", "sand", 
           "intergraph", "doParallel", "foreach","Bergm") 
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, library, character.only=TRUE)
data("faux.mesa.high")
#load data
####################################
### check the replicability of the optimal setting for karate network 
####################################
target_mod <- readRDS("faux_mod3.rds")
source("krergm_1.R")
dist_type = "mahalanobis"
mc_cores = 30
####################################
prior_sd_edges = 2.236
prior_sd_nodematch = 2.236
prior_sd_gwesp = 2.236
####################################
prior_edges = -2
prior_nodematch = 0.5
prior_gwesp = 0.5
transf_id = 1:3

est_type = "N-D"
imp_rounds <- 2
samp_size <- c(24000,96000)
samp_df <- c(4,4)
samp_scale_fact <- c(4,2)
# #####################################
# for(i in 1:imp_rounds){
#  samp_size[i] <- readline(prompt = paste("Enter the sample size for the ", i, "-th importance sampling step : ", sep="") ) %>% as.numeric() 
#  samp_df[i] <- readline(prompt = paste("Enter the df for the ", i, "-th importance sampling step : ", sep="") ) %>% as.numeric()
#  samp_scale_fact[i] <- readline(prompt = paste("Enter the scale factor for the ", i, "-th importance sampling step : ", sep="") ) %>% as.numeric()
# }
#samp_size <- readline(prompt = "Enter the sample size for importance sampling step : ") %>% as.numeric()
MPLE_type <- "glm"
bw_type <- "single"
mcmc.burnin <- 5*10^4
########################################
#bw_type <- "single"
# if(MPLE_type == "NULL"){ MPLE_type = NULL }
# if(bw_type == "elementwise"){ k=readline(prompt = "Enter the k for KNN CV : ") %>% as.numeric() }
# if(bw_type == "percentage"){ percentage=readline(prompt = "Enter the percentage for bandwidth: ") %>% as.numeric() }
transf <- "sqrt1"
kernel_fun = "dnorm"
kernel_type = "gaussian"
#kernel.type=c("gaussian","epanechnikov")
bw_fun = "bw.nrd0"
splinefun_method = "monoH.FC"
M <- 20
percentage = 0.01
############################################################################################
kbergm_list <- vector(mode="list", length = M)
for(i in 1:M){
 start_time <- Sys.time()
 kbergm_list[[i]] <- kbergm(target_mod$formula,prior.mean= c(prior_edges, prior_nodematch, prior_gwesp ), prior.sd=c(prior_sd_edges,prior_sd_nodematch,prior_sd_gwesp), samp.size=samp_size, 
                            samp.df = samp_df, samp.scale.fact = samp_scale_fact , est.type = est_type, dist.type = dist_type, bw.fun = get(bw_fun),
                            kernel.type=kernel_type, kernel.fun = get(kernel_fun),bw.type = bw_type, 
                             transf=transf,transf.id=transf_id,percentage = percentage,mc.cores = mc_cores, 
                            mcmc.burnin=mcmc.burnin, splinefun.method = splinefun_method,resample = TRUE,MPLE.type = MPLE_type,threshold=NULL, 
                            imp.rounds = imp_rounds,k=k, cov.estimate = TRUE,seed.base=201812+i)
 end_time <- Sys.time()
 kbergm_list[[i]]$time = end_time - start_time
 print( paste("This run takes ", end_time-start_time, " to fit", sep="") )
 print( paste("The estimated posterior mean is", kbergm_list[[i]]$post.mean) )
 cat("\n")
}
###########################################################################################
saveRDS(kbergm_list,
        "faux_kbergm_list_4442_24000_96000_sqrt1.rds")
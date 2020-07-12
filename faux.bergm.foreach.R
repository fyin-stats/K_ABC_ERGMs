wants <- c("stringr", "parallel", "ergm",
           "mvtnorm", "network", "NetData", "sand", 
           "intergraph", "doParallel", "foreach","Bergm") 
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, library, character.only=TRUE)
data("faux.mesa.high")
###################################
###################################
target_mod <- readRDS("faux_mod3.rds")
#start_time <- Sys.time()
burn_in = 1000
main_iters = 4000
aux_iters = 5*10^4
prior_sd_edges = 2.236
prior_sd_nodematch = 2.236
prior_sd_gwesp = 2.236

prior_edges = -2
prior_nodematch = 0.5
prior_gwesp = 0.5
 
M <- 20
num_cores <- detectCores()
####################################
registerDoParallel(cores=num_cores)
faux_bergm_list <- foreach(i = 1:M)%dopar%{
 start_time <- Sys.time()
 set.seed(201812+i)
 faux_bergm <- bergm(target_mod$formula, burn.in = burn_in, main.iters = main_iters, aux.iters = aux_iters, 
                        prior.mean = c(prior_edges, prior_nodematch, prior_gwesp ),  
                        prior.sigma = diag( c(prior_sd_edges,prior_sd_nodematch,prior_sd_gwesp)^2 ), nchains = length(target_mod$coef)*2 , gamma = 0.5) 
 #karate_bergm_default <- bergm(karate_mod1$formula)
 end_time <- Sys.time()
 #faux_bergm$time <- end_time-start_time
 print( paste("This run takes ", difftime(end_time,start_time), " to fit", sep="") )
 faux_bergm
}
#####################################
saveRDS(faux_bergm_list,
        "faux_bergm_list_standard_1000_4000_50000.rds")
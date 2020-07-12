wants <- c("stringr", "parallel", "ergm",
           "mvtnorm", "network", "NetData", "sand", 
           "intergraph", "doParallel", "foreach","Bergm") 
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, library, character.only=TRUE)
####################################################
data("faux.mesa.high")
###################################
target_mod <- readRDS("faux_mod3.rds")
burn_in = 4000
main_iters = 16000
aux_iters = 5*10^5
prior_sd = 2.236
#prior_edges = -2
####################################
set.seed(123456)
start_time <- Sys.time()
faux_bergm <- bergm(target_mod$formula, burn.in = burn_in, main.iters = main_iters, aux.iters = aux_iters, 
                        prior.mean = c(-2,0.5,0.5),  
                        prior.sigma = diag(prior_sd^2, length(target_mod$coef)), nchains = length(target_mod$coef)*2, gamma = 0.5) 
end_time <- Sys.time()
print( paste("This run takes ", end_time-start_time, " to fit", sep="") )
saveRDS(faux_bergm,
        "faux_bergm_list_20505_sqrt5_super_high_quality_long.rds")
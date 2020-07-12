wants <- c("stringr", "parallel", "ergm",
           "mvtnorm", "network", "NetData", "sand", 
           "intergraph", "doParallel", "foreach","Bergm") 
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, library, character.only=TRUE)
###################################
data("karate")
###################################
#### consider a small network, say lazega, or karate
karate_network <- intergraph::asNetwork(karate)
###################################
karate_mod1 <- readRDS("karate_mod1.rds")
###################################
burn_in = 2500
main_iters = 12500
aux_iters = 10^5
prior_sd = 10
####################################
start_time <- Sys.time()
set.seed(123456)
karate_bergm <- bergm(karate_mod1$formula, burn.in = aux_iters, main.iters = main_iters,
                      aux.iters = aux_iters, prior.mean = rep(0,length(karate_mod1$coef)),
                      prior.sigma = diag(prior_sd^2, length(karate_mod1$coef)), nchains = length(target_mod$coef)*2, gamma=0.5)
end_time <- Sys.time()
####################################
saveRDS(karate_bergm,
        "karate_bergm_list_high_quality.rds")
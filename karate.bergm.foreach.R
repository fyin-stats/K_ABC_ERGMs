wants <- c("stringr", "parallel", "ergm",
           "mvtnorm", "network", "NetData", "sand", 
           "intergraph", "doParallel", "foreach","Bergm") 
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, library, character.only=TRUE)
###################################
data("karate")
###################################
karate_network <- intergraph::asNetwork(karate)
# karate_mod1 <- readRDS("/home/fyin/ABC_for_ERGMs/karate_mod1.rds")
target_mod <- readRDS("karate_mod1.rds")
#start_time <- Sys.time()
####################################
###################################
burn_in = 500
main_iters = 1500
aux_iters = 10^4
prior_sd = 10
prior_edges = 0 
####################################
M = 20
num_cores <- detectCores()
####################################
registerDoParallel(cores=num_cores)
####################################
karate_bergm_list <- foreach(i = 1:M)%dopar%{
 start_time <- Sys.time()
 set.seed(201812+i)
 karate_bergm <- bergm(target_mod$formula, burn.in = burn_in, main.iters = main_iters, aux.iters = aux_iters, 
                        prior.mean = c(prior_edges, rep(0,length(target_mod$coef)-1) ),  
                        prior.sigma = diag(prior_sd^2, length(target_mod$coef) ), gamma = 0.5) 
 end_time <- Sys.time()
 karate_bergm$time <- end_time-start_time
 print( paste("This run takes ", end_time-start_time, " to fit", sep="") )
 karate_bergm
}
## calculate posterior mean
saveRDS(karate_bergm_list,
        "karate_bergm_list_500_1500_10000.rds")
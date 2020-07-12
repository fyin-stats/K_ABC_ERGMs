wants <- c( "parallel",  "dplyr",  "ergm", 
            "mvtnorm", "network", "NetData", "sand", "intergraph", "doParallel", "foreach", "Bergm") 
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, library, character.only=TRUE)

data("karate")
data("faux.mesa.high")
###################################
#### consider a small network, say lazega, or karate
karate_network <- intergraph::asNetwork(karate)
####################################
### check the replicability of the optimal setting for karate network 
####################################
target_mod <- readRDS("/faux_mod3.rds")
source("krergm_1.R")
dist_type = "mahalanobis"
mc_cores = rev(c(1,2,5,8,10,15,20,25,30,35))
prior_sd = 2.236
est_type = "N-D"
imp_rounds <- 2
samp_size <- c(24000,96000)
samp_df <- c(4,4)
samp_scale_fact <- c(4,2)
MPLE_type <- "glm"
bw_type <- "single"
mcmc.burnin <- 50000

if(MPLE_type == "NULL"){ MPLE_type = NULL }
if(bw_type == "elementwise"){ k=readline(prompt = "Enter the k for KNN CV : ") %>% as.numeric() }
if(bw_type == "percentage"){ percentage=readline(prompt = "Enter the percentage for bandwidth: ") %>% as.numeric() }
transf <- "sqrt1"
kernel_fun = "dnorm"
kernel_type = "gaussian"
bw_fun = "bw.nrd0"
splinefun_method = "monoH.FC"
M <- 1
percentage = 0.01

# computational time
mpletime <- rep(NA, length(mc_cores))
samptime <- rep(NA, length(mc_cores))
runtime <- rep(NA, length(mc_cores))
qtime <- rep(NA, length(mc_cores))

kbergm_list <- vector(mode="list", length = length(mc_cores))
for(j in 1:length(mc_cores)){
 kbergm_list[[j]] <- vector(mode="list", length = M)
 for(i in 1:M){
  start_time <- Sys.time()
  kbergm_list[[j]][[i]] <- kbergm(target_mod$formula,prior.mean= c(-2,0.5,0.5), prior.sd=c(2.236,2.236,2.236), samp.size=samp_size, samp.df = samp_df, samp.scale.fact = samp_scale_fact , est.type = est_type, dist.type = dist_type,transf=transf, bw.fun = get(bw_fun),kernel.type=kernel_type, kernel.fun = get(kernel_fun) ,bw.type = bw_type, percentage = percentage,
                            mc.cores = mc_cores[j], mcmc.burnin=mcmc.burnin,splinefun.method = splinefun_method,resample = TRUE,MPLE.type = MPLE_type,threshold=NULL, imp.rounds = imp_rounds,k=k)
  end_time <- Sys.time()
  #kbergm_list[[i]]$time = end_time - start_time
  print( paste("This run takes ", end_time-start_time, " to fit", sep="") )
  print( paste("The estimated posterior mean is", kbergm_list[[j]][[i]]$post.mean) )
  cat("\n")
 }
 mpletime[j] <- mean(do.call(c,lapply(kbergm_list[[j]], function(x) x$fittime)))
 samptime[j] <- mean(do.call(c,lapply(kbergm_list[[j]], function(x) x$samptime)))
 runtime[j] <- mean(do.call(c,lapply(kbergm_list[[j]], function(x) x$runtime)))
 qtime[j] <- mean(do.call(c,lapply(kbergm_list[[j]], function(x) x$qtime)))
}

faux_time <- cbind(mpletime, samptime, runtime, qtime, mc_cores)
colnames(faux_time) <- c("MPLE", "Sample time", "Total run time", "Quantile finding time", "Number of cores")
saveRDS(faux_time, "faux_time.rds")
#
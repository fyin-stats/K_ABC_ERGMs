#Another take on kernel regression for approximate Bayesian ERGM inference
#(may be redundant with my other stuff, but oh well).  This file has some
#pretty simple/straightforward tools that are just focused on kernel Bayes
#inference by importance sampling from an MPLE-centered distribution.
#
#CTB, 12/1/18
wants <- c("ergm", "parallel", "data.table", "dplyr", "doParallel",
           "ggplot2", "foreach", "quantreg",
           "splines2", "splines", "logspline", "scam", "abind", "Bergm", "mvtnorm") 
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, library, character.only=TRUE)
#if(any(!has)) try(install.packages(wants[!has]))
#lapply(wants, function(x) try(library(x)), character.only=TRUE)
####################################################
rmvt<-function(n,mu=0,cov=diag(1),df=1){
  require(MASS)
  draws<-matrix(nrow=n,ncol=length(mu))
  for(i in 1:n)
    draws[i,]<-mvrnorm(n=1,mu=rep(0,length(mu)),Sigma=cov)/sqrt(rchisq(1,df=df)/df)+mu
  if(n>1)
    draws
  else
    draws[1,]
}

dmvt<-function(x,mu=0,cov=diag(1),df=1,log=FALSE){
  if(is.null(dim(x)))
    x<-matrix(x,nrow=1)
  x<-sweep(x,2,mu,"-")
  p<-length(mu)
  n<-NROW(x)
  lp<-rep(lgamma((df+p)/2)-lgamma(df/2)-p/2*log(df)-p/2*log(pi)-0.5*determinant(cov)$modulus, n)
  icov<-solve(cov)
  for(i in 1:n){
    lp[i]<-lp[i]-(df+p)/2*log1p(t(x[i,])%*%icov%*%x[i,]/df)
  }
  if(log)
    lp
  else
    exp(lp)
}

dmvnorm<-function(x,mu=0,cov=diag(1),df=1,log=FALSE){
  if(is.null(dim(x)))
    x<-matrix(x,nrow=1)
  x<-sweep(x,2,mu,"-")
  p<-length(mu)
  n<-NROW(x)
  lp<-rep(-p/2*log(2*pi)-0.5*determinant(cov)$modulus, n)
  icov<-solve(cov)
  for(i in 1:n){
    lp[i]<-lp[i]-0.5*t(x[i,])%*%icov%*%x[i,]
  }
  if(log)
    lp
  else
    exp(lp)
}
##################################################################################################
################# kbergm currently supports the estimation based on Nadaraya-Watson estimator, GAM,... (potentially local linear, local polynomial...)
################# GAM not satisfying
##################################################################################################
kbergm<-function(form,prior.mean=0,prior.sd=10,samp.size=1e3,samp.df=4,
                 samp.scale.fact=4,est.type="N-D",transf=c("none", "log1p", "sqrt1"),
                 transf.id=NULL,dist.type=c("euclidean","mahalanobis","scaled_euclidean"),
                 bw.fun=bw.nrd0,kernel.fun = dnorm,kernel.type=c("gaussian","epanechnikov"),
                 bw.type=c("single","elementwise","percentage"),
                 percentage=0.01,mc.cores=1,splinefun.method = "monoH.FC",
                 resample=FALSE,resample.size=500,resample.replace=FALSE,mcmc.burnin=NULL, 
                 mcmc.interval=1,MPLE.type=c("glm","penalized"),threshold=0.1,imp.rounds=1,k=50,
                 cv.method="optim",cov.estimate=FALSE,seed.base=123,clipping=NULL,...){
  require(parallel)
  require(sna)
  require(ergm)
  require(doParallel)
  require(foreach)
  RNGkind("L'Ecuyer-CMRG")
  #Fit the MPLE
  runtime1 <- Sys.time()
  cat("Performing initial MPLE fit...\n")
  if(!is.null(MPLE.type)){
   ifit<-ergm(form,estimate="MPLE", control=control.ergm(MPLE.type = MPLE.type)) }
  else{
   ifit<-ergm(form,control=control.ergm(main.method = "MCMLE"))
  }
  runtime2 <- Sys.time()
  samp.size.1 <- samp.size[1] # sample size for first round
  samp.scale.fact.1 <- samp.scale.fact[1]
  samp.df.1 <- samp.df[1]
  ostats<-ifit$target.stats  #Get target stats
  #Draw the importance sample
  cat("Drawing the importance sample\n")
  sscale<-ifit$covar*samp.scale.fact.1
  smu<-ifit$coef
  #parsamp<-rmvt(samp.size.1,mu=smu,cov=sscale,df=samp.df.1)
  set.seed(seed.base)
  parsamp<-mvtnorm::rmvt(samp.size.1,sigma = sscale,df=samp.df.1,delta=smu)  
  #Obtain the importance ratio
  cat("Getting importance ratios\n")
  p<-length(smu)
  samp.mu <- matrix(0, nrow=imp.rounds, ncol=p)
  samp.scale.array <- array(0, dim=c(p,p,imp.rounds)) 
  samp.mu[1,] <- smu
  samp.scale.array[,,1] <- sscale
  slprob<- mvtnorm::dmvt(parsamp, sigma = sscale, df = samp.df.1, delta = smu, log=TRUE)
  #slprob<-mvtnorm::dmvt(parsamp,mu=smu,cov=sscale,df=samp.df.1,log=TRUE)
  if( length(prior.mean) == 1){
    prior.mean<-rep(prior.mean,length=p) }
  if( length(prior.sd) == 1 ){
    prior.sd<-rep(prior.sd,length=p) }
  if(p>1)
    prior.cov<-diag(prior.sd^2)
  else
    prior.cov<-diag(1)*prior.sd^2
  plprob<-mvtnorm::dmvnorm(parsamp,mean=prior.mean,sigma=prior.cov,log=TRUE)
  #plprob<-dmvnorm(parsamp,mu=prior.mean,cov=prior.sdm,log=TRUE)
  iw<-plprob-slprob
  #Draw stats
  #time1 <- Sys.time()
  cat("Drawing the statistics at the sampled parameters\n")
  
  y<-ergm.getnetwork(form)
  mod <- ergm_model(form, y) # internal representation of an ergm network model
  sy <- summary_formula(form)
  if( any(is.na(as.matrix.network(y)))){
    print("Network has missing data.")
  }
  # ergm.Cprepare; collates the information in the given object into a form suitable for
  # being passed to the C routines
  Clist <- ergm.Cprepare(y, mod)
  control <- control.ergm(MCMC.burnin = ifelse(is.null(mcmc.burnin), 10000, as.numeric(mcmc.burnin)), 
                          MCMC.interval = ifelse(is.null(mcmc.interval), 1, as.numeric(mcmc.interval)), MCMC.samplesize = 1)
  proposal <- ergm_proposal(object = ~., constraint = ~., arguments = control$MCMC.prop.args, nw=y)
  registerDoParallel(cores = mc.cores)
  #set.seed(seed.base*377)
  stats <- t(matrix(unlist(mclapply(1:samp.size.1,function(z){
                        ergm_MCMC_slave(Clist = Clist, proposal = proposal, eta=parsamp[z,], control = control,
                  verbose = FALSE)$s},mc.cores=mc.cores, mc.preschedule = TRUE,mc.set.seed = TRUE)),ncol=samp.size))
  stats <- sweep(stats,2,sy,"+")
  samp.stats <- stats
  target.stats<-ostats
  # which element do you want to transform? by default, if transf.id = NULL, we transform all dimensions...
  if(is.null(transf.id)){
    transf.id <- 1:p
  }
  
  if(match.arg(transf) == "log1p"){
    stats <- log1p(stats)
    target.stats <- log1p(ostats)
  } else if(match.arg(transf)=="sqrt1"){
    stats <- sqrt(stats+1)
    target.stats <- sqrt(ostats+1)
  }
  
  cat("Obtaining distances to target, and kernel weights\n")
  if(match.arg(dist.type)=="euclidean"){
    d<-rowSums(sweep(stats,2,target.stats,"-")^2)^0.5
  }else if(match.arg(dist.type) == "mahalanobis"){
    cov_stats <- cov(stats)
    d<-sqrt(mahalanobis(stats,target.stats,cov_stats))
  } else{
     #scaled distance as suggested by Prangle (2017), use MAD
     # cov_stats <- cov(stats)
     median_stats <- apply(stats, 2, median)
     MAD <- apply(abs(sweep(stats, 2, median_stats)),2,median)
     S <- diag(MAD^2)
     d <- sqrt(mahalanobis(stats,target.stats,S))
   }
  if(!is.null(threshold)){
    # only want to keep the samples that are sufficiently close to the observed stats, threshold corresponds to quantile
    keep_index <- (d < quantile(d,threshold)) # d is the mahalanobis distance between observed stats and simulated stats
    iw<-iw[keep_index]
    parsamp<-parsamp[keep_index,]    
    stats<-stats[keep_index,]
    d<-d[keep_index]
  }
  # heuristic way to determine the bandwidth
  Eqhat <- rep(NA, length = p)
  Eq2hat <- rep(NA, length = p) #just in case needed
  if(match.arg(bw.type) %in% c("single", "percentage")){
    #using univariate kernel density estimate as a heuristic
    if(match.arg(bw.type)=="percentage"){
     # by default 1%
     bw <- sort(d,decreasing=FALSE)[ floor(samp.size.1 * percentage) ] 
     kw <- kernel.fun(d/bw, log=FALSE)
     kw <- kw/sum(kw)
     w <- exp(iw) * kw
     } else{
     bw<-density(d, bw="nrd0",kernel = kernel.type)$bw 
     kw<-kernel.fun(d/bw,log=TRUE)
     kw<-kw-logSum(kw)
     # by default, bw.nrd0 
     w<-exp(iw+kw)
    }
    # clipping method
    if(!is.null(clipping)){
      clipping_threshold <- sort(w,decreasing = TRUE)[floor(samp.size.1*clipping)]
      w <- sapply(w, function(x) min(x,clipping_threshold))
      w[w<clipping_threshold] = 0 # clip those that do not qualify
    }
     #Calculate the kernel density estimate of
    #E_s q p(q)/s(q) | Y=y
    cat("Calculating kernel density estimates\n")
    sw<-sum(w) #sum of weights
    nw<-w/sw
    Eqhat<-colSums(sweep(parsamp,1,w,"*"))/sw
    Eq2hat<-colSums(sweep(parsamp^2,1,w,"*"))/sw   #locally constant estimate of the function of posterior samples
    bw<-rep(bw, p) # make it a vector of repeated values
    nw<-matrix(rep(nw,p),nrow=samp.size[1],ncol=p)
   } else if(bw.type=="elementwise"){
    bw<-rep(NA,p)
    kw<-matrix(0,nrow=samp.size[1],ncol=p) #may consider more efficient way to store kw and nw
    nw<-matrix(0,nrow=samp.size[1],ncol=p)
   for(i in 1:p){
    # CV.dist.bw : a function that computes the bandwidth to scale the distance for Nadaraya-watson kernel estimator
    # kNN.CV.bw <- function(k=50,d, params, hmin, hmax, nbh)
    bw0 <- density(d, kernel = kernel.type)$bw 
    #browser()
    bw[i]<-kNN.CV.bw(k=k,d,kernel.fun=kernel.fun,params=parsamp[,i],method=cv.method,hmin=min(d)/2,hmax=sort(d,decreasing=FALSE)[ floor(samp.size.1 * percentage) ],nbh=200,num_cores = mc.cores)$hopt
    kw[,i]<-kernel.fun(d/bw[i],log=TRUE)
    kw[,i]<-kw[,i]-logSum(kw[,i])
    #Calculate the kernel density estimate of
    #E_s q p(q)/s(q) | Y=y
    cat("Calculating kernel density estimates\n")
    w<-exp(iw+kw[,i])
    sw<-sum(w) #sum of weights
    nw[,i]<-w/sw
    Eqhat[i]<-sum(parsamp[,i] * nw[,i])
    Eq2hat[i]<-sum(parsamp[,i]^2 * nw[,i])
    }
  }
  #Varqhat<-Eq2hat-Eqhat^2
  # if multiple rounds of importance sampling
  if(imp.rounds>1){
   cat("Multiple rounds of importance sampling needed : ", imp.rounds, "\n")
   for(i in 2:imp.rounds){
     if( length(samp.size) == 1){
       temp.samp.size <- samp.size
      }
     else if(length(samp.size) == imp.rounds){
       temp.samp.size <- samp.size[i]
      }
     if(length(samp.scale.fact) == 1){
       temp.samp.scale.fact <- samp.scale.fact
     }
     else if( length(samp.scale.fact) == imp.rounds ){
       temp.samp.scale.fact <- samp.scale.fact[i]
     }
     if(length(samp.df)==1){
        temp.samp.df = samp.df
     }
     else if(length(samp.df) == imp.rounds){
        temp.samp.df = samp.df[i]
     }
    cat("Now starting round", i, "\n")
    new_mu <- Eqhat
    new_sscale <- (t(sweep(sweep(parsamp,2,Eqhat),1,w,"*")) %*% (sweep(parsamp,2,Eqhat)) / sw) * temp.samp.scale.fact
    samp.mu[i,] <- new_mu
    samp.scale.array[,,i] <- new_sscale
    print(new_mu)
    print(new_sscale)
    #set.seed(seed.base*i)
    parsamp<-mvtnorm::rmvt(temp.samp.size,sigma = new_sscale,df=temp.samp.df,delta=new_mu)  
    cat("Getting importance ratios\n")
    p<-length(smu)
    slprob<- mvtnorm::dmvt(parsamp, sigma = new_sscale, df = temp.samp.df, delta = new_mu, log=TRUE)
    #slprob<-dmvt(parsamp,mu=new_mu,cov=new_sscale,df=temp.samp.df,log=TRUE)
    if( length(prior.mean) == 1){
     prior.mean<-rep(prior.mean,length=p) }
    if( length(prior.sd) == 1 ){
     prior.sd<-rep(prior.sd,length=p) }
    if(p>1)
     prior.cov<-diag(prior.sd^2)
    else
     prior.cov<-diag(1)*prior.sd^2
    #plprob<-dmvnorm(parsamp,mu=prior.mean,cov=prior.sdm,log=TRUE)
    plprob<-mvtnorm::dmvnorm(parsamp,mean=prior.mean,sigma=prior.cov,log=TRUE)
    iw<-plprob-slprob
    # truncated importance sampling
    #Draw stats
    cat("Drawing the statistics at the sampled parameters\n") 
    set.seed(seed.base*i+i)
    stats <- t(matrix(unlist(mclapply(1:temp.samp.size,function(z){
                      ergm_MCMC_slave(Clist = Clist, proposal = proposal, eta=parsamp[z,], control = control,
                      verbose = FALSE)$s},mc.cores=mc.cores, mc.preschedule = TRUE,mc.set.seed = TRUE)),ncol=temp.samp.size))

    stats <- sweep(stats,2,sy,"+")
    samp.stats <- stats # this corresponds to the actual sampled stats
    # stats is the working stats used for reweighting
    # transformation
    if(match.arg(transf) == "log1p"){
      #transform the given dimensions only
      stats <- log1p(stats)
      #target.stats <- log1p(ostats)
    } else if(match.arg(transf)=="sqrt1"){
      #transform the given dimensions only
      stats <- sqrt(stats+1)
      #target.stats <- sqrt(ostats+1)
    }
    #Get distances from the target, and kernel weights
    cat("Obtaining distances to target, and kernel weights\n")
    if(match.arg(dist.type)=="euclidean"){
      d<-rowSums(sweep(stats,2,target.stats,"-")^2)^0.5
    }else if(match.arg(dist.type) == "mahalanobis"){
      cov_stats <- cov(stats)
      d<-sqrt(mahalanobis(stats,target.stats,cov_stats))
    } else{
      #scaled distance as suggested by Prangle (2017)
      median_stats <- apply(stats, 2, median)
      MAD <- apply(abs(sweep(stats, 2, median_stats)),2,median)
      S <- diag(MAD^2)
      d <- sqrt(mahalanobis(stats,target.stats,S))
    }
    Eqhat <- rep(NA, length = p)
    Eq2hat <- rep(NA, length = p) #just in case needed
    if(match.arg(bw.type) %in% c("single", "percentage")){
      #using univariate kernel density estimate as a heuristic
      if(match.arg(bw.type)=="percentage"){
        # by default 1%
        bw <- sort(d,decreasing=FALSE)[ floor(samp.size.1 * percentage) ] 
        kw <- kernel.fun(d/bw, log=FALSE)
        kw <- kw/sum(kw)
        w <- exp(iw) * kw
       }
       else{ bw<-density(d, bw="nrd0",kernel = kernel.type)$bw # by default, bw.nrd0
             kw<-kernel.fun(d/bw,log=TRUE)
             kw<-kw-logSum(kw)
             # by default, bw.nrd0 
             w<-exp(iw+kw)
       }     
      # kernel weights and importance weights
      # use clipping method
      if(!is.null(clipping)){
        clipping_threshold <- sort(w,decreasing = TRUE)[floor(temp.samp.size*clipping)]
        w <- sapply(w, function(x) min(x,clipping_threshold))
        w[w<clipping_threshold] = 0 # clip those that do not qualify
      }
      #Calculate the kernel density estimate of
      #E_s q p(q)/s(q) | Y=y
      cat("Calculating kernel density estimates\n")
      #w<-exp(iw+kw)
      sw<-sum(w) #sum of weights
      nw<-w/sw
      Eqhat<-colSums(sweep(parsamp,1,w,"*"))/sw
      Eq2hat<-colSums(sweep(parsamp^2,1,w,"*"))/sw   #locally constant estimate of the function of posterior samples
      bw<-rep(bw, p) # make it a vector of repeated values
      nw<-matrix(rep(nw,p),nrow=temp.samp.size,ncol=p)
    } else if(bw.type=="elementwise"){
      bw<-rep(NA,p)
      kw<-matrix(0,nrow=temp.samp.size,ncol=p) #may consider more efficient way to store kw and nw
      nw<-matrix(0,nrow=temp.samp.size,ncol=p)
      for(i in 1:p){
        # CV.dist.bw : a function that computes the bandwidth to scale the distance for Nadaraya-watson kernel estimator
        # kNN.CV.bw <- function(k=50,d, params, hmin, hmax, nbh)
        bw0 <- density(d, kernel = kernel.type)$bw
        bw[i]<-kNN.CV.bw(k=k,d,kernel.fun=kernel.fun,params=parsamp[,i],method=cv.method,hmin=min(d)/2,hmax=sort(d,decreasing=FALSE)[ floor(temp.samp.size * percentage) ],
                         nbh=200,num_cores = mc.cores)$hopt
        kw[,i]<-kernel.fun(d/bw[i],log=TRUE)
        kw[,i]<-kw[,i]-logSum(kw[,i])
        #Calculate the kernel density estimate of
        #E_s q p(q)/s(q) | Y=y
        cat("Calculating kernel density estimates\n")
        w<-exp(iw+kw[,i])
        sw<-sum(w) #sum of weights
        nw[,i]<-w/sw
        Eqhat[i]<-sum(parsamp[,i] * nw[,i])
        Eq2hat[i]<-sum(parsamp[,i]^2 * nw[,i])
      }
    }
   }
  }
  Varqhat<-Eq2hat-Eqhat^2
  # estimate the covariance matrix?
  runtime3 <- Sys.time()
  if(cov.estimate){
    Covqhat <- (t(sweep(sweep(parsamp,2,Eqhat),1,w,"*")) %*% (sweep(parsamp,2,Eqhat)) / sw)
  } else{
    Covqhat <- NA
  }
  #Often, we want to know how well-determined the sign of an effect is (Bayesian
  #equivalent of a 2-tailed test).  We can just estimate that directly...
  # pge0<-colSums(sweep(parsamp>=0,1,w,"*"))/sw
  #START Obsoleted bound-based intervals-----
  #Get some estimated posterior intervals, using conservative estimates; we can
  #safely assume unimodality, and hence the Vysochanskii-Petunin inequality:
  #  Pr(|X-E(X)|>k*sigma) <= 4/9 k^-2
  #which implies the following k values:
  #  99%: 6.666667
  #  95%: 2.981424
  #  50%: 0.942809
  #(these are all based on sqrt((4/z)/9), where z is the tail weight.
  #The V-P intervals aren't all that great, but they are quite safe...
  #pi99<-rbind(Eqhat-6.666667*sqrt(Varqhat),Eqhat+6.666667*sqrt(Varqhat))
  #pi95<-rbind(Eqhat-2.981424*sqrt(Varqhat),Eqhat+2.981424*sqrt(Varqhat))
  #pi50<-rbind(Eqhat-0.942809*sqrt(Varqhat),Eqhat+0.942809*sqrt(Varqhat))
  ##############################
  # V-P intervals
  #pi99_vp<-rbind(Eqhat-6.666667*sqrt(Varqhat),Eqhat+6.666667*sqrt(Varqhat))
  #pi95_vp<-rbind(Eqhat-2.981424*sqrt(Varqhat),Eqhat+2.981424*sqrt(Varqhat))
  #pi50_vp<-rbind(Eqhat-0.942809*sqrt(Varqhat),Eqhat+0.942809*sqrt(Varqhat))
  #END Obsoleted bound-based intervals-----
  #Estimate posterior marginals using an approximate ECDF technique; we can
  #directly estimate the CDF at each draw by recognizing that each draw has
  #weight proportional to exp(iw+kw).  Here, we we use spline rather than a
  #true ECDF, in order to smooth things out.  We can then take derivatives
  #in order to estimate the marginal density.
  #timeQ1 <- Sys.time()
  findQ<-function(q,theta,pr,cdf){
    #Since optimize sucks, first narrow things down manually
    if(pr[1]>=q){
      lb<-theta[1]
    }else{
      lb<-theta[max(which(pr<=q))]
    }
    if(pr[length(pr)]<=q){
      ub<-theta[length(pr)]
    }else{
      ub<-theta[min(which(pr>=q))]
    }
    #Now, interpolate within this interval
    #optimize(function(z){(cdf(x=z)-q)^2},interval=c(lb,ub))$minimum
    #optimx::optimx(par=(lb+ub)/2,fn=function(z){(cdf(x=z)-q)^2},lower=lb, upper=ub,method="Nelder-Mead")$p1
    optimize(function(z){(cdf(x=z)-q)^2},interval=c(lb,ub))$minimum
  }
  postmarcdf<-vector(mode="list",length=p)
  postmarden<-vector(mode="list",length=p)
  pi99<-matrix(nrow=2,ncol=p)
  pi95<-matrix(nrow=2,ncol=p)
  pi50<-matrix(nrow=2,ncol=p)
  eff_samp_size <- rep(NA, p)
  browser()
  for(i in 1:p){
    ord<-order(parsamp[,i])
    y<-cumsum(nw[ord,i])
    x<-parsamp[ord,i]
    sf<-splinefun(x=x,y=y,method="monoH.FC")
    postmarcdf[[i]]<-list(theta=x,F=y)
    postmarden[[i]]<-list(theta=x,f=sf(x=seq(from=min(x),to=max(x),length=1e3),deriv=1))
    pi99[1,i]<-findQ(0.005,x,y,sf)
    pi99[2,i]<-findQ(0.995,x,y,sf)
    pi95[1,i]<-findQ(0.025,x,y,sf)
    pi95[2,i]<-findQ(0.975,x,y,sf)
    pi50[1,i]<-findQ(0.25,x,y,sf)
    pi50[2,i]<-findQ(0.75,x,y,sf)
    eff_samp_size[i] <- 1/(sum(nw[,i]^2))
  }
  runtime4 <- Sys.time()
   # effective sample size for importance samples
  # eff_samp_size <- 1/(sum(nw^2)) 
  if(match.arg(bw.type) %in% c("single", "percentage")){
    nw <- nw[,1] #convert nw back if univariate kernel
    #kw <- kw[,1]
    #Eqhat1<-colSums(sweep(parsamp,1,w,"*"))/sw  #importance weights attached, is this nadaraya-watson estimator?
    #Eq2hat1<-colSums(sweep(parsamp^2,1,w,"*"))/sw #importance weights attached, each draw has weight proportional to exp(iw+kw)
    #Varqhat1<-Eq2hat-Eqhat^2 #variance, a single value, an estimate of the variance b  
  }
  # convert weighted samples to unweighted samples by resampling
  if(resample & match.arg(bw.type) %in% c("single", "percentage")){
     #resample.size <- min(c(resample.size,eff_samp_size/2))
     resample_id <- resample.IS(nw,m=resample.size)
     parsamp_resample <- parsamp[resample_id,]
  }
  else if(!resample){
     parsamp_resample <- NA
  }
  else if(resample & bw.type=="elementwise"){
    parsamp_resample <- matrix(0, nrow = resample.size, ncol = p)
    for(i in 1:p){
      resample_id <- resample.IS(nw[,i], m=resample.size)
      parsamp_resample[,i] <- parsamp[resample_id,i]
    }
  }
  #Copula based tool for joint posterior sampling - now 
  #First, find pseudo-observations
  
  #psobs<-matrix(nrow=dim(parsamp)[1],ncol=p)
  #for(i in 1:p)
  #  psobs[,i]<-postmarcdf[[i]]$F[match(parsamp[,i],postmarcdf[[i]]$theta)] #pseudo-observations are between 0 and 1
  
  ######################################  
  #for the purpose of comparison
  #calculate the post mean and post sd using the old way
  bw1<-density(d, kernel = kernel.type)$bw # by default, bw.nrd0
  kw1<-kernel.fun(d/bw1,log=TRUE)
  #bw1<-bw.fun(d)
  #kw1<-dnorm(d/bw1,log=TRUE)
  kw1<-kw1-logSum(kw1)  #log-normalizing factor
  #Calculate the kernel density estimate of
  #E_s q p(q)/s(q) | Y=y
  cat("Calculating kernel density estimates\n")
  w1<-exp(iw+kw1)
  sw1<-sum(w1)
  nw1<-w1/sw1  #normalize the weight
  Eqhat1<-colSums(sweep(parsamp,1,w1,"*"))/sw1  #importance weights attached, is this nadaraya-watson estimator?
  Eq2hat1<-colSums(sweep(parsamp^2,1,w1,"*"))/sw1 #importance weights attached, each draw has weight proportional to exp(iw+kw)
  Varqhat1<-Eq2hat1-Eqhat1^2 #variance, a single value, an estimate of the variance b  
  postmarcdf1<-vector(mode="list",length=p)
  postmarden1<-vector(mode="list",length=p)
  pi991<-matrix(nrow=2,ncol=p)
  pi951<-matrix(nrow=2,ncol=p)
  pi501<-matrix(nrow=2,ncol=p)
  for(i in 1:p){
    ord<-order(parsamp[,i])
    y<-cumsum(nw1[ord])
    x<-parsamp[ord,i]
    sf<-splinefun(x=x,y=y,method="monoH.FC")
    postmarcdf1[[i]]<-list(theta=x,F=y)
    postmarden1[[i]]<-list(theta=x,f=sf(x=seq(from=min(x),to=max(x),length=1e3),deriv=1))
    pi991[1,i]<-findQ(0.005,x,y,sf)
    pi991[2,i]<-findQ(0.995,x,y,sf)
    pi951[1,i]<-findQ(0.025,x,y,sf)
    pi951[2,i]<-findQ(0.975,x,y,sf)
    pi501[1,i]<-findQ(0.25,x,y,sf)
    pi501[2,i]<-findQ(0.75,x,y,sf)
  }
  cat("Returning results.\n")
  vnam<-names(smu)
  names(Eqhat)<-vnam
  names(Varqhat)<-vnam
  names(Eqhat1)<-vnam
  names(Varqhat1)<-vnam
  names(prior.mean)<-vnam
  names(prior.sd)<-vnam
  names(postmarcdf)<-vnam
  names(postmarden)<-vnam
  colnames(pi99)<-vnam
  rownames(pi99)<-c("Q0.005","Q0.995")
  colnames(pi95)<-vnam
  rownames(pi95)<-c("Q0.025","Q0.975")
  colnames(pi50)<-vnam
  rownames(pi50)<-c("Q0.25","Q0.75")
  colnames(pi991)<-vnam
  rownames(pi991)<-c("Q0.005","Q0.995")
  colnames(pi951)<-vnam
  rownames(pi951)<-c("Q0.025","Q0.975")
  colnames(pi501)<-vnam
  rownames(pi501)<-c("Q0.25","Q0.75")
  out<-list()
  out$post.mean<-Eqhat
  out$post.sd<-sqrt(Varqhat)
  out$post.cov<-Covqhat
  out$post.mean.old<-Eqhat1
  out$post.sd.old<-sqrt(Varqhat1)
  out$post.interval.99<-pi99
  out$post.interval.95<-pi95
  out$post.interval.50<-pi50
  out$post.interval.991<-pi991
  out$post.interval.951<-pi951
  out$post.interval.501<-pi501
  #out$post.interval.90<-pi90
  #out$post.interval.99.vp <- pi99_vp
  #out$post.interval.95.vp <- pi95_vp
  #out$post.interval.50.vp <- pi50_vp
  out$post.marginal.cdf<-postmarcdf
  out$post.marginal.density<-postmarden
  #out$prob.ge.0<-pge0
  out$log.kernel.weights<-kw
  out$log.imp.weights<-iw
  out$nw <- nw
  out$nw1 <- nw1
  out$sw <- sw
  out$d <- d
  out$IS.eff.samp.size <- eff_samp_size
  out$param.samp<-parsamp
  out$param.samp.resample <- parsamp_resample
  out$resample<-resample
  out$target.stats<-target.stats #ostats equals to target stats if transf=none
  out$ostats <- ostats
  out$samp.dist<-d 
  out$samp.stats<-samp.stats
  out$stats <- stats
  out$transf <- transf
  out$initial.fit<-ifit
  out$samp.scale.array<-samp.scale.array
  out$samp.mu<-samp.mu
  out$samp.df<-samp.df
  out$samp.size<-samp.size
  out$bandwidth<-bw
  out$bandwidth.fun<-as.character(expression(bw.fun))
  out$prior.mean<-prior.mean
  out$prior.sd<-prior.sd
  out$formula<-form
  out$fittime <- difftime(runtime2, runtime1)
  out$samptime <- difftime(runtime3, runtime2)
  out$runtime <- difftime(runtime4, runtime1)
  out$qtime <- difftime(runtime4,runtime3)  
  out$seed.base <- seed.base
  #out$post.samples.copula <- post_samples_copula
  class(out)<-"kbergm"
  out
}

print.kbergm<-function(x,...){
  cat("Kernel Bayes ERGM Fit\n\n")
  cat("Posterior mean estimates:\n")
  print(x$post.mean)
}

summary.kbergm<-function(object,...){
  class(object)<-c("summary.kbergm",class(object))
  object
}

print.summary.kbergm<-function(x,...){
  cat("=================================\n")
  cat("Bayesian ERGM w/Kernel Estimation\n")
  cat("=================================\n\n")
  cat("Formula: ",as.character(x$formula),"\n\n")
  cat("Kernel estimates:\n")
  ctab<-cbind(x$post.mean,x$post.sd,x$prob.ge.0,1-x$prob.ge.0, pmin(x$prob.ge.0,1-x$prob.ge.0))
  rownames(ctab)<-names(x$post.mean)
  colnames(ctab)<-c("Post.Mean","Post.SD","Pr(q>=0)","Pr(q<0)","Pr(sign(q)!=signhat)")
  printCoefmat(ctab)
  cat("\nPrior means:\n")
  print(x$prior.mean)
  cat("Prior SDs:\n")
  print(x$prior.sd)
  cat("\nSample draws:",x$samp.size,"Effective sample size:",1/sum(exp(x$log.kernel.weights+x$log.imp.weights)),"\n")
  cat("Sample mean parameters:\n")
  print(x$samp.mu)
  cat("Stat distance quantiles (vs. target statistics):\n")
  print(quantile(x$samp.dist))
  cat("Bandwidth:",x$bandwidth,"\n\n")
}

plot.kbergm<-function(object,...){
  op<-par(no.readonly=TRUE)
  ask<-devAskNewPage()
  n<-length(object$post.mean)  #dimension of the model parameters
  vnam<-names(object$post.mean)
  par(mfrow=n2mfrow(n))
  for(i in 1:n){
    d<-abs(object$param.samp[,i]-object$post.mean[i])
    d<-d/max(d)
    plot(object$param.samp[,i],object$log.imp.weights,cex=0.2,col=rgb(0,0,0,0.5),xlab=expression(hat(theta)), ylab="Log Importance Weights",main=vnam[i])
    points(object$param.samp[,i],object$log.imp.weights,pch=19,cex=0.2,col=hsv(d*0.65, alpha=exp(object$log.kernel.weights-max(object$log.kernel.weights))^0.5))
    abline(v=object$post.mean[i],lty=1,col=1)
    abline(v=object$samp.mu[i],lty=2,col=2)
    abline(v=object$prior.mean[i],lty=3,col=3)
  }
  devAskNewPage(TRUE)
  for(i in 1:n){
    q<-object$post.marginal.cdf[[i]]$theta
    f<-object$post.marginal.cdf[[i]]$F
    plot(0,0,type="n",xlim=range(q),ylim=c(0,1),xlab=expression(theta), ylab=expression(paste("Pr",group("(",paste(theta<=x,group("|",Y,".")),")"))),main=vnam[i])
    abline(v=object$post.interval.99[,i],col=2,lty=3)
    abline(v=object$post.interval.95[,i],col=2,lty=2,lwd=2)
    abline(v=object$post.mean[i],col=2,lwd=2)
    segments(min(q)-diff(range(q)),0,q[1],f[1])
    segments(max(q)+diff(range(q)),1,q[length(q)],f[length(q)])
    lines(q,f,lwd=2)
  }
  devAskNewPage(TRUE)
  for(i in 1:n){
    x<-object$param.samp[,i]
    nw<-object$nw
    plot(seq(min(x), max(x), length=length(object$post.marginal.density[[i]]$f)), wkde(x,seq(min(x), max(x), length=length(object$post.marginal.density[[i]]$f)),w=nw),
     xlab=expression(theta), ylab="kernel density estimate of weigted sample", main=vnam[i],type="l")    
  }
  if(object$resample){
  devAskNewPage(TRUE)
  for(i in 1:n){
    x<-object$param.samp.resample[,i]
    plot(density(x), xlab=expression(theta), ylab="kernel density estimate of resampling importance sampling", main=vnam[i],type="l")
    #nw<-object$nw
    #plot(seq(min(x), max(x), length=length(object$post.marginal.density[[i]]$f)), density(x),
    #     xlab=expression(theta), ylab="kernel density estimate of resampling importance sampling", main=vnam[i],type="l")    
  } 
 }


  devAskNewPage(ask)
  par(op)
}

#######################################################################################################
############
##################################
#####input : dat : which should be a network object
#####output: edge_id matrix, null_id matrix
##################################
edge.null.id <- function(dat)
{
  #check data, if it is a network object
  if(! (network::is.network(dat))){
    stop("Not a network object!")
  }
  
  #pull out the adjacency matrix
  adj_mat <- dat[,]
  diag(adj_mat) <- rep(NA, length(diag(adj_mat))) #set the diagonal elements to be zero so that, we assume those networks are simple

  #edge_id_matrix <- igraph::get.edgelist(asIgraph(dat)) #sorry, this function is sooo convenient   
  adj_mat <- dat[,]
  diag(adj_mat) <- rep(NA, length(diag(adj_mat)))
      
      if(network::is.directed(dat)){
        null_id_matrix <- which(adj_mat == 0, arr.ind = T)  #pull out the id of absent edges
        edge_id_matrix <- which(adj_mat == 1, arr.ind = T)  #pull out the id of present edges
      }
      else {
        adj_mat_upper <- lower.tri.remove(adj_mat) #the way we set up the 
        #dyad_id, upper triangle of the adjacency matrix
        edge_id_matrix <- which(adj_mat_upper == 1, arr.ind = T)
        null_id_matrix <- which(adj_mat_upper == 0, arr.ind = T) }


   return(list(edge_id_matrix = edge_id_matrix, null_id_matrix = null_id_matrix))
}
#############
########################################################################################################
#######################################################################################################
########### abcergm : approximate bayesian computation for ERGMs, IS
########### need to define metrics...
#######################################################################################################
abcergm <- function(form,prior.mean=0,prior.sd=10,samp.size=1e3,samp.df=4,samp.scale.fact=4, dist.type=c("euclidean","mahalanobis"),bw.fun=bw.nrd0,mc.cores=1,splinefun.method = "monoH.FC",MPLE.type=c("glm","penalized"),threshold = 0.01,...){
  require(doParallel)
  require(sna)
  require(ergm)
  #
  #Fit the MPLE
  cat("Performing initial MPLE fit...\n")
  if(!is.null(MPLE.type)){
   ifit<-ergm(form,estimate="MPLE", control=control.ergm(MPLE.type = MPLE.type)) }
  else{
   ifit<-ergm(form,control=control.ergm(main.method = "MCMLE"))
  }
  samp.size.1 <- samp.size[1] # sample size for first round
  ostats<-ifit$target.stats  #Get target stats
  #Draw the importance sample
  cat("Drawing the importance sample\n")
  sscale<-ifit$covar*samp.scale.fact
  smu<-ifit$coef
  parsamp<-rmvt(samp.size.1,mu=smu,cov=sscale,df=samp.df)
  #Obtain the importance ratio
  cat("Getting importance ratios\n")
  p<-length(smu)
  slprob<-dmvt(parsamp,mu=smu,cov=sscale,df=samp.df,log=TRUE)
  prior.mean<-rep(prior.mean,length=p)
  prior.sd<-rep(prior.sd,length=p)
  if(p>1)
    prior.sdm<-diag(prior.sd)
  else
    prior.sdm<-diag(1)*prior.sd
  plprob<-dmvnorm(parsamp,mu=prior.mean,cov=prior.sdm,log=TRUE)
  iw<-plprob-slprob
  #Draw stats
  cat("Drawing the statistics at the sampled parameters\n")
  if(is.null(mcmc.burnin)){
   stats<-t(matrix(unlist(mclapply(1:samp.size.1,function(z){
    simulate(form,coef=parsamp[z,],statsonly=TRUE,nsim=1,...)},mc.cores=mc.cores)),ncol=samp.size.1))  } else{
      stats<-t(matrix(unlist(mclapply(1:samp.size.1,function(z){
        simulate(form,coef=parsamp[z,],statsonly=TRUE,nsim=1,control=control.simulate.formula(MCMC.burnin = mcmc.burnin,MCMC.interval=mcmc.interval),...)  },mc.cores=mc.cores)),ncol=samp.size.1)) }
  #Get distances from the target, and kernel weights
  cat("Obtaining distances to target, and kernel weights\n")
  if(match.arg(dist.type)=="euclidean"){
    d<-rowSums(sweep(stats,2,ostats,"-")^2)^0.5
  }else{
    d<-rep(0,samp.size.1)
    iscov<-solve(cov(stats))
    for(i in 1:samp.size.1)
      d[i]<-sqrt(t(stats[i,]-ostats)%*%iscov%*%(stats[i,]-ostats))
  }
  if(!is.null(threshold)){
    # only want to keep the samples that are sufficiently close to the observed stats, threshold corresponds to quantile
    keep_index <- (d <= quantile(d,threshold)) # d is the mahalanobis distance between observed stats and simulated stats
    iw<-iw[keep_index]
    parsamp<-parsamp[keep_index,]    
    stats<-stats[keep_index,]
    d<-d[keep_index]
  }

  ###
  niw <- exp(iw)/sum(exp(iw))
  ###
  out$post.mean <- parsamp %>% apply(2, function(x) reldist::wtd.mean(x, niw) )
  out$post.sd <- parsamp %>% apply(2, function(x) sqrt(reldist::wtd.var(x, niw)) )
}
##########################################################################################################
###########
###########
########################################################################################################
stats.param.top <- function(x,top=20, ...){
  # top weights
  top_nw <- head(sort(x$nw,decreasing=TRUE), top)
  top_index <- match(top_nw,x$nw)
  
  #print weights
  cat("The top", top, "weights are \n")
  print(top_nw)
   
  #print stats
  cat("The sampled statistics associated with the top", top, "weights\n")
  print(x$samp.stats[top_index,])
  cat("\n")

  cat("The sampled parameter associated with the top", top, "weights\n")
  print(x$param.samp[top_index,])
  cat("\n")

  return(list(top_nw = top_nw, top_samp_stats = x$samp.stats[top_index,], top_param_samp = x$param.samp[top_index,] ))
}
###########################################################################################################

########################################
################## CV.dist.bw
################## Using cross validation for bandwidth selection of kernel regression 
################## for kernel bayes rule of ERGMs
########################################

########################################
### input, d, y, kernel function (by default a probability density function)
### hmin, hmax
### nbh (number of candidate bandwidths)
### importance weights
CV.dist.bw <- function(d, y, log_iw = rep(0, length(d)), hmin = 0.001, hmax = 0.3, nbh = 50){
  library(sna)
  n <- length(d)
  vecth <- seq(from = hmin, to = hmax, length = nbh)
  matCV <- cbind(vecth, rep(0,nbh))
  log_iw<-log_iw-logSum(log_iw)
  # go over all possible selections
  for(h in 1:nbh){
    ypred <- rep(0,n)
    log_kw<-dnorm(d/vecth[h],log=TRUE) #bandwidth does not change in the inner loop
    log_kw<-log_kw-logSum(log_kw)
    for(i in 1:n){
      nw<-exp(log_iw[-i]+log_kw[-i])#normalized weight for kernel regression
      nw<-nw/sum(nw) # normalized weight
      ypred[i] <- sum(nw * y[-i])
    }
    matCV[h,2] <- sum( (y-as.matrix(ypred,ncol=1))^2 )
  }
  hopt<-matCV[which(matCV[,2]==min(matCV[,2],na.rm=TRUE))]
  return(list(matCV = matCV, hopt=hopt))
}

################################################################################################################
#############################
##### function for improved SIR with replacement (Skare 2003)
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
##################################################################################################################
###############################
######## kNN CV bandwidth selection
######## input : distance, target stats, simulated parameters, simulated statistics, hmin, hmax, nbh 
######## 
##################################################################################################################
###############################
kNN.CV.bw <- function(k=50,d,kernel.fun, method=c("optim","grid"),params, hmin, hmax, nbh,num_cores=1){
  library(foreach)
  library(doParallel)
  registerDoParallel(cores = num_cores)
  # kNN cross validation bandwidth selection, because we only care about the prediction in the
  # region is close to the observed statistics
  kNN_id <- which( d <= sort(d,decreasing = FALSE)[k] )
  if(match.arg(method) == "grid"){
   vecth <- seq(hmin, hmax, length=nbh)
   matCV <- cbind(vecth, rep(0,nbh))  # 
   # loop over all possible combinations of 
   for(h in 1:nbh){
     #ypred <- rep(NA, k) #only care about the prediction near observed stats
     log_kw<-kernel.fun(d/vecth[h], log = TRUE) # assume that mahalanobis distance does not change for leave-1-out
     #browser()
     ypred <- foreach(i = 1:k) %dopar% {
       kw<-exp(log_kw[-kNN_id[i]])
       nw<-kw/sum(kw)
       sum(nw * params[-kNN_id[i]]) #params is a vector
     }
     ypred <- do.call(c,ypred) #convert vector of lists to vector
     matCV[h,2] <- sum( (ypred - params[kNN_id])^2 )
   }
   hopt = vecth[which(matCV[,2] == min(matCV[,2], na.rm=TRUE))] }
  else{
    # need a function to define the loss given a bandwidth
    kNN.CV.loss <- function(h){ # h is the bandwidth
       log_kw<-kernel.fun(d/h, log = TRUE) # assume that mahalanobis distance does not change for leave-1-out
       #browser()
       ypred <- foreach(i = 1:k) %dopar% {
        kw<-exp(log_kw[-kNN_id[i]])
        nw<-kw/sum(kw)
        sum(nw * params[-kNN_id[i]]) #params is a vector
       }
       ypred <- do.call(c,ypred) #convert vector of lists to vector
       return(sum( (ypred - params[kNN_id])^2 ))
    }
    hopt=optimize(kNN.CV.loss, interval=c(hmin,hmax))$minimum }
   return(list(hopt=hopt))
}
##############
#### epanechnikov function
##############
epanechnikov <- function(x,log=FALSE){
  y <- 3/4 * (1 - x^2) * (abs(x)<=1)
  if(log==TRUE){
    y<-log(y)
  }
  return(y)
}
###########################################################
###########################################################
# total.var.dist, function to calculate the total variation distance
# between two random samples from univariate distributions
###########################################################
# total.var.dist <- function(s1, s2, nbins=100){
#   
#   
#   
#   
# }
##############################################################
# using fitCopula function to find a good estimate of the dependence parameter
# mvdc object
# temp_paramMargins = vector(mode="list", length = p)
# for(i in 1:p){
#  temp_paramMargins[[i]] <- list(mean = 0, sd = 1)
# }
# Mvd<-mvdc(copula=ellipCopula(family="normal", param=0.5), margins = rep("norm", p), paramMargins = temp_paramMargins )  # we only need the copula slot to be filled
# starting value 
# a.0 <- matrix(1, nrow = p, ncol = p)
# for(i in 1:(p-1)){
#  for(j in (i+1):p){
#     a.0[i,j] <- weightedCorr(parsamp[,i], parsamp[,j], weights = nw) } } #starting value of fit coupla
# assign values to the lower triangle of a.0
# pseduo observation
###############################################################
## assume each marginal is Gaussian and let's glue them together using t-copula
# mvdc object
#temp_paramMargins = vector(mode="list", length = p)
#for(i in 1:p){
#  temp_paramMargins[[i]] <- list(mean = Eqhat[i], sd = sqrt(Varqhat)[i])
#}
#a.0 <- matrix(1, nrow = p, ncol = p)
#for(i in 1:p){
#  for(j in 1:p){
#     a.0[i,j] <- wtd.cor(parsamp[,i], parsamp[,j], nw)[1] } } #starting value of fit coupla
#Mvd<-mvdc(copula=ellipCopula(family="t", dim = p, dispstr = "un", param = a.0[upper.tri(a.0)], df=4), margins = rep("norm", p), paramMargins = temp_paramMargins )  # we only need the copula slot to be filled
# starting value     
#post_samples_copula <- rMvdc(samp.size, Mvd)
#Convert pseudo-obs to logit scale
############ old codes
# lpsobs<-pmax(pmin(psobs,1e-6),1-1e-6)
# lpsobs<-log(lpsobs/(1-lpsobs))


#
# if(est.type == "N-D"){
#   Eqhat<-colSums(sweep(parsamp,1,w,"*"))/sw
#   Eq2hat<-colSums(sweep(parsamp^2,1,w,"*"))/sw   #locally constant estimate of the function of posterior samples
# }
# else if(est.type == "GAM"){
#   for(i in 1:p){
#     temp_dat<-data.frame(y=parsamp[,i],stats)
#     temp_dat2<-data.frame(y=parsamp[,i]^2,stats)
#     predictorvariables <- paste("s(X", 1:p, ")", sep="")
#     temp_formula <- formula(paste("y ~ ", 
#                                   paste(predictorvariables, collapse=" + ")))
#     temp_Eqhat_gam <- mgcv::gam(temp_formula, family=gaussian(), data = temp_dat,weights=exp(iw)/mean(exp(iw)))
#     temp_Eq2hat_gam <- mgcv::gam(temp_formula, family=gaussian(), data = temp_dat2,weights=exp(iw)/mean(exp(iw)))
#     temp_newdat <- data.frame(y=0,matrix(ostats,nrow=1))
#     Eqhat[i] <- predict(temp_Eqhat_gam, temp_newdat)
#     Eq2hat[i] <- predict(temp_Eq2hat_gam, temp_newdat)
#   }
# }
################################################
############# gof based on joint posterior samples
################################################
gof.kbergm<-function(x,directed=FALSE,sample.size=100,aux.iters=10000,n.deg = NULL, n.dist = NULL, n.esp = NULL, n.ideg = NULL, 
                               n.odeg = NULL,...){
  #class(object)<-c("summary.kbergm",class(object))
  #object
  cat("Goodness of fit for Kbergm, Joing samples obtained based on sampling importance sampling\n\n")
  # function (x, directed = FALSE, sample.size = 100, aux.iters = 10000, 
  #           n.deg = NULL, n.dist = NULL, n.esp = NULL, n.ideg = NULL, 
  #           n.odeg = NULL, ...) 
  #{
    if (class(x) == "kbergm") {
      #if (x$nchains > 1) 
        x$param.samp.resample <- apply(x$param.samp.resample, 2, cbind)
      F <- as.matrix(x$param.samp.resample[sample(dim(x$param.samp.resample)[1], sample.size), 
                             ])
    #}
    #else {
      #F <- x$param.samp.resample
    }
    if (directed == FALSE) {
      for (i in 1:sample.size) {
        a <- gof(x$formula, coef = F[i, ], verbose = FALSE, 
                 control = control.gof.formula(nsim = 1, MCMC.burnin = aux.iters))
        if (i == 1) 
          A <- as.vector(a$pobs.deg)
        A <- cbind(A, as.vector(a$psim.deg))
        if (i == 1) 
          B <- as.vector(a$pobs.dist)
        B <- cbind(B, as.vector(a$psim.dist))
        if (i == 1) 
          C <- as.vector(a$pobs.espart)
        C <- cbind(C, as.vector(a$psim.espart))
      }
      if (is.null(n.deg)) 
        n.deg <- dim(A)[1]
      if (is.null(n.dist)) 
        n.dist <- dim(B)[1] - 1
      if (is.null(n.esp)) 
        n.esp <- dim(C)[1]
      a5 <- apply(A[1:n.deg, -1], 1, quantile, probs = 0.05)
      b5 <- apply(B[-(n.dist:(dim(B)[1] - 1)), -1], 1, quantile, 
                  probs = 0.05)
      c5 <- apply(C[1:n.esp, -1], 1, quantile, probs = 0.05)
      a95 <- apply(A[1:n.deg, -1], 1, quantile, probs = 0.95)
      b95 <- apply(B[-(n.dist:(dim(B)[1] - 1)), -1], 1, quantile, 
                   probs = 0.95)
      c95 <- apply(C[1:n.esp, -1], 1, quantile, probs = 0.95)
      par(mfrow = c(1, 3), oma = c(0, 0, 3, 0), mar = c(4, 
                                                        3, 1.5, 1))
      boxplot(as.data.frame(t(A[1:n.deg, -1])), xaxt = "n", 
              xlab = "degree", ylab = "proportion of nodes")
      axis(1, seq(1, n.deg), seq(0, n.deg - 1))
      lines(A[1:n.deg, 1], lwd = 2, col = 2)
      lines(a5, col = "darkgray")
      lines(a95, col = "darkgray")
      title("Bayesian goodness-of-fit diagnostics, Kernel ABC", outer = TRUE)
      boxplot(as.data.frame(t(B[-(n.dist:(dim(B)[1] - 1)), 
                                -1])), xaxt = "n", xlab = "minimum geodesic distance", 
              ylab = "proportion of dyads")
      axis(1, seq(1, n.dist), labels = c(seq(1, (n.dist - 1)), 
                                         "NR"))
      lines(B[-(n.dist:(dim(B)[1] - 1)), 1], lwd = 2, col = 2)
      lines(b5, col = "darkgray")
      lines(b95, col = "darkgray")
      boxplot(as.data.frame(t(C[1:n.esp, -1])), xaxt = "n", 
              xlab = "edge-wise shared partners", ylab = "proportion of edges")
      axis(1, seq(1, n.esp), seq(0, n.esp - 1))
      lines(C[1:n.esp, 1], lwd = 2, col = 2)
      lines(c5, col = "darkgray")
      lines(c95, col = "darkgray")
      out = list(sim.degree = A[, -1], sim.dist = B[, -1], 
                 sim.esp = C[, -1], obs.degree = A[, 1], obs.dist = B[, 
                                                                      1], obs.esp = C[, 1], fun = F)
    }
    else {
      for (i in 1:sample.size) {
        a <- gof(x$formula, coef = F[i, ], verbose = FALSE, 
                 GOF = ~idegree + odegree + espartners + distance, 
                 control = control.gof.formula(nsim = 1, MCMC.burnin = aux.iters))
        if (i == 1) 
          A <- as.vector(a$pobs.ideg)
        A <- cbind(A, as.vector(a$psim.ideg))
        if (i == 1) 
          AA <- as.vector(a$pobs.odeg)
        AA <- cbind(AA, as.vector(a$psim.odeg))
        if (i == 1) 
          B <- as.vector(a$pobs.dist)
        B <- cbind(B, as.vector(a$psim.dist))
        if (i == 1) 
          C <- as.vector(a$pobs.espart)
        C <- cbind(C, as.vector(a$psim.espart))
      }
      if (is.null(n.ideg)) 
        n.ideg <- dim(A)[1]
      if (is.null(n.odeg)) 
        n.odeg <- dim(AA)[1]
      if (is.null(n.dist)) 
        n.dist <- dim(B)[1] - 1
      if (is.null(n.esp)) 
        n.esp <- dim(C)[1]
      a5 <- apply(A[1:n.ideg, -1], 1, quantile, probs = 0.05)
      aa5 <- apply(AA[1:n.odeg, -1], 1, quantile, probs = 0.05)
      b5 <- apply(B[-(n.dist:(dim(B)[1] - 1)), -1], 1, quantile, 
                  probs = 0.05)
      c5 <- apply(C[1:n.esp, -1], 1, quantile, probs = 0.05)
      a95 <- apply(A[1:n.ideg, -1], 1, quantile, probs = 0.95)
      aa95 <- apply(AA[1:n.odeg, -1], 1, quantile, probs = 0.95)
      b95 <- apply(B[-(n.dist:(dim(B)[1] - 1)), -1], 1, quantile, 
                   probs = 0.95)
      c95 <- apply(C[1:n.esp, -1], 1, quantile, probs = 0.95)
      par(mfrow = c(2, 2), oma = c(0, 0, 3, 0), mar = c(4, 
                                                        3, 1.5, 1))
      boxplot(as.data.frame(t(A[1:n.ideg, -1])), xaxt = "n", 
              xlab = "in degree", ylab = "proportion of nodes")
      axis(1, seq(1, n.ideg), seq(0, n.ideg - 1))
      lines(A[1:n.ideg, 1], lwd = 2, col = 2)
      lines(a5, col = "darkgray")
      lines(a95, col = "darkgray")
      title("Bayesian goodness-of-fit diagnostics, Kernel ABC", outer = TRUE)
      boxplot(as.data.frame(t(AA[1:n.odeg, -1])), xaxt = "n", 
              xlab = "out degree", ylab = "proportion of nodes")
      axis(1, seq(1, n.odeg), seq(0, n.odeg - 1))
      lines(AA[1:n.odeg, 1], lwd = 2, col = 2)
      lines(aa5, col = "darkgray")
      lines(aa95, col = "darkgray")
      boxplot(as.data.frame(t(B[-(n.dist:(dim(B)[1] - 1)), 
                                -1])), xaxt = "n", xlab = "minimum geodesic distance", 
              ylab = "proportion of dyads")
      axis(1, seq(1, n.dist), labels = c(seq(1, (n.dist - 1)), 
                                         "NR"))
      lines(B[-(n.dist:(dim(B)[1] - 1)), 1], lwd = 2, col = 2)
      lines(b5, col = "darkgray")
      lines(b95, col = "darkgray")
      boxplot(as.data.frame(t(C[1:n.esp, -1])), xaxt = "n", 
              xlab = "edge-wise shared partners", ylab = "proportion of edges")
      axis(1, seq(1, n.esp), seq(0, n.esp - 1))
      lines(C[1:n.esp, 1], lwd = 2, col = 2)
      lines(c5, col = "darkgray")
      lines(c95, col = "darkgray")
      out = list(sim.idegree = A[, -1], sim.odegree = AA[, 
                                                         -1], sim.dist = B[, -1], sim.esp = C[, -1], obs.degree = A[, 
                                                                                                                    1], obs.dist = B[, 1], obs.esp = C[, 1])
    }
}


















###############################################################
########### kbergm without full access to sufficient statistics
###############################################################
###########
########### Using degree distribution as surrogate
# kbergm.ego <- function(form,prior.mean=0,prior.sd=10,samp.size=1e3,samp.df=4,samp.scale.fact=4,est.type="N-D",
#                        transf=c("none", "log1p", "sqrt1"),dist.type=c("euclidean","mahalanobis"),
#                        bw.fun=bw.nrd0,kernel.fun = dnorm,kernel.type=c("gaussian","epanechnikov"),
#                        bw.type=c("single","elementwise","percentage"),percentage=0.01,mc.cores=1,
#                        splinefun.method = "monoH.FC",resample=FALSE,resample.size=500,resample.replace=FALSE,mcmc.burnin=NULL, 
#                        mcmc.interval=NULL,MPLE.type=c("glm","penalized"),
#                        threshold=0.1,imp.rounds=1,k=50,
#                        cv.method="optim",...){
#   
#   
#   
#   
#   
#   
#   
#   
# }
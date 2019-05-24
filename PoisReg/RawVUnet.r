library(rjags)
library(coda)
library(png)
library(rgdal)
library(raster)
library(magick)
library(bayesplot)
library(doMC); registerDoMC(4)
library(rjags); load.module("glm"); load.module("lecuyer")
library(random)
library(ggplot2)
load.module("dic")

RMSE <- function(x, xhat){
  return(sqrt(mean((x-xhat)^2)))
}
z_standardise <- function(x,xbar, s){
  return((x-xbar)/s)
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
MCMC_mode = function(samps){
  d <- density(samps)
  return(d$x[which.max(d$y)])
}

jagsresults <- function(x, params, regex=FALSE, invert=FALSE, 
                        probs=c(0.025, 0.25, 0.5, 0.75, 0.975), signif, ...) {
  if(!regex) {
    params <- paste0(gsub('(?=\\.|\\[|\\])', '\\1\\\\', params, perl=TRUE),
                     '(\\[.*\\])?', collapse='|')
    params <- paste("^", gsub("\\|", "$|^", params), '$', sep = "")
  } else if(length(params) > 1) {
    stop("If 'regex' is TRUE, 'params' must be a single regex string.",
         call.=FALSE)
  }
  if(any(is(x) %in% c('rjags.parallel', 'rjags'))) {
    nm <- dimnames(x$BUGSoutput$sims.array)[[3]]
    i <- grep(params, nm, invert=invert, ...)
    if(length(i) == 0) stop("No parameters match 'params'", call.=FALSE)
    samp <- x$BUGSoutput$sims.array[, , i, drop=FALSE] 
    rhat_neff <- x$BUGSoutput$summary[i, c('Rhat', 'n.eff'), drop=FALSE]
    mode <- MCMC_mode(samp)
    out <- cbind(t(apply(
      samp, 3, function(x) 
        c(mean=mean(x), mode,sd=sd(x), quantile(x, probs=probs)))), rhat_neff)
  } else if(any(is(x)=='mcmc.list')) {
    nm <- colnames(x[[1]])
    i <- grep(params, nm, invert=invert, ...)
    if(length(i) == 0) stop("No parameters match 'params'", call.=FALSE)
    out <- t(apply(do.call(rbind, x), 2, function(z) {
      c(mean=mean(z), sd=sd(z), quantile(z, probs))
    }))[i, , drop=FALSE]
  } else {
    stop("x must be an 'mcmc.list' or 'rjags'  object.")
  }
  if(!missing(signif)) {
    out[, colnames(out) != 'n.eff'] <- 
      signif(out[, colnames(out) != 'n.eff'], signif)  
    out
  } else out
}


HOME_DIR = 'Documents/ssc_challenge/ssc-case-study-2019/'
TRAIN_DIR = paste0(HOME_DIR, 'ssc-data/CompTrain/masks/')

setwd(HOME_DIR)
##############################################################################################################################################
########################### RAW DATA #######################################################################
##############################################################################################################################################
train_data_raw = read.table('./ssc-data/CompTrain/labels.csv', sep = ',', header = T, stringsAsFactors = F)
train_data_raw$sum_px = rep(0, nrow(train_data_raw))
train_data_raw$sum_px_s = rep(0, nrow(train_data_raw))

NTRAIN = nrow(train_data_raw)
for(i in 1:NTRAIN){
  img_file = train_data_raw[i,1]
  #img_file = gsub('.TIF', '.png', img_file, ignore.case = T)
  img = image_read(paste0('./ssc-data/CompTrain/', img_file))
  #img = image_scale(img, "128x128!")
  img = as.integer(img[[1]])/255
  #img = as.matrix(raster(as.character(paste0('./ssc-data/CompTrain/', img_file)))[,,1])/255
  train_data_raw$sum_px[i] = sum(img)
  
}
train_data_raw$stain = train_data_raw$stain-1
train_data_raw= train_data_raw[1:NTRAIN, ]

mean_sum_px = mean(train_data_raw$sum_px)
sd_sum_px = sd(train_data_raw$sum_px)
train_data_raw$sum_px_s <- sapply(train_data_raw$sum_px,function(x){z_standardise(x, mean_sum_px, sd_sum_px)})
train_data_raw$blur_lv <- 1
train_data_raw$blur_lv[train_data_raw$blur==23] <- 2
train_data_raw$blur_lv[train_data_raw$blur==48] <- 3


train_data_raw$group = 0
for(i in 1:NTRAIN){
  if(train_data_raw$stain[i]==0){
    if(train_data_raw$blur_lv[i]==1){
      train_data_raw$group[i] = 1
    }else if(train_data_raw$blur_lv[i]==2){
      train_data_raw$group[i]=2
    }else if(train_data_raw$blur_lv[i]==3){
      train_data_raw$group[i] = 3
    }
  }else{
    if(train_data_raw$blur_lv[i]==1){
      train_data_raw$group[i] = 4
    }else if(train_data_raw$blur_lv[i]==2){
      train_data_raw$group[i]=5
    }else if(train_data_raw$blur_lv[i]==3){
      train_data_raw$group[i] = 6
    }
  }
}


test_data_raw = read.table('./ssc-data/CompValid/labels.csv', sep = ',', header = T, stringsAsFactors = F)
test_data_raw$sum_px = rep(0, nrow(test_data_raw))
test_data_raw$sum_px_s = rep(0, nrow(test_data_raw))
NTEST = nrow(test_data_raw)
for(i in 1:NTEST){
  
  img_file = test_data_raw[i,1]
  #img_file = gsub('.TIF', '.png', img_file, ignore.case = T)
  img = image_read(paste0('./ssc-data/CompValid/', img_file))
  #img = image_scale(img, "128x128!")
  img = as.integer(img[[1]])/255
  #img = as.matrix(raster(as.character(paste0('./ssc-data/CompTrain/', img_file)))[,,1])/255
  test_data_raw$sum_px[i] = sum(img)
}
test_data_raw$stain = test_data_raw$stain-1
test_data_raw= test_data_raw[1:NTEST, ]

test_data_raw$sum_px_s <- sapply(test_data_raw$sum_px,function(x){z_standardise(x, mean_sum_px, sd_sum_px)})
test_data_raw$blur_lv <- 1
test_data_raw$blur_lv[test_data_raw$blur==23] <- 2
test_data_raw$blur_lv[test_data_raw$blur==48] <- 3

test_data_raw$group = 0
for(i in 1:NTEST){
  if(test_data_raw$stain[i]==0){
    if(test_data_raw$blur_lv[i]==1){
      test_data_raw$group[i] = 1
    }else if(test_data_raw$blur_lv[i]==2){
      test_data_raw$group[i]=2
    }else if(test_data_raw$blur_lv[i]==3){
      test_data_raw$group[i] = 3
    }
  }else{
    if(test_data_raw$blur_lv[i]==1){
      test_data_raw$group[i] = 4
    }else if(test_data_raw$blur_lv[i]==2){
      test_data_raw$group[i]=5
    }else if(test_data_raw$blur_lv[i]==3){
      test_data_raw$group[i] = 6
    }
  }
}


##############################################################################################################################################
########################### UNET MASKS DATA #######################################################################
##############################################################################################################################################
train_data_unet = read.table('./ssc-data/CompTrain/labels.csv', sep = ',', header = T, stringsAsFactors = F)
train_data_unet$sum_px = rep(0, nrow(train_data_unet))
train_data_unet$sum_px_s = rep(0, nrow(train_data_unet))

NTRAIN = nrow(train_data_unet)
for(i in 1:NTRAIN){
  img_file = train_data_unet[i,1]
  img_file = gsub('.TIF', '.png', img_file, ignore.case = T)
  img = image_read(paste0('./ssc-data/CompTrain/masks/', img_file))
  #img = image_scale(img, "128x128!")
  img = as.integer(img[[1]])/255
  #img = as.matrix(raster(as.character(paste0('./ssc-data/CompTrain/', img_file)))[,,1])/255
  train_data_unet$sum_px[i] = sum(img)
  
}
train_data_unet$stain = train_data_unet$stain-1
train_data_unet= train_data_unet[1:NTRAIN, ]

mean_sum_px = mean(train_data_unet$sum_px)
sd_sum_px = sd(train_data_unet$sum_px)
train_data_unet$sum_px_s <- sapply(train_data_unet$sum_px,function(x){z_standardise(x, mean_sum_px, sd_sum_px)})
train_data_unet$blur_lv <- 1
train_data_unet$blur_lv[train_data_unet$blur==23] <- 2
train_data_unet$blur_lv[train_data_unet$blur==48] <- 3


train_data_unet$group = 0
for(i in 1:NTRAIN){
  if(train_data_unet$stain[i]==0){
    if(train_data_unet$blur_lv[i]==1){
      train_data_unet$group[i] = 1
    }else if(train_data_unet$blur_lv[i]==2){
      train_data_unet$group[i]=2
    }else if(train_data_unet$blur_lv[i]==3){
      train_data_unet$group[i] = 3
    }
  }else{
    if(train_data_unet$blur_lv[i]==1){
      train_data_unet$group[i] = 4
    }else if(train_data_unet$blur_lv[i]==2){
      train_data_unet$group[i]=5
    }else if(train_data_unet$blur_lv[i]==3){
      train_data_unet$group[i] = 6
    }
  }
}


test_data_unet = read.table('./ssc-data/CompValid/labels.csv', sep = ',', header = T, stringsAsFactors = F)
test_data_unet$sum_px = rep(0, nrow(test_data_unet))
test_data_unet$sum_px_s = rep(0, nrow(test_data_unet))
NTEST = nrow(test_data_unet)
for(i in 1:NTEST){
  
  img_file = test_data_unet[i,1]
  img_file = gsub('.TIF', '.png', img_file, ignore.case = T)
  img = image_read(paste0('./ssc-data/CompValid/masks/', img_file))
  #img = image_scale(img, "128x128!")
  img = as.integer(img[[1]])/255
  #img = as.matrix(raster(as.character(paste0('./ssc-data/CompTrain/', img_file)))[,,1])/255
  test_data_unet$sum_px[i] = sum(img)
}
test_data_unet$stain = test_data_unet$stain-1
test_data_unet= test_data_unet[1:NTEST, ]

test_data_unet$sum_px_s <- sapply(test_data_unet$sum_px,function(x){z_standardise(x, mean_sum_px, sd_sum_px)})
test_data_unet$blur_lv <- 1
test_data_unet$blur_lv[test_data_unet$blur==23] <- 2
test_data_unet$blur_lv[test_data_unet$blur==48] <- 3

test_data_unet$group = 0
for(i in 1:NTEST){
  if(test_data_unet$stain[i]==0){
    if(test_data_unet$blur_lv[i]==1){
      test_data_unet$group[i] = 1
    }else if(test_data_unet$blur_lv[i]==2){
      test_data_unet$group[i]=2
    }else if(test_data_unet$blur_lv[i]==3){
      test_data_unet$group[i] = 3
    }
  }else{
    if(test_data_unet$blur_lv[i]==1){
      test_data_unet$group[i] = 4
    }else if(test_data_unet$blur_lv[i]==2){
      test_data_unet$group[i]=5
    }else if(test_data_unet$blur_lv[i]==3){
      test_data_unet$group[i] = 6
    }
  }
}



##############################################################################################################################################
########################### Model 1 Simple Polyn Reg... #######################################################################
##############################################################################################################################################

#1.a - RAW DATA
jags_input_raw = list(N = nrow(train_data_raw), x1 = train_data_raw$sum_px_s, g = train_data_raw$group, y = train_data_raw$count,
                  NTest = nrow(test_data_raw), xtest1 = test_data_raw$sum_px_s, gtest = test_data_raw$group)
jags_file_name = './PoisReg/mdl1_SimplePolynReg.jags'
params_to_monitor =c('y_hat', 'eta_hat', 'beta0', 'beta1','beta2','sigma','loglik', 'WAIC')
jags_inits = list(y_true = jags_input_raw$y, beta0 = rnorm(1,0,10), beta1 = rnorm(1,0,10), beta2 = rnorm(1,0,10))
jags.parsamples_raw <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
  model.jags <- jags.model(file = jags_file_name, data =jags_input_raw, inits = jags_inits,
                           n.chain=1, n.adapt=1000)
  result <- coda.samples(model = model.jags, variable.names = params_to_monitor,n.iter = 8000)
  return(result)
}
samps_m_raw = as.matrix(jags.parsamples_raw)
#param = "beta2"
#hist(samps_m_raw[,colnames(samps_m_raw) == param ], breaks = 100, xlab = param, main = param, freq = F)
#mcmc_trace(jags.parsamples_raw, pars=c('sigma'))
results_raw = data.frame(jagsresults(jags.parsamples_raw, 'y_hat'))
#RMSE(x = test_data_raw$count, results$mean)
print(paste('RMSE for Raw data = ', RMSE(x = test_data_raw$count, results_raw$mean)))
test_data_raw$yhat = results_raw$mean

loglik <- data.frame(samps_m_raw) %>% select(contains('loglik'))
loglik <- as.matrix(loglik)
loo(loglik, cores = getOption("mc.cores", 4))
#1.b - Unet
jags_input_unet = list(N = nrow(train_data_unet), x1 = train_data_unet$sum_px_s, g = train_data_unet$group, y = train_data_unet$count,
                       NTest = nrow(test_data_unet), xtest1 = test_data_unet$sum_px_s, gtest = test_data_unet$group)
jags_file_name = './PoisReg/mdl1_SimplePolynReg.jags'
params_to_monitor =c('y_hat', 'eta_hat', 'beta0', 'beta1','beta2','sigma','loglik', 'WAIC')
jags_inits = list(y_true = jags_input_unet$y, beta0 = rnorm(1,0,10), beta1 = rnorm(1,0,10), beta2 = rnorm(1,0,10))
jags.parsamples_unet <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
  model.jags <- jags.model(file = jags_file_name, data =jags_input_unet, inits = jags_inits,
                           n.chain=1, n.adapt=1000)
  result <- coda.samples(model = model.jags, variable.names = params_to_monitor,n.iter = 8000)
  return(result)
}
samps_m_unet = as.matrix(jags.parsamples_unet)
#param = "y_hat[142]"
#int.hist(samps_m_unet[,colnames(samps_m_unet) == param ], xlab = param, main = param, freq = F)
#mcmc_trace(jags.parsamples_unet, pars=c('sigma'))
results_unet = data.frame(jagsresults(jags.parsamples_unet, 'y_hat'))
print(paste('RMSE for UNET data = ', RMSE(x = test_data_unet$count, results_unet$mean)))
test_data_unet$yhat = results_unet$mean

loglik <- data.frame(samps_m_unet) %>% select(contains('loglik'))
loglik <- as.matrix(loglik)
loo(loglik, cores = getOption("mc.cores", 4))

##############################################################################################################################################
########################### Model 2 Varying SLopes... #######################################################################
##############################################################################################################################################

#1.a - RAW DATA
jags_input_raw = list(N = nrow(train_data_raw), x1 = train_data_raw$sum_px_s, g = train_data_raw$group, y = train_data_raw$count,
                      NTest = nrow(test_data_raw), xtest1 = test_data_raw$sum_px_s, gtest = test_data_raw$group)
jags_file_name = './PoisReg/mdl2_VarySlope.jags'
params_to_monitor =c('y_hat', 'eta_hat', 'beta0', 'beta1','beta2','sigma','loglik', 'WAIC')
jags_inits = list(y_true = jags_input_raw$y, beta0 = rnorm(1,0,10), beta1 = rnorm(6,0,10), beta2 = rnorm(6,0,10))
jags.parsamples_raw <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
  model.jags <- jags.model(file = jags_file_name, data =jags_input_raw, inits = jags_inits,
                           n.chain=1, n.adapt=1000)
  result <- coda.samples(model = model.jags, variable.names = params_to_monitor,n.iter = 8000)
  return(result)
}
samps_m_raw = as.matrix(jags.parsamples_raw)
#param = "beta2"
#hist(samps_m_raw[,colnames(samps_m_raw) == param ], breaks = 100, xlab = param, main = param, freq = F)
#mcmc_trace(jags.parsamples_raw, pars=c('sigma'))
results_raw = data.frame(jagsresults(jags.parsamples_raw, 'y_hat'))
#RMSE(x = test_data_raw$count, results$mean)
print(paste('RMSE for Raw data = ', RMSE(x = test_data_raw$count, results_raw$mean)))
test_data_raw$yhat = results_raw$mean
loglik <- data.frame(samps_m_raw) %>% select(contains('loglik'))
loglik <- as.matrix(loglik)
loo(loglik, cores = getOption("mc.cores", 4))

#2.b - Unet
jags_input_unet = list(N = nrow(train_data_unet), x1 = train_data_unet$sum_px_s, g = train_data_unet$group, y = train_data_unet$count,
                       NTest = nrow(test_data_unet), xtest1 = test_data_unet$sum_px_s, gtest = test_data_unet$group)
jags_file_name = './PoisReg/mdl2_VarySlope.jags'
params_to_monitor =c('y_hat', 'eta_hat', 'beta0', 'beta1','beta2','sigma','loglik', 'WAIC')
jags_inits = list(y_true = jags_input_unet$y, beta0 = rnorm(1,0,10), beta1 = rnorm(6,0,10), beta2 = rnorm(6,0,10))
jags.parsamples_unet <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
  model.jags <- jags.model(file = jags_file_name, data =jags_input_unet, inits = jags_inits,
                           n.chain=1, n.adapt=1000)
  result <- coda.samples(model = model.jags, variable.names = params_to_monitor,n.iter = 8000)
  return(result)
}
samps_m_unet = as.matrix(jags.parsamples_unet)
#param = "y_hat[142]"
#int.hist(samps_m_unet[,colnames(samps_m_unet) == param ], xlab = param, main = param, freq = F)
#mcmc_trace(jags.parsamples_unet, pars=c('sigma'))
results_unet = data.frame(jagsresults(jags.parsamples_unet, 'y_hat'))
print(paste('RMSE for UNET data = ', RMSE(x = test_data_unet$count, results_unet$mean)))
test_data_unet$yhat = results_unet$mean


loglik <- data.frame(samps_m_unet) %>% select(contains('loglik'))
loglik <- as.matrix(loglik)
loo(loglik, cores = getOption("mc.cores", 4))

+ ##############################################################################################################################################
########################### Model 3 Varying Slopes and Intcpt... #######################################################################
##############################################################################################################################################

#3.a - RAW DATA
jags_input_raw = list(N = nrow(train_data_raw), x1 = train_data_raw$sum_px_s, g = train_data_raw$group, y = train_data_raw$count,
                  NTest = nrow(test_data_raw), xtest1 = test_data_raw$sum_px_s, gtest = test_data_raw$group)
jags_file_name = './PoisReg/mdl3_VarySlopVaryIntcpt.jags'
params_to_monitor =c('y_hat', 'eta_hat', 'beta0', 'beta1','beta2','sigma','loglik', 'WAIC')
jags_inits = list(y_true = jags_input_raw$y, beta0 = rnorm(6,0,10), beta1 = rnorm(6,0,10), beta2 = rnorm(6,0,10))
jags.parsamples_raw <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
  model.jags <- jags.model(file = jags_file_name, data =jags_input_raw, inits = jags_inits,
                           n.chain=1, n.adapt=1000)
  result <- coda.samples(model = model.jags, variable.names = params_to_monitor,n.iter = 8000)
  return(result)
}
samps_m_raw = as.matrix(jags.parsamples_raw)
param = "beta1[1]"
hist(samps_m_raw[,colnames(samps_m_raw) == param ], breaks = 100, xlab = param, main = param, freq = F)
mcmc_trace(jags.parsamples_raw, pars=c('sigma'))
results = data.frame(jagsresults(jags.parsamples_raw, 'y_hat'))
RMSE(x = test_data_raw$count, results$mean)
test_data_raw$yhat = results$mean
test_data_raw$lwr = results$X2.5.
test_data_raw$upr = results$X97.5.  

loglik <- data.frame(samps_m_raw) %>% select(contains('loglik'))
loglik <- as.matrix(loglik)
loo(loglik, cores = getOption("mc.cores", 4))

plot(test_data_raw$count, test_data_raw$count-test_data_raw$yhat)

# check coverage
true_val_contained = (test_data_raw$count>=test_data_raw$lwr) & (test_data_raw$count<=test_data_raw$upr)
sprintf('95percent coverage is %.4f', mean(true_val_contained))

mean(test_data_raw$upr - test_data_raw$lwr)
#3.b - Unet
  jags_input_unet = list(N = nrow(train_data_unet), x1 = train_data_unet$sum_px_s, g = train_data_unet$group, y = train_data_unet$count,
                        NTest = nrow(test_data_unet), xtest1 = test_data_unet$sum_px_s, gtest = test_data_unet$group)
  jags_file_name = './PoisReg/mdl3_VarySlopVaryIntcpt.jags'
  params_to_monitor =c('y_hat', 'eta_hat', 'beta0', 'beta1','beta2','sigma', 'loglik', 'WAIC')
  jags_inits = list(y_true = jags_input_unet$y, beta0 = rnorm(6,0,10), beta1 = rnorm(6,0,10), beta2 = rnorm(6,0,10))
  jags.parsamples_unet <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
    model.jags <- jags.model(file = jags_file_name, data =jags_input_unet, inits = jags_inits,
                             n.chain=1, n.adapt=1000)
    result <- coda.samples(model = model.jags, variable.names = params_to_monitor, n.iter = 8000)
    return(result)
  }
  samps_m_unet = as.matrix(jags.parsamples_unet)
  param = "WAIC"
  int.hist(samps_m_unet[,colnames(samps_m_unet) == param ], xlab = param, main = param, freq = F)
  mcmc_trace(jags.parsamples_unet, pars=c('sigma'))
  results = data.frame(jagsresults(jags.parsamples_unet, 'y_hat'))
  RMSE(x = test_data_unet$count, results$mean)
  test_data_unet$yhat = results$mean


loglik <- data.frame(samps_m_unet) %>% select(contains('loglik'))
loglik <- as.matrix(loglik)
loo(loglik)

colnames(samps_m_unet)


### Posterior Predictive Plot
int.hist = function(x,ylab="Frequency",...) {
  barplot(table(factor(x,levels=min(x):max(x))),space=0,xaxt="n",ylab=ylab,...);axis(1)
}
img_name = 'H04_C14_F23_s16_w2.TIF'
idx =which(test_data_unet$image_name==img_name)
samps_of_interest_unet = samps_m_unet[,colnames(samps_m_unet) == paste0('y_hat[',idx,']')]
CI_unet = quantile(samps_of_interest_unet, probs = c(0.025, 0.975))

samps_of_interest_raw = samps_m_raw[,colnames(samps_m_raw) == paste0('y_hat[',idx,']')]
CI_raw = quantile(samps_of_interest_raw, probs = c(0.025, 0.975))

hist(samps_of_interest_unet, xlab = 'Cell Count', ylab = 'Density', main = NA, freq =F, probability = T,
     xaxs = "i",yaxs="i", zero.line = F, border = rgb(227/255,12/255,0/255,0.5), col = rgb(244/255,135/255,32/255,0.5),xlim = c(6,20))
segments(CI_unet[1], 0.01, CI_unet[2], 0.01, lwd=3, col=rgb(141/255,130/255,122/255,1))
points(x = mean(samps_of_interest_unet), y = 0.01, cex=2, pch = 15, col = rgb(141/255,130/255,122/255,1))

hist(samps_of_interest_raw, xlab = 'Cell Count', ylab = 'Density', main = NA, freq =F, probability = T,
     xaxs = "i",yaxs="i", zero.line = F, border = rgb(227/255,12/255,0/255,0.5), col = rgb(244/255,135/255,32/255,0.5),xlim = c(6,20))
segments(CI_raw[1], 0.01, CI_raw[2], 0.01, lwd=3, col=rgb(141/255,130/255,122/255,1))
points(x = mean(samps_of_interest_raw), y = 0.01, pch = 15,cex=2, col = rgb(141/255,130/255,122/255,1))


### Posterior Predictive Plot
par(mfrow=c(1,2))
param_results_raw <- data.frame(jagsresults(jags.parsamples_raw, c('beta0', 'beta1','beta2')))

plot(train_data_raw$sum_px_s[train_data_raw$group==1], train_data_raw$count[train_data_raw$group==1],
     ylim = c(0, 100), xlim = c(min(train_data_raw$sum_px_s),max(train_data_raw$sum_px_s)), xaxt='n',
     pch = 16, col = rgb(31/255,119/255,180/255,0.2), frame.plot=F,
     xlab = 'Sum of Pixels', ylab = 'Cell Count')
at =seq(min(train_data_raw$sum_px_s),max(train_data_raw$sum_px_s), by=0.5)
labels <- at*sd_sum_px + mean_sum_px
axis(side=1, at = at, labels = round(labels, 0))
points(train_data_raw$sum_px_s[train_data_raw$group==2], train_data_raw$count[train_data_raw$group==2],
    pch = 16, col =  rgb(255/255,127/255,14/255,0.2))
points(train_data_raw$sum_px_s[train_data_raw$group==3], train_data_raw$count[train_data_raw$group==3],
       pch = 16, col =  rgb(44/255,160/255,44/255,0.2))

points(train_data_raw$sum_px_s[train_data_raw$group==4], train_data_raw$count[train_data_raw$group==4],
     pch=17, col = rgb(31/255,119/255,180/255,0.2))
points(train_data_raw$sum_px_s[train_data_raw$group==5], train_data_raw$count[train_data_raw$group==5],
       pch = 17, col =  rgb(255/255,127/255,14/255,0.2))
points(train_data_raw$sum_px_s[train_data_raw$group==6], train_data_raw$count[train_data_raw$group==6],
       pch = 17, col =  rgb(44/255,160/255,44/255,0.2))

curve(param_results_raw[1,1] + param_results_raw[7,1]*x + param_results_raw[13,1]*x^2, 
      from = min(train_data_raw$sum_px_s), to =max(train_data_raw$sum_px_s),add = T,
      col =  rgb(31/255,119/255,180/255,1), lwd=2)
curve(param_results_raw[2,1] + param_results_raw[8,1]*x + param_results_raw[14,1]*x^2, 
      from = min(train_data_raw$sum_px_s), to =max(train_data_raw$sum_px_s),add = T,
      col =  rgb(255/255,127/255,14/255,1), lwd=2)
curve(param_results_raw[3,1] + param_results_raw[9,1]*x + param_results_raw[15,1]*x^2, 
      from = min(train_data_raw$sum_px_s), to =max(train_data_raw$sum_px_s),add = T,
      col =  rgb(44/255,160/255,44/255,1), lwd=2)

curve(param_results_raw[4,1] + param_results_raw[10,1]*x + param_results_raw[16,1]*x^2, 
      from = min(train_data_raw$sum_px_s), to =max(train_data_raw$sum_px_s),add = T,
      col =  rgb(31/255,119/255,180/255,1), lwd=2, lty=2)
curve(param_results_raw[5,1] + param_results_raw[11,1]*x + param_results_raw[17,1]*x^2, 
      from = min(train_data_raw$sum_px_s), to =max(train_data_raw$sum_px_s),add = T,
      col =  rgb(255/255,127/255,14/255,1), lwd=2, lty=2)
curve(param_results_raw[6,1] + param_results_raw[12,1]*x + param_results_raw[18,1]*x^2, 
      from = min(train_data_raw$sum_px_s), to =max(train_data_raw$sum_px_s),add = T,
      col =  rgb(44/255,160/255,44/255,1), lwd=2, lty=2)

param_results_unet <- data.frame(jagsresults(jags.parsamples_unet, c('beta0', 'beta1','beta2')))

plot(train_data_unet$sum_px_s[train_data_unet$group==1], train_data_unet$count[train_data_unet$group==1],
     ylim = c(0, 100), xlim = c(min(train_data_unet$sum_px_s),max(train_data_unet$sum_px_s)), xaxt='n',
     pch = 16, col = rgb(31/255,119/255,180/255,0.2), frame.plot=F,
     xlab = 'Sum of Pixels', ylab = 'Cell Count')
at =seq(min(train_data_unet$sum_px_s),max(train_data_unet$sum_px_s), by=0.5)
labels <- at*sd_sum_px + mean_sum_px
axis(side=1, at = at, labels = round(labels, 0))
points(train_data_unet$sum_px_s[train_data_unet$group==2], train_data_unet$count[train_data_unet$group==2],
       pch = 16, col =  rgb(255/255,127/255,14/255,0.2))
points(train_data_unet$sum_px_s[train_data_unet$group==3], train_data_unet$count[train_data_unet$group==3],
       pch = 16, col =  rgb(44/255,160/255,44/255,0.2))

points(train_data_unet$sum_px_s[train_data_unet$group==4], train_data_unet$count[train_data_unet$group==4],
       pch=17, col = rgb(31/255,119/255,180/255,0.2))
points(train_data_unet$sum_px_s[train_data_unet$group==5], train_data_unet$count[train_data_unet$group==5],
       pch = 17, col =  rgb(255/255,127/255,14/255,0.2))
points(train_data_unet$sum_px_s[train_data_unet$group==6], train_data_unet$count[train_data_unet$group==6],
       pch = 17, col =  rgb(44/255,160/255,44/255,0.2))

curve(param_results_unet[1,1] + param_results_unet[7,1]*x + param_results_unet[13,1]*x^2, 
      from = min(train_data_unet$sum_px_s), to =max(train_data_unet$sum_px_s),add = T,
      col =  rgb(31/255,119/255,180/255,1), lwd=2)
curve(param_results_unet[2,1] + param_results_unet[8,1]*x + param_results_unet[14,1]*x^2, 
      from = min(train_data_unet$sum_px_s), to =max(train_data_unet$sum_px_s),add = T,
      col =  rgb(255/255,127/255,14/255,1), lwd=2)
curve(param_results_unet[3,1] + param_results_unet[9,1]*x + param_results_unet[15,1]*x^2, 
      from = min(train_data_unet$sum_px_s), to =max(train_data_unet$sum_px_s),add = T,
      col =  rgb(44/255,160/255,44/255,1), lwd=2)

curve(param_results_unet[4,1] + param_results_unet[10,1]*x + param_results_unet[16,1]*x^2, 
      from = min(train_data_unet$sum_px_s), to =max(train_data_unet$sum_px_s),add = T,
      col =  rgb(31/255,119/255,180/255,1), lwd=2, lty=2)
curve(param_results_unet[5,1] + param_results_unet[11,1]*x + param_results_unet[17,1]*x^2, 
      from = min(train_data_unet$sum_px_s), to =max(train_data_unet$sum_px_s),add = T,
      col =  rgb(255/255,127/255,14/255,1), lwd=2, lty=2)
curve(param_results_unet[6,1] + param_results_unet[12,1]*x + param_results_unet[18,1]*x^2, 
      from = min(train_data_unet$sum_px_s), to =max(train_data_unet$sum_px_s),add = T,
      col =  rgb(44/255,160/255,44/255,1), lwd=2, lty=2)







param_results_unet <- data.frame(jagsresults(jags.parsamples_unet, c('beta0', 'beta1','beta2')))
ggplot(data = train_data_unet, aes(x = sum_px_s, y = count, color = as.factor(blur), size = as.factor(stain))) + 
  geom_point() +
  stat_function(fun = function(x) param_results_unet[1,1]+param_results_unet[6+1,1]*x + param_results_unet[12+1,1]*x^2,size = 1, color = 'red') +
  stat_function(fun = function(x)param_results_unet[2,1]+param_results_unet[6+2,1]*x + param_results_unet[12+2,1]*x^2, size = 1, color = 'green') +
  stat_function(fun = function(x) param_results_unet[3,1]+param_results_unet[6+3,1]*x + param_results_unet[12+3,1]*x^2, size = 1, color = 'blue') +
  stat_function(fun = function(x) param_results_unet[4,1]+param_results_unet[6+4,1]*x + param_results_unet[12+4,1]*x^2,size = 2, color = 'red') +
  stat_function(fun = function(x)param_results_unet[5,1]+param_results_unet[6+5,1]*x + param_results_unet[12+5,1]*x^2, size = 2, color = 'green') +
  stat_function(fun = function(x) param_results_unet[6,1]+param_results_unet[6+6,1]*x + param_results_unet[12+6,1]*x^2, size = 2, color = 'blue') +
  ylim(0,101)

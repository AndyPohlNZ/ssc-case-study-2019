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

train_data = read.table('./ssc-data/CompTrain/labels.csv', sep = ',', header = T, stringsAsFactors = F)
train_data$sum_px = rep(0, nrow(train_data))
train_data$sum_px_s = rep(0, nrow(train_data))

NTRAIN = nrow(train_data)
for(i in 1:NTRAIN){
  img_file = train_data[i,1]
  img_file = gsub('.TIF', '.png', img_file, ignore.case = T)
  img = image_read(paste0('./ssc-data/CompTrain/masks/', img_file))
  #img = image_scale(img, "128x128!")
  img = as.integer(img[[1]])/255
  #img = as.matrix(raster(as.character(paste0('./ssc-data/CompTrain/', img_file)))[,,1])/255
  train_data$sum_px[i] = sum(img)
  
}
train_data$stain = train_data$stain-1
train_data= train_data[1:NTRAIN, ]

mean_sum_px = mean(train_data$sum_px)
sd_sum_px = sd(train_data$sum_px)
train_data$sum_px_s <- sapply(train_data$sum_px,function(x){z_standardise(x, mean_sum_px, sd_sum_px)})


test_data = read.table('./ssc-data/CompValid/labels.csv', sep = ',', header = T, stringsAsFactors = F)
test_data$sum_px = rep(0, nrow(test_data))
test_data$sum_px_s = rep(0, nrow(test_data))

NTEST = nrow(test_data)
for(i in 1:NTEST){
  
  img_file = test_data[i,1]
  img_file = gsub('.TIF', '.png', img_file, ignore.case = T)
  img = image_read(paste0('./ssc-data/CompValid/masks/', img_file))
  #img = image_scale(img, "128x128!")
  img = as.integer(img[[1]])/255
  #img = as.matrix(raster(as.character(paste0('./ssc-data/CompTrain/', img_file)))[,,1])/255
  test_data$sum_px[i] = sum(img)
}
test_data$stain = test_data$stain-1
test_data= test_data[1:NTEST, ]

test_data$sum_px_s <- sapply(test_data$sum_px,function(x){z_standardise(x, mean_sum_px, sd_sum_px)})


#################################################################
### ### Cells = circles 
#################################################################
jinits <- function() {
  ### all the other params ###
  .Rng.name="lecuyer::RngStream"
  .Rng.seed=randomNumbers(n = 1, min = 1, max = 1e+06,col=1)
}
jags_input = list(N = nrow(train_data), W = max(dim(img)), y = train_data$count, sum_px = train_data$sum_px,
                  NTest = nrow(test_data), sum_px_test = test_data$sum_px)
jags_file_name = './PoisReg/Circle_model/simpleCircle.jags'
params_to_monitor =c('y_hat', 'lambda_hat', 'gamma')
jags.parsamples <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
  model.jags <- jags.model(file = jags_file_name, data =jags_input,
                           n.chain=1, n.adapt=1000)
  result <- coda.samples(model = model.jags, variable.names = params_to_monitor,n.iter = 8000)
  return(result)
}

samps_m = as.matrix(jags.parsamples)
param = "gamma"
hist(samps_m[,colnames(samps_m) == param ], breaks = 100, xlab = param, main = param, freq = F)
mcmc_trace(jags.parsamples, pars=c('lambda_hat[476]'))

results = data.frame(jagsresults(jags.parsamples, 'y_hat'))
RMSE(x = test_data$count, results$mean)

test_data$yhat = results$mean

#################################################################
### Cells = circles adjusting for different radius between stain1 and 2
#################################################################
jinits <- function() {
  ### all the other params ###
  .Rng.name="lecuyer::RngStream"
  .Rng.seed=randomNumbers(n = 1, min = 1, max = 1e+06,col=1)
}
jags_input = list(N = nrow(train_data), W = max(dim(img)), y = train_data$count, sum_px = train_data$sum_px, stain = train_data$stain,
                  NTest = nrow(test_data), sum_px_test = test_data$sum_px, stain_test = test_data$stain)
jags_file_name = './PoisReg/Circle_model/simpleCircle2.jags'
params_to_monitor =c('y_hat', 'lambda_hat', 'gamma', 'gamma_s')
jags.parsamples <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
  model.jags <- jags.model(file = jags_file_name, data =jags_input,
                           n.chain=1, n.adapt=1000)
  result <- coda.samples(model = model.jags, variable.names = params_to_monitor,n.iter = 8000)
  return(result)
}

samps_m = as.matrix(jags.parsamples)
param = "lambda_hat[479]"
hist(samps_m[,colnames(samps_m) == param ], breaks = 100, xlab = param, main = param, freq = F)
mcmc_trace(jags.parsamples, pars=param)

results = data.frame(jagsresults(jags.parsamples, 'y_hat'))
RMSE(x = test_data$count, results$mean)

test_data$yhat = results$mean

#################################################################
### Cells = circles adjusting for overlap
#################################################################
jinits <- function() {
  ### all the other params ###
  .Rng.name="lecuyer::RngStream"
  .Rng.seed=randomNumbers(n = 1, min = 1, max = 1e+06,col=1)
}
jags_input = list(N = nrow(train_data), W = max(dim(img)), y = train_data$count, sum_px = train_data$sum_px, stain = train_data$stain,
                  NTest = nrow(test_data), sum_px_test = test_data$sum_px, stain_test = test_data$stain)
jags_file_name = './PoisReg/Circle_model/simpleCircle3.jags'
params_to_monitor =c('y_hat', 'lambda_hat', 'gamma', 'gamma_s', 'beta0', 'beta1')
jags.parsamples <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
  model.jags <- jags.model(file = jags_file_name, data =jags_input,
                           n.chain=1, n.adapt=1000)
  result <- coda.samples(model = model.jags, variable.names = params_to_monitor,n.iter = 8000)
  return(result)
}

samps_m = as.matrix(jags.parsamples)
param = "beta0"
hist(samps_m[,colnames(samps_m) == param ], breaks = 100, xlab = param, main = param, freq = F)
mcmc_trace(jags.parsamples, pars=param)

results = data.frame(jagsresults(jags.parsamples, 'y_hat'))
RMSE(x = test_data$count, results$mean)

test_data$yhat = results$mean

#################################################################
### Cells = circles adjusting for overlap deterministically
#################################################################
jinits <- function() {
  ### all the other params ###
  .Rng.name="lecuyer::RngStream"
  .Rng.seed=randomNumbers(n = 1, min = 1, max = 1e+06,col=1)
}
jags_input = list(N = nrow(train_data), W = max(dim(img)), y = train_data$count, sum_px = train_data$sum_px, stain = train_data$stain,
                  NTest = nrow(test_data), sum_px_test = test_data$sum_px, stain_test = test_data$stain)
jags_file_name = './PoisReg/Circle_model/simpleCircle4.jags'
params_to_monitor =c('y_hat', 'lambda_hat', 'gamma', 'gamma_s', 'beta', 'theta')
jags.parsamples <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
  model.jags <- jags.model(file = jags_file_name, data =jags_input,
                           n.chain=1, n.adapt=1000)
  result <- coda.samples(model = model.jags, variable.names = params_to_monitor,n.iter = 8000)
  return(result)
}

samps_m = as.matrix(jags.parsamples)
param = "theta"
hist(samps_m[,colnames(samps_m) == param ], breaks = 100, xlab = param, main = param, freq = F)
mcmc_trace(jags.parsamples, pars=param)

results = data.frame(jagsresults(jags.parsamples, 'y_hat'))
RMSE(x = test_data$count, results$mean)

test_data$yhat = results$mean

#################################################################
### Cells = circles adjusting for overlap 
#################################################################
jinits <- function() {
  ### all the other params ###
  .Rng.name="lecuyer::RngStream"
  .Rng.seed=randomNumbers(n = 1, min = 1, max = 1e+06,col=1)
}
jags_input = list(N = nrow(train_data), W = max(dim(img)), y = train_data$count, sum_px = train_data$sum_px, stain = train_data$stain,
                  NTest = nrow(test_data), sum_px_test = test_data$sum_px, stain_test = test_data$stain)
jags_file_name = './PoisReg/Circle_model/simpleCircle5.jags'
params_to_monitor =c('y_hat', 'lambda_hat', 'gamma', 'gamma_s', 'beta')
jags.parsamples <- foreach(i=1:getDoParWorkers(),  .combine='c', .final=mcmc.list) %dopar% {
  model.jags <- jags.model(file = jags_file_name, data =jags_input,
                           n.chain=1, n.adapt=1000)
  result <- coda.samples(model = model.jags, variable.names = params_to_monitor,n.iter = 1000)
  return(result)
}

samps_m = as.matrix(jags.parsamples)
param = "beta"
hist(samps_m[,colnames(samps_m) == param ], breaks = 100, xlab = param, main = param, freq = F)
mcmc_trace(jags.parsamples, pars=param)

results = data.frame(jagsresults(jags.parsamples, 'y_hat'))
RMSE(x = test_data$count, results$mean)

test_data$yhat = results$mean



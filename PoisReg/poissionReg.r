library(rjags)
library(coda)
library(png)

RMSE <- function(x, xhat){
  return(sqrt(mean((x-xhat)^2)))
}
z_standardise <- function(x,xbar, s){
  return((x-xbar)/s)
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
TRAIN_DIR = paste0(HOME_DIR, 'ssc-data/CompTrain/')

setwd(HOME_DIR)

train_data = read.table('./ssc-data/CompTrain/labels.csv', sep = ',', header = T, stringsAsFactors = F)
train_data$mask_pixels = rep(0, nrow(train_data))
for(i in 1:nrow(train_data)){
  img_file = train_data[i,1]
  mask_file = gsub('.TIF', '.png',img_file, ignore.case = T)
  mask = readPNG(paste0('./ssc-data/CompTrain/masks/', mask_file))
  train_data$mask_pixels[i] = sum(mask)/(256*256)
}
maskpx_mean = mean(train_data$mask_pixels)
maskpx_sd = sd(train_data$mask_pixels)

train_data$mask_pixels_s <- sapply(train_data$mask_pixels,function(x){z_standardise(x, maskpx_mean, maskpx_sd)})
train_data$stain <- train_data$stain-1


test_data = read.table('./ssc-data/CompValid/labels.csv', sep = ',', header = T, stringsAsFactors = F)
test_data$mask_pixels = rep(0, nrow(train_data))
for(i in 1:nrow(test_data)){
  img_file = test_data[i,1]
  mask_file = gsub('.TIF', '.png',img_file, ignore.case = T)
  mask = readPNG(paste0('./ssc-data/CompValid/masks/', mask_file))
  test_data$mask_pixels[i] = sum(mask)/(256*256)
}

test_data$mask_pixels_s <- sapply(test_data$mask_pixels,function(x){z_standardise(x, maskpx_mean, maskpx_sd)})
test_data$stain <- test_data$stain-1


#################################################################
### Simple pixel count regression...####
#################################################################
jags_input = list(N = nrow(train_data), x = train_data$mask_pixels_s, y = train_data$count,
                  NTest = nrow(test_data), xtest = test_data$mask_pixels_s)
jags_file_name = './PoisReg/pxCnt.jags'
jagsmdl <- jags.model(file = jags_file_name, data = jags_input, 
                      n.chains = 4, n.adapt = 2000) # does n.adapt equate to burnin...?
samps <- coda.samples(jagsmdl, variable.names = c('beta0', 'beta1','y_hat', 'lambda_hat'), n.iter = 8000, thin = 1)
samps_m = as.matrix(samps)
colnames(samps_m)
param = "beta0"
hist(samps_m[,colnames(samps_m) == param ], breaks = 100, xlab = param, main = param, freq = F)

results = data.frame(jagsresults(samps, 'y_hat'))

RMSE(x = test_data$count, results$mean)


#################################################################
### Simple pixel count regression via 0...####
#################################################################
jags_input = list(N = nrow(train_data), x = train_data$mask_pixels_s, x2 = train_data$stain, y = train_data$count,
                  NTest = nrow(test_data), xtest = test_data$mask_pixels_s, xtest2 = test_data$stain)
jags_file_name = './PoisReg/pxCnt_zerointcpt.jags'
jagsmdl <- jags.model(file = jags_file_name, data = jags_input, 
                      n.chains = 4, n.adapt = 2000) # does n.adapt equate to burnin...?
samps <- coda.samples(jagsmdl, variable.names = c('beta1','beta2','y_hat', 'lambda_hat'), n.iter = 8000, thin = 1)
samps_m = as.matrix(samps)
colnames(samps_m)
param = "y_hat[142]"
hist(samps_m[,colnames(samps_m) == param ], breaks = 100, xlab = param, main = param, freq = F)

results = data.frame(jagsresults(samps, 'y_hat'))

RMSE(x = test_data$count, results$mean)

test_data$image_name[142]
#################################################################
### Simple pixel count  + stain regression...####
#################################################################
jags_input = list(N = nrow(train_data), x1 = train_data$mask_pixels_s, x2 = train_data$stain, y = train_data$count,
                  NTest = nrow(test_data), xtest1 = test_data$mask_pixels_s, xtest2 = test_data$stain)
jags_file_name = './PoisReg/pxCnt_Stain.jags'
jagsmdl <- jags.model(file = jags_file_name, data = jags_input, 
                      n.chains = 4, n.adapt = 2000) # does n.adapt equate to burnin...?
samps <- coda.samples(jagsmdl, variable.names = c('beta0', 'beta1','beta2', 'y_hat', 'lambda_hat'), n.iter = 8000, thin = 1)
samps_m = as.matrix(samps)
colnames(samps_m)
param = "beta2"
hist(samps_m[,colnames(samps_m) == param ], breaks = 100, xlab = param, main = param, freq = F)

results = data.frame(jagsresults(samps, 'y_hat'))
RMSE(x = test_data$count, results$mean)

#################################################################
### Simple pixel count*stain regression...####
#################################################################
jags_input = list(N = nrow(train_data), x1 = train_data$mask_pixels_s, x2 = train_data$stain, y = train_data$count,
                  NTest = nrow(test_data), xtest1 = test_data$mask_pixels_s, xtest2 = test_data$stain)
jags_file_name = './PoisReg/full_interaction.jags'
jagsmdl <- jags.model(file = jags_file_name, data = jags_input, 
                      n.chains = 4, n.adapt = 2000) # does n.adapt equate to burnin...?
samps <- coda.samples(jagsmdl, variable.names = c('beta0', 'beta1','beta2', 'beta3', 'y_hat', 'lambda_hat'), n.iter = 8000, thin = 1)
samps_m = as.matrix(samps)
colnames(samps_m)
param = "beta3"
hist(samps_m[,colnames(samps_m) == param ], breaks = 100, xlab = param, main = param, freq = F)

results = data.frame(jagsresults(samps, 'y_hat'))
RMSE(x = test_data$count, results$mean)

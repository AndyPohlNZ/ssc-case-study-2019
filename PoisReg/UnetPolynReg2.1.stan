//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.

data {
  int<lower=0> N;
  //int<lower=0> W;
  int<lower=0> NTest;
  int<lower = 0>NGroup;
  real<lower=0, upper=100>  y[N];
  int<lower=0, upper=NGroup> g[N];
  real<lower = 0>[N]  x;
  real<lower = 0>[NTest]  xtest;
  int<lower=0, upper=NGroup> gtest[NTest];

}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[NGroup] beta0;
  vector[NGroup] beta1;
  vector[NGroup] beta2;
  real<lower =0> sigma;
}



// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for(i in 1:N){
    y[i] ~ normal();

  }
    //prior
  gamma_s ~ normal(20, 10);
  gamma ~ normal(30, 10);
  theta ~ uniform(0.7, 1);
  beta ~ uniform(0,0.05);
}

  generated quantities {
    int y_hat[NTest] ;
    vector[NTest] lambda_hat;
    vector[NTest] llambda_hat;
    vector[NTest] gamma_adj_hat;
    for(i in 1:NTest){
      y_hat[i] = poisson_rng(llambda_hat[i]);
      llambda_hat[i] = log(lambda_hat[i]);
       gamma_adj_hat[i] = gamma + gamma_s*stain_test[i];

      lambda_hat[i] = sum_px[i]/ ((theta*pi() * gamma_adj_hat[i]*gamma_adj_hat[i]) + (1-theta)*(pi()*gamma_adj_hat[i]*gamma_adj_hat[i] - (2*beta*pi()*gamma_adj_hat[i]*gamma_adj_hat[i])));
    }
  }


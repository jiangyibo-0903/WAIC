functions {  
  real dGORHPieces_lpdf(real x, real gam, real reg, vector lambda, int v, vector S) {  
    real logProb;  
    if (v == 1) {  
      logProb = reg - (1 + 1 / gam) * log(1 + gam * x * exp(reg));  
    } else {  
      real M = lambda[v] * (x - S[v-1]);  
      for (g in 1:(v-1)) {  
        real K;  
        if (g == 1) {  
          K = lambda[1] * S[1];  
        } else {  
          K = lambda[g] * (S[g] - S[g-1]);  
        }  
        M = M + K;  
      }  
      logProb = log(lambda[v]) + reg - (1 + 1 / gam) * log(1 + gam * exp(reg) * M);  
    }  
    return logProb;  
  }  
}  
  
data {  
  int<lower=0> N;   
  int<lower=0> J; 
  int<lower=0> V;   
  int<lower=0, upper=1> resp[N, J];  
  matrix[N, J] rt; 
  vector[V] s; 
}  
  
transformed data {
  int<lower=0> OR[N,J];
  {  
    for (i in 1:N) {  
      for (j in 1:J) {  
        int v = 1;  
        while (v < V && rt[i, j] > s[v]) {  
          v = v + 1;  
        }  
        OR[i, j] = v;  
      }  
    }  
  }  
}  
  
parameters {  
  vector[N] theta;  
  vector[N] tau;  
  vector[J] a;  
  vector[J] b;  
  vector<lower=0>[J] gamma;  
  vector[J] phi;  
  vector[J] zeta;  
  real<lower=-pi()/2, upper=pi()/2> varphi;  
  real mu_b;  
  real<lower=0> sigma2_b;  
  vector<lower=0>[V-1] lambda_free; 
}  
transformed parameters {  
  vector<lower=0>[V] lambda;    
  lambda[1] = 1; 
  for (g in 2:V) {  
    lambda[g] = lambda_free[g-1];   
  }  
}    
model {  
  for (j in 1:J) {  
    a[j] ~ lognormal(0, 1);  
    phi[j] ~ lognormal(0, 1);  
    gamma[j] ~ inv_gamma(1, 1);  
    b[j] ~ normal(mu_b, sqrt(sigma2_b));  
    zeta[j] ~ normal(0, 10);  
  }
  for (i in 1:N) {  
    theta[i] ~ normal(0, 1);  
    tau[i] ~ normal(0, 1);  
  }
  mu_b ~ normal(0, sqrt(10 * sigma2_b));  
  sigma2_b ~ inv_gamma(0.01, 0.01);
  for (g in 1:(V-1)) {  
    lambda_free[g] ~ inv_gamma(0.01, 0.01);  
  }
  for (i in 1:N) {  
    for (j in 1:J) {  
      real ind = a[j] * (theta[i] - b[j]);  
      real p = 1 / (1 + exp(-ind)) * (ind > 0) + exp(ind) / (1 + exp(ind)) * (ind <= 0);  
      resp[i, j] ~ bernoulli(p);
      target += dGORHPieces_lpdf(rt[i, j]| gamma[j], phi[j] * (theta[i] * sin(varphi) + tau[i] * cos(varphi) - zeta[j]),lambda, OR[i, j], s);  
    }  
  }  
}

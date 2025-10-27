library(readr)
library(rstan)
library(loo)
library(dplyr)  
library(openxlsx)
DETA_ijv <- function(T, S){
  n_student <- nrow(T)
  n_item <- ncol(T)
  V <- length(S)
  DETA <- array(0, dim = c(n_student, n_item, V))
  order_interval <- matrix(, n_student, n_item)
  for(i in 1:n_student){
    for(j in 1:n_item){
      v <- which(T[i, j] - S <= 0)[1]
      order_interval[i, j] <- v
      DETA[i, j, v] <- 1
    }}
  return(list(order = order_interval, DETA = DETA))
}

res <- read_csv("D:\\新建文件夹\\数据\\反应数据\\res.csv")
rt <- read_csv("D:\\新建文件夹\\数据\\反应时间数据\\rt.csv")
rt<-as.matrix(rt[2:length(rt)])
res<-as.matrix(res[2:length(res)])
dim(rt)
dim(res)
size = dim(rt)
N = size[1]
J = size[2]
V=5
set.seed(V)
print(getwd())
print(paste0('V=', V))
S <- quantile(rt, probs=seq(0, 1, length = V + 1))[-c(1, V + 1)]

length(S)
S
s<-c(S,max(rt)+0.1)
s
write.csv(s,"s.csv")
write.csv(S,"SS.csv")
OR<-DETA_ijv(rt,s)$order
write.csv(OR,"OR_Pieces.csv")
data <- list(N=N,J=J,V=V,rt = rt, resp = res, s=s)
data
init <- list(list(theta = rnorm(N),   
                  a = runif(J, 0.5, 1.5),  
                  b = rnorm(J),  
                  mu_b = 0,   
                  sigma2_b = 1,  
                  phi = runif(J, 0.5, 1.5),   
                  zeta = rnorm(J),  
                  tau = rnorm(N),   
                  varphi = 0.5,  
                  gamma = rep(1, J),   
                  lambda = rep(1, V)))

init
fir = Sys.time()
fir
mcmc.out<-stan(file="m.stan",data = data,chains = 1,iter = 6000,init = init,warmup = 2000)
end = Sys.time()
end - fir
samples <- extract(mcmc.out) 
write.csv(samples,paste0('samples','-',N,'-',J,'.csv'))
par.summary<-summary(mcmc.out)
all.params<-as.data.frame(par.summary)
write.csv(all.params,paste0('allparams','-',N,'-',J,'.csv'))
samples <- read_csv(paste0('samples','-',N,'-',J,'.csv'))
mcmc_theta <- samples[, 1:N+1] 
mcmc_tau <- samples[, (N+2):(2*N+1)]
mcmc_a <- samples[, (2*N+2):(2*N+J+1)]  
mcmc_b <- samples[, (2*N+J+2):(2*N+2*J+1)]  
mcmc_gamma <- samples[, (2*N+2*J+2):(2*N+3*J+1)] 
mcmc_phi <- samples[, (3*J+2*N+2):(2*N+4*J+1)]
mcmc_zeta <- samples[, (2*N+4*J+2):(5*J+2*N+1)]
mcmc_varphi<-samples[,2*N+5*J+2]
mcmc_lambda <- samples[, (ncol(samples) - V):(ncol(samples)-1)]
#条件抽样

eap_a <- colMeans(mcmc_a)  
eap_b <- colMeans(mcmc_b)  
eap_th <- colMeans(mcmc_theta)  
eap_gamma <- colMeans(mcmc_gamma)  
eap_lambda <- colMeans(mcmc_lambda)  
eap_phi <- colMeans(mcmc_phi)  
eap_tau <- colMeans(mcmc_tau)  
eap_zeta <- colMeans(mcmc_zeta)
eap_varphi <- colMeans(mcmc_varphi)
param_names <- c(
  paste0("a_", 1:J),
  paste0("b_", 1:J),
  paste0("theta_", 1:N),
  paste0("gamma_", 1:J),
  paste0("lambda_", 1:V),
  paste0("phi_", 1:J),
  paste0("tau_", 1:N),
  paste0("zeta_", 1:J),
  "varphi"
)

param_values <- c(
  eap_a,
  eap_b,
  eap_th,
  eap_gamma,
  eap_lambda,
  eap_phi,
  eap_tau,
  eap_zeta,
  eap_varphi
)

# 创建数据框
eap_df <- data.frame(
  parameter = param_names,
  value = param_values
)

# 保存为单个CSV文件
write.csv(eap_df, 
          paste0('all_eap_results-', N, '-', J, '.csv'), 
          row.names = FALSE)
#给定rt
data1 <- list(N=N,J=J,V=V,rt = rt, resp = res, s=s,
              constant_gamma=eap_gamma,constant_phi=eap_phi,
              constant_zeta=eap_zeta,constant_varphi=eap_varphi,constant_lambda=eap_lambda)
init1 <- list(list(theta = rnorm(N),   
                   a = runif(J, 0.5, 1.5),  
                   b = rnorm(J),  
                   mu_b = 0,   
                   sigma2_b = 1,  
                   tau = rnorm(N)))
mcmc.out1<-stan(file="m_resp_give_rt.stan",data = data1,chains = 1,iter = 4000,init = init1,warmup = 2000)
samples1<-extract(mcmc.out1)
traceplot(mcmc.out1,pars=c("b","a"), inc_warmup = FALSE)
write.csv(samples1,"samples_res_give_rt.csv")
par.summary1<-summary(mcmc.out1)
all.params2<-as.data.frame(par.summary1)
write.csv(all.params2,"allparams_res_give_rt.csv")

#给定resp
data2<-list(N=N,J=J,V=V,rt = rt, resp = res, s=s,constant_a=eap_a,constant_b=eap_b)
init2<-list(list(theta = rnorm(N),   
                 mu_b = 0,   
                 sigma2_b = 1,  
                 phi = runif(J, 0.5, 1.5),   
                 zeta = rnorm(J),  
                 tau = rnorm(N),   
                 varphi = 0.5,  
                 gamma = rep(1, J),   
                 lambda = rep(1, V)))
mcmc.out2<-stan(file="m_rt_give_resp.stan",data = data2,chains = 1,iter = 8000,init = init2,warmup = 2000)

samples2<-extract(mcmc.out2)

traceplot(mcmc.out2,pars=c("gamma"), inc_warmup = FALSE)
write.csv(samples2,"samples_rt_give_resp.csv")
par.summary2<-summary(mcmc.out2)
all.params3<-as.data.frame(par.summary2)
write.csv(all.params3,"allparams_rt_give_resp.csv")
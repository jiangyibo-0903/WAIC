clear
clc

Model=Model1();
V=65;
resp=csvread(['resp.csv'],1,1);
rt=csvread(['rt.csv'],1,1);
jk=size(resp);
N=jk(1);
J=jk(2); 
S=csvread(['s.csv'],1,1);
OR=csvread(['OR_Pieces.csv'],1,1);

samples=csvread(['samples.csv'],1,1); 
dim=size(samples);
niter=dim(1);

mcmc_theta=samples(:,1:N);
mcmc_tau=samples(:,(N+1):(2*N));
mcmc_a=samples(:,(2*N+1):(2*N+J));
mcmc_b=samples(:,(2*N+J+1):(2*N+2*J));
mcmc_gamma=samples(:,(2*N+2*J+1):(2*N+3*J));
mcmc_phi=samples(:,(2*N+3*J+1):(2*N+4*J));
mcmc_zeta=samples(:,(2*N+4*J+1):(2*N+5*J));
mcmc_varphi=samples(:,2*N+5*J+1);
mcmc_lambda=samples(:,(dim(2)-65:(dim(2)-1)));

eap_a=mean(mcmc_a,1);
eap_b=mean(mcmc_b,1);
eap_gamma=mean(mcmc_gamma,1);
eap_lambda=mean(mcmc_lambda,1);
eap_phi=mean(mcmc_phi,1);
eap_th=mean(mcmc_theta,1);
eap_tau=mean(mcmc_tau,1);
eap_varphi=mean(mcmc_varphi);
eap_zeta=mean(mcmc_zeta,1);

%P的有关计算
for i=1:N
   f1=normpdf(eap_th(i),0,1);
   y=resp(i,:);
   kj=(eap_a*eap_th(i)-eap_a.*eap_b).*y-log(1+exp(eap_a*eap_th(i)-eap_a.*eap_b));
   Ln_hat(i)=sum(kj)+log(f1);
end
logf_p_resp=zeros(niter,N);
parfor r=1:(niter)
    a=mcmc_a(r,:);
    b=mcmc_b(r,:);   
    k1=zeros(1,N);
    for i=1:N
       y=resp(i,:);
       fun = @ (x) exp(log(normpdf(x,0,1))+sum((a*x-a.*b).*y-log(1+exp(a*x-a.*b)))-Ln_hat(i));
       h=integral(fun,-Inf,Inf,'ArrayValued',true);
       k1(i)=log(h)+Ln_hat(i);
    end
    logf_p_resp(r,:)=k1;
end
%rt
HH=zeros(V,1);
HH(1)=S(1);
for m=2:(V-1)
    HH(m)=HH(m-1)+eap_lambda(m)*(S(m)-S(m-1));
end
Ln_hat=zeros(N,1);
for i=1:N
    taui_star=eap_th(i)*sin(eap_varphi)+eap_tau(i)*cos(eap_varphi);
    f1=normpdf(taui_star,0,1);
    t=rt(i,:);
    ORi=OR(i,:);
    ein1=eap_phi.*(taui_star-eap_zeta);
    ht=zeros(1,J);
    for j=1:J
        k=ORi(j);
        if k==1
            ht(j)=t(j);
        else
            ht(j)=eap_lambda(k)*(t(j)-S(k-1))+HH(k-1);
        end
    end
    ein2=(1+1./eap_gamma).*log(1+eap_gamma.*exp(ein1).*ht);
    Ln_hat(i)=sum(ein1-ein2+log(eap_lambda(ORi)))+log(f1);
end

logf_p_rt=zeros(niter,N);

parfor r=1:(niter)
    
    gamma=mcmc_gamma(r,:);
    lambda=mcmc_lambda(r,:);
    phi=mcmc_phi(r,:);
    zeta=mcmc_zeta(r,:);
    HH=zeros(V,1);
    HH(1)=S(1);
    for m=2:(V-1)
        HH(m)=HH(m-1)+lambda(m)*(S(m)-S(m-1));
    end
    kl=zeros(N,1);
    for i=1:N
        t=rt(i,:);
        ORi=OR(i,:);
        ht=zeros(1,J);
        for j=1:J
            k=ORi(j);
            if k==1
                ht(j)=t(j);
            else
                ht(j)=lambda(k)*(t(j)-S(k-1))+HH(k-1);
            end
        end
        fun = @ (x)  exp(log(normpdf(x,0,1))+sum(phi.*(x-zeta)+log(lambda(ORi))-(1+1./gamma).*log(1+gamma.*exp(phi.*(x-zeta)).*ht)));
        h=integral(fun,-Inf,Inf,'ArrayValued',true);
        kl(i)=log(h);
        
    end
    logf_p_rt(r,:)=kl;
end
%AC|rt
samples1=csvread('samples_res_give_rt.csv',1,1); 
dim=size(samples1);
niter=dim(1);
dim=size(samples1);
niter=dim(1);
mcmc_a=samples1(:,(2*N+1):(2*N+J));
mcmc_b=samples1(:,(2*N+J+1):(2*N+2*J));
mcmc_theta=samples1(:,1:N);
mcmc_tau=samples1(:,N+1:2*N);
eap_a=mean(mcmc_a,1);
eap_b=mean(mcmc_b,1);
eap_th=mean(mcmc_theta,1);
eap_tau=mean(mcmc_tau,1);

CPOR=zeros(niter,N);
logf_res_give_rt=zeros(niter,N);
rng(2)
parfor r=1:niter
    a=mcmc_a(r,:);
    b=mcmc_b(r,:);
    tau=mcmc_tau(r,:);
    th=mcmc_theta(r,:);
    fk=zeros(N,1);
    for i=1:N
        res=resp(i,:);
        t=rt(i,:);
        x0=[th(i),tau(i)];
        [mode,fval]=fminunc(@(x) NLogf(x,t,res,a,b,phi,zeta,gamma,varphi,lambda,S,OR(i,:)),x0);
        sigma=Var_Cov(mode,t,a,b,phi,zeta,gamma,varphi,lambda,S,OR(i,:));
        %R=[th(i),tau(i)];
        R=mvnrnd(mode,sigma,1000);
        kb=mvnpdf(R,mode,sigma);
        fb=JointpdfMatrix(R,t,res,a,b,phi,zeta,gamma,varphi,lambda,S,OR(i,:));
        %CPOR(r,i)=mean(kb./fb);
        fk(i)=mean(fb./kb);
    end
    logf_res_give_rt(r,:)=log(fk./fm);
    
end


%RT|res

samples2=csvread(['samples_rt_give_resp.csv'],1,1);
dim=size(samples2);
niter=dim(1);
mcmc_theta=samples2(:,1:N);
mcmc_tau = samples2(:, N+1:2*N);  
mcmc_gamma = samples2(:, (2*N+1):(2*N+J));  
mcmc_phi = samples2(:, (2*N+J+1):(2*N+2*J));  
mcmc_zeta = samples2(:, (2*N+2*J+1):(3*J+2*N));
mcmc_lambda=samples2(:,(dim(2)-65:(dim(2)-1)));
mcmc_varphi=samples2(:,2*N+3*J+1);
eap_gamma=mean(mcmc_gamma,1);
eap_lambda=mean(mcmc_lambda,1);
eap_phi=mean(mcmc_phi,1);
eap_theta=mean(mcmc_theta,1);
eap_tau=mean(mcmc_tau,1);
eap_zeta=mean(mcmc_zeta,1);
eap_varphi=mean(mcmc_varphi,1);
HH=zeros(V,1);
HH(1)=S(1);
for m=2:(V-1)
    HH(m)=HH(m-1)+eap_lambda(m)*(S(m)-S(m-1));
end
Ln_hat=zeros(N,1);
for i=1:N
    taui_star=eap_th(i)*sin(eap_varphi)+eap_tau(i)*cos(eap_varphi);
    f1=normpdf(taui_star,0,1);
    t=rt(i,:);
    ORi=OR(i,:);
    ein1=eap_phi.*(taui_star-eap_zeta);
    ht=zeros(1,J);
    for j=1:J
        k=ORi(j);
        if k==1
            ht(j)=t(j);
        else
            ht(j)=eap_lambda(k)*(t(j)-S(k-1))+HH(k-1);
        end
    end
    ein2=(1+1./eap_gamma).*log(1+eap_gamma.*exp(ein1).*ht);
    Ln_hat(i)=sum(ein1-ein2+log(eap_lambda(ORi)))+log(f1);
end
logf_rt_give_res=zeros(niter,N);
parfor r=1:niter
    
    gamma=mcmc_gamma(r,:);
    lambda=mcmc_lambda(r,:);
    phi=mcmc_phi(r,:);
    zeta=mcmc_zeta(r,:);
    HH=zeros(V,1);
    HH(1)=S(1);
    for m=2:(V-1)
        HH(m)=HH(m-1)+lambda(m)*(S(m)-S(m-1));
    end
    kl=zeros(N,1);
    for i=1:N
        t=rt(i,:);
        ORi=OR(i,:);
        ht=zeros(1,J);
        for j=1:J
            k=ORi(j);
            if k==1
                ht(j)=t(j);
            else
                ht(j)=lambda(k)*(t(j)-S(k-1))+HH(k-1);
            end
        end
        fun = @ (x)  exp(log(normpdf(x,0,1))+sum(phi.*(x-zeta)+log(lambda(ORi))-(1+1./gamma).*log(1+gamma.*exp(phi.*(x-zeta)).*ht))-Ln_hat(i));
        h=integral(fun,-Inf,Inf,'ArrayValued',true);
        kl(i)=log(h)+Ln_hat(i);
        
    end
    logf_rt_give_res(r,:)=kl;
end
%joint
dim=size(samples);
niter=dim(1);
mcmc_theta=samples(:,1:N);
mcmc_tau=samples(:,(N+1):(2*N));
mcmc_a=samples(:,(2*N+1):(2*N+J));
mcmc_b=samples(:,(2*N+J+1):(2*N+2*J));
mcmc_gamma=samples(:,(2*N+2*J+1):(2*N+3*J));
mcmc_phi=samples(:,(2*N+3*J+1):(2*N+4*J));
mcmc_zeta=samples(:,(2*N+4*J+1):(2*N+5*J));
mcmc_varphi=samples(:,2*N+5*J+1);
mcmc_lambda=samples(:,(dim(2)-65:(dim(2)-1)));

eap_a=mean(mcmc_a,1);
eap_b=mean(mcmc_b,1);
eap_gamma=mean(mcmc_gamma,1);
eap_lambda=mean(mcmc_lambda,1);
eap_phi=mean(mcmc_phi,1);
eap_th=mean(mcmc_theta,1);
eap_tau=mean(mcmc_tau,1);
eap_varphi=mean(mcmc_varphi);
eap_zeta=mean(mcmc_zeta,1);
logf=zeros(niter,N);
logf3=logf;
logf4=logf;
rng(1)
parfor r=1:niter
    
    a=mcmc_a(r,:);
    b=mcmc_b(r,:);
    gamma=mcmc_gamma(r,:);
    lambda=mcmc_lambda(r,:);
    phi=mcmc_phi(r,:);
    tau=mcmc_tau(r,:);
    th=mcmc_theta(r,:);
    varphi=mcmc_varphi(r);
    zeta=mcmc_zeta(r,:);
    fyt=zeros(N,1);
    for i=1:N
        
        res=resp(i,:);
        t=rt(i,:);
        x0=[th(i),tau(i)];
        options = optimoptions(@fminunc,'Display','none');
        [mode,fval]=fminunc(@(x) NLogf(x,t,res,a,b,phi,zeta,gamma,varphi,lambda,S,OR(i,:)),x0,options);
        sigma=Var_Cov(mode,t,a,b,phi,zeta,gamma,varphi,lambda,S,OR(i,:));
        R=mvnrnd(mode,sigma,1000); %B*2
        kb=mvnpdf(R,mode,sigma);%B*1
        fb=JointpdfMatrix(R,t,res,a,b,phi,zeta,gamma,varphi,lambda,S,OR(i,:));
        fk=fb./kb;
        fyt(i)=mean(fk); 
    end
    logf(r,:)=log(fyt);%joint
    logf3(r,:)=logf(r,:)-logf_p_rt(r,:);
    logf4(r,:)=logf(r,:)-logf_p_resp(r,:);
end

H = exp(logf);
logPPO = log(mean(H,1));
lppd = sum(logPPO);
p_total = 2 * sum((logPPO - mean(logf,1)))
WAIC_total = -2 * (lppd - p_total)

H = exp(logf_res_give_rt);
logPPO = log(mean(H,1));
lppd = sum(logPPO);
p_res = 2 * sum((logPPO - mean(logf_p_resp,1)))
WAIC_res = -2 * (lppd - p_res)

H = exp(logf_rt_give_res);
logPPO = log(mean(H,1));
lppd = sum(logPPO);
p_rt= 2 * sum((logPPO - mean(logf_p_rt,1)))
WAIC_rt = -2 * (lppd - p_rt)

H = exp(logf);
Ht=exp(logf_res_give_rt);
logPPO = log(mean(H,1))-log(mean(Ht,1));
lppd = sum(logPPO);
p_rt_give_res = 2 * sum((logPPO - mean(logf4,1)))
WAIC_rt_give_res = -2 * (lppd - p_rt_give_res)

H = exp(logf);
Ht=exp(logf_rt_give_res);
logPPO = log(mean(H,1))-log(mean(Ht,1));
lppd = sum(logPPO);
p_resp_give_rt = 2 * sum((logPPO - mean(logf3,1)))
WAIC_resp_give_rt = -2 * (lppd - p_resp_give_rt)


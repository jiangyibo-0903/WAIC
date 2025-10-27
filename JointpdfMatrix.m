function  fn= JointpdfMatrix(X,t,y,a,b,phi,zeta,gamma,varphi,lambda,S,ORi)
%joint density of X=[th,tau],t,y given Omega.
dim=size(X);
KL=dim(1);
th=X(:,1);
tau=X(:,2);
ein=th*a-ones(KL,1)*(a.*b);
respi=ones(KL,1)*y;
ein1=sin(varphi)*(th*phi)+cos(varphi)*(tau*phi)-ones(KL,1)*(phi.*zeta);
M=zeros(1,length(lambda));
M(1)=lambda(1)*S(1);
for k=2:(length(lambda)-1)
   M(k)=M(k-1)+lambda(k)*(S(k)-S(k-1)); 
end
H=zeros(length(t),1);
for j=1:length(t)
    v=ORi(j);
    if(v==1)
        H(j)=t(j);
    else
        H(j)=lambda(v)*(t(j)-S(v-1))+M(v-1);
    end
end
fn=respi.*ein-log(1+exp(ein))+ones(KL,1)*log(lambda(ORi))+ein1-(1+ones(KL,1)*(1./gamma)).*log(1+(ones(KL,1)*(gamma.*H')).*exp(ein1));
fn=normpdf(tau).*normpdf(th).*exp(sum(fn,2));
end


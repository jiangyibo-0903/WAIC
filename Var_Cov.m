function fn = Var_Cov(x,t,a,b,phi,zeta,gamma,varphi,lambda,S,ORi)
%Variance-Covariance matrix based on the -inverse of the 
%  second derivatives of logf
fn=zeros(2,2);
z=zeros(2,2);
z(1,1)=1;
z(2,2)=1;
eta=sin(varphi);
ap=cos(varphi);
ein=a*x(1)-a.*b;
ex=a.^2.*exp(ein)./((1+exp(ein)).^2);
ein1=phi*(x(1)*eta+x(2)*ap)-phi.*zeta;
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
ep=phi.^2.*(1+gamma).*H'.*exp(ein1)./((1+gamma.*H'.*exp(ein1)).^2);
z(1,1)=z(1,1)+sum(ex)+eta^2*sum(ep);
z(2,2)=z(2,2)+ap^2*sum(ep);
h=eta*ap*sum(ep);
z(1,2)=h;
z(2,1)=z(1,2);
det=z(1,1)*z(2,2)-z(1,2)*z(2,1);
fn(1,1)=z(2,2)/det;
fn(2,2)=z(1,1)/det;
fn(1,2)=-z(1,2)/det;
fn(2,1)=fn(1,2);
end


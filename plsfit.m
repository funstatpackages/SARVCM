%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [beta,gamma,Yhat,sse]=plsfit(B,Q2,K,lambda,X,Y)

n=length(Y);
np=size(X,2);
nq=size(Q2,2);
BQ2=B*Q2;
Hmtx=0;
Svec=0;
for(i=1:n)
    %Bi=B(i,:);
    Xi=X(i,:);
    Yi=Y(i);
    %BiQ2=Bi*Q2;
    BiQ2=BQ2(i,:);
    XiBiQ2=kron(Xi,BiQ2);
    Hmtx=Hmtx+XiBiQ2'*XiBiQ2;
    Svec=Svec+XiBiQ2'*Yi;
end;
Lam=diag(lambda);
Dlam=kron(Lam,Q2'*K*Q2);
theta=inv(Hmtx+Dlam)*Svec;
%gamma=[];
%beta=[];
%for (ip=1:np)
%    thetak=theta(((ip-1)*length(theta)/np+1):(ip*length(theta)/np));
%    gammak=Q2*thetak;
%    gamma=[gamma gammak];
%    betak=B*gammak;
%    beta=[beta betak];
%end;

thetamtx=reshape(theta,nq,np);
gamma=Q2*thetamtx;
beta=B*Q2*thetamtx;

Yhat=sum(X.*beta,2);
sse=sum((Y-Yhat).^2);
end
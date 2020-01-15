%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [beta,gamma,Yhat,sse]=plsfitCV(B,Q2,K,lambda,X,Y)

n=length(Y);
np=size(X,2);
Hmtx=0;
Svec=0;
for (i=1:n)
    Bi=B(i,:);
    Xi=X(i,:);
    Yi=Y(i);
    BiQ2=Bi*Q2;
    XiBiQ2=kron(Xi,BiQ2);
    Hmtx=Hmtx+XiBiQ2'*XiBiQ2;
    Svec=Svec+XiBiQ2'*Yi;
end;

nl=size(lambda,1);
gamma_all=[];
beta_all=[];
Yhat_all=[];
sse_all=[];

for(il=1:nl)
    Lam=diag(lambda(il,:));
    Dlam=kron(Lam,Q2'*K*Q2);
    theta=inv(Hmtx+Dlam)*Svec;
    gamma=[];
    beta=[];
    for (ip=1:np)
        thetak=theta(((ip-1)*length(theta)/np+1):(ip*length(theta)/np));
        gammak=Q2*thetak;
        gamma=[gamma gammak];
        betak=B*gammak;
        beta=[beta betak];
    end;
    Yhat=sum(X.*beta,2);
    sse=sum((Y-Yhat).^2);
    sse_all=[sse_all sse];
end;
end
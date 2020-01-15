%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta gamma beta Yhat sse gcv df lambda_c]=...
    plsfitGCV(B,Q2,K,X_l,X_nl,Y,lambda);
% grid search for best penalty parameter lambda
J=size(Q2,2);
BQ2=B*Q2;

n=length(Y);
n_l=size(X_l,2);
n_nl=size(X_nl,2);

XB=[];
for(i=1:n)
	Xi=X_nl(i,:);
	Yi=Y(i);
	BiQ2=BQ2(i,:);
	XBi=kron(Xi,BiQ2);
	XB=[XB; XBi];
end;
W=[X_l XB];
%rhs=W'*Y;
WW=W'*W;
P=Q2'*K*Q2;

[~,flag]=chol(WW);
if(flag==0)
    Ainv=chol(WW,'upper');
    % test=Ainv'*Ainv;
    A=inv(Ainv');
    D=zeros((n_l+n_nl*J),(n_l+n_nl*J));
    D((n_l+1):(n_l+n_nl*J),(n_l+1):(n_l+n_nl*J))=kron(eye(n_nl),P);
    ADA=A*D*A';
    C=eig(ADA);
end;

lam_mtx=repmat(lambda',1,n_nl);
if(n_nl==0)
    lam_mtx=zeros(1,n_nl);
end;
nlam=size(lam_mtx,1);

theta_all=[];
gamma_all=[];
beta_all=[];
Yhat_all=[];
sse_all=[];
df_all=[];
gcv_all=[];
for(il=1:nlam)
    Lam=diag(lam_mtx(il,:));
    Dlam=zeros((n_l+n_nl*J),(n_l+n_nl*J));
    Dlam((n_l+1):(n_l+n_nl*J),(n_l+1):(n_l+n_nl*J))=kron(Lam,P);
    lhs=WW+Dlam;
    tmp=inv(lhs)*W';
    theta=tmp*Y;
    theta_l=theta(1:n_l);
    theta_all=[theta_all theta_l];
    theta_nl=theta((n_l+1):(n_l+n_nl*J));
    theta_mtx=reshape(theta_nl,J,n_nl);
    gamma=Q2*theta_mtx;
    gamma_all=[gamma_all gamma];
    beta=B*gamma;
    beta_all=[beta_all beta];
    Yhat=X_l*theta_l+sum(X_nl.*beta,2);
    Yhat_all=[Yhat_all Yhat];
    sse=sum((Y-Yhat).^2);
    sse_all=[sse_all sse];
    if(flag==0)
        if(n_nl>0)
            df=sum(1./(1+C*mean(lam_mtx(il,:))));
        end;
        if(n_nl==0)
            df=sum(1./(1+C));
        end;
    end;
    if(flag>0)
        Hmtx=W*tmp;
        df=sum(trace(Hmtx));
    end;
    df_all=[df_all df];
    gcv=n*sse/(n-df)^2;
    gcv_all=[gcv_all gcv];
end;
[gcv j]=min(gcv_all);
lambda_c=lam_mtx(j);
theta=theta_all(:,j);
gamma=gamma_all(:,(n_nl*(j-1)+1):(n_nl*j));
beta=beta_all(:,(n_nl*(j-1)+1):(n_nl*j));
Yhat=Yhat_all(:,j);
sse=sse_all(j);
gcv=gcv_all(j);
df=df_all(j);
end
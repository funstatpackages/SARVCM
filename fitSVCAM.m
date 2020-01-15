function [theta_l gamma beta alpha_hat...
    Ystar_hat Yhat mse mle_max aic bic]=...
    fitSVCAM(X_l,X_nl,Y,Wg,B,Q2,K,lambda,penalty)
n=length(Y);
n_l=size(X_l,2);
n_nl=size(X_nl,2);
alpha_init=(1:9)./10;
nalpha=length(alpha_init);

if ~penalty
    lambda=0;
end

% Step 0. Setup
if n_nl>0
    J=size(Q2,2);
    BQ2=B*Q2;
    XB=[];
    for i=1:n 
        Xi=X_nl(i,:);
        BiQ2=BQ2(i,:);
        XBi=kron(Xi,BiQ2);
        XB=[XB; XBi];
        P=Q2'*K*Q2;
    end
end
if n_nl==0
    XB=[]; BQ2=[]; Q2=[]; P=[];
    J=0;
end
W=[X_l XB];
WW=W'*W;

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
if(flag>0)
    VV=WW+diag([zeros(n_l,1); ones(n_nl*J,1)*(1e-12)]);
    [~,flag2]=chol(VV);
    if(flag2==0)
    Ainv=chol(VV,'upper');
    end;
    if(flag2>0)
    VV=WW+diag([zeros(n_l,1); ones(n_nl*J,1)*(1e-10)]);
    Ainv=chol(VV,'upper');
    end;
    A=inv(Ainv');
    D=zeros((n_l+n_nl*J),(n_l+n_nl*J));
    D((n_l+1):(n_l+n_nl*J),(n_l+1):(n_l+n_nl*J))=kron(eye(n_nl),P);
    ADA=A*D*A';
    C=eig(ADA);
end;

lam_mtx=repmat(lambda',1,n_nl);
if(n_nl==0)
    lambda=0;
    lam_mtx=zeros(1,n_nl);
end;
nlam=size(lam_mtx,1);

Ystar_all=[];
detA_all=[];
for ia=1:nalpha
	Ai=eye(n)-alpha_init(ia)*Wg;
    detA=log(det(Ai));
	Ystar_all=[Ystar_all Ai*Y];
    detA_all=[detA_all detA];
end;
rhs_all=W'*Ystar_all;

% Step 1. First Round Selection
sse_all=zeros(nlam,nalpha);
df_all=zeros(nlam,nalpha);
for il=1:nlam
    Lam=diag(lam_mtx(il,:));
    Dlam=zeros((n_l+n_nl*J),(n_l+n_nl*J));
    Dlam((n_l+1):(n_l+n_nl*J),(n_l+1):(n_l+n_nl*J))=kron(Lam,P);
    lhs=WW+Dlam;
    if n_nl>0
        df=sum(1./(1+C*mean(lam_mtx(il,:))));
    end;
    if n_nl==0
        df=sum(1./(1+C));
    end;    
    theta=lhs\rhs_all;
    theta_l=theta(1:n_l,:); 
    theta_nl=theta((n_l+1):(n_l+n_nl*J),:);
    
    for ia=1:nalpha
        theta_mtx=reshape(theta_nl(:,ia),J,n_nl);
        beta=BQ2*theta_mtx;
        if n_nl>0
            Yhat=X_l*theta_l(:,ia)+sum(X_nl.*beta,2);
        else
            Yhat=X_l*theta_l(:,ia);
        end;
        sse=sum((Ystar_all(:,ia)-Yhat).^2);
        sse_all(il,ia)=sse;
    end;
    df_all(il,:)=df;   
end;
if n_nl>0
    gcv_all=n*sse_all./(n-df_all).^2;
    if ~penalty
        ind=repmat(1,1,nalpha);
    else
        [gcv ind]=min(gcv_all);
    end
    sse=[];
    for ii=1:nalpha
        sse=[sse sse_all(ind(ii),ii)];
    end;
else
    sse=sse_all; gcv=[]; 
end;
mle_all=detA_all-n*log(sqrt(sse/n));
[mle_max ind_alpha]=max(mle_all);
if n_nl>0
    ind_lam=ind(ind_alpha);
else
    ind_lam=1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2. Second Round Selection
alpha_hat0=alpha_init(ind_alpha)
st=alpha_hat0-0.09;
ed=alpha_hat0+0.09;
alpha_init=[st:0.01:ed];
nalpha=length(alpha_init);

Ystar_all=[];
detA_all=[];
for(ia=1:nalpha)
	Ai=eye(n)-alpha_init(ia)*Wg;
    detA=log(det(Ai));
	Ystar_all=[Ystar_all Ai*Y];
    detA_all=[detA_all detA];
end;
rhs_all=W'*Ystar_all;

sse_all=zeros(nlam,nalpha);
df_all=zeros(nlam,nalpha);
for il=1:nlam
    Lam=diag(lam_mtx(il,:));
    Dlam=zeros((n_l+n_nl*J),(n_l+n_nl*J));
    Dlam((n_l+1):(n_l+n_nl*J),(n_l+1):(n_l+n_nl*J))=kron(Lam,P);
    lhs=WW+Dlam;
    if n_nl>0
        df=sum(1./(1+C*mean(lam_mtx(il,:))));
    end;
    if n_nl==0
        df=sum(1./(1+C));
    end;  
    theta=lhs\rhs_all;
    theta_l=theta(1:n_l,:); 
    theta_nl=theta((n_l+1):(n_l+n_nl*J),:);
    
    for ia=1:nalpha
        theta_mtx=reshape(theta_nl(:,ia),J,n_nl);
        beta=BQ2*theta_mtx;
        if n_nl>0
            Yhat=X_l*theta_l(:,ia)+sum(X_nl.*beta,2);
        else
            Yhat=X_l*theta_l(:,ia);
        end;
        sse=sum((Ystar_all(:,ia)-Yhat).^2);
        sse_all(il,ia)=sse;
    end;
    df_all(il,:)=df;   
end;
if n_nl>0
    gcv_all=n*sse_all./(n-df_all).^2;
    if ~penalty
        ind=repmat(1,1,nalpha);
    else
        [gcv ind]=min(gcv_all);
    end
    sse=[];
    for ii=1:nalpha
        sse=[sse sse_all(ind(ii),ii)];
    end;
else
    sse=sse_all; gcv=[];
end;
mle_all=detA_all-n*log(sqrt(sse/n));
[mle_max ind_alpha]=max(mle_all);
if n_nl>0
    ind_lam=ind(ind_alpha);
else
    ind_lam=1;
end;

% Step 3. Final Fitting
alpha_hat=alpha_init(ind_alpha)
lam_c=lambda(ind_lam);
Ystar=Ystar_all(:,ind_alpha);
Ahat=eye(n)-alpha_hat*Wg;
detA=detA_all(ind_alpha);
Lam=diag(lam_mtx(ind_lam,:));
Dlam=zeros((n_l+n_nl*J),(n_l+n_nl*J));
Dlam((n_l+1):(n_l+n_nl*J),(n_l+1):(n_l+n_nl*J))=kron(Lam,P);
lhs=WW+Dlam;
rhs=W'*Ystar;
if n_nl>0
    df=sum(1./(1+C*mean(lam_mtx(ind_lam,:))));
end;
if n_nl==0
    df=sum(1./(1+C));
end;
theta=lhs\rhs;
theta_l=theta(1:n_l); 
theta_nl=theta((n_l+1):(n_l+n_nl*J));
theta_mtx=reshape(theta_nl,J,n_nl);
gamma=Q2*theta_mtx;
beta=BQ2*theta_mtx;
if n_nl>0
    Ystar_hat=X_l*theta_l+sum(X_nl.*beta,2);
else
    Ystar_hat=X_l*theta_l;
end;
sse=sum((Ystar-Ystar_hat).^2); 
mse_star=sse/n;
aic=df-mle_max;
bic=df*log(n)/2-mle_max;
Yhat=Ahat\Ystar_hat;
mse=mean((Y-Yhat).^2);
end
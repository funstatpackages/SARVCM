%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta_hat_all=beta_pred(Z,V,T,d,gamma_hat,theta_hat,ind_nl,ind_l)
npop=size(Z,1);
np=length(ind_nl)+length(ind_l);
% Generate Berstein Basis;
uu=Z(:,1);
vv=Z(:,2);
[Ind,cnt,Lam]=HQgetInd(V,T,uu,vv);
t=size(T,1); %number of triangles
m=(d+1)*(d+2)/2;
c=repmat(eye(m),t,1);
A=zeros(length(Ind),m);
Ind2=(1:length(Ind))';
for j=1:m
    A(:,j)=seval2(V,T,c(:,j),length(Ind),Ind2,cnt,Lam);
end
B=HQblkdiag(A,cnt);
B=sparse(B);
beta_r=B*gamma_hat;
beta_hat=NaN(npop,size(gamma_hat,2));
beta_hat(Ind,:)=beta_r;
beta_hat_all=zeros(npop,np);
beta_hat_all(:,ind_nl)=beta_hat;
beta_hat_all(:,ind_l)=repmat(theta_hat',npop,1);
end
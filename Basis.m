function [B B0 Ind]=Basis(Z,V,T,d);
% Generate Berstein Basis;
n=size(Z,1);
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
B0=HQblkdiag(A,cnt);
B=zeros(n,m*t);
B(Ind,:)=B0;
B0=sparse(B0);
B=sparse(B);
end
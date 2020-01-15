%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [beta_hat,gamma_hat,Y_hat,sse_hat,lambda_c,Ind]=...
    GCVlambda(Z,X,Y,V,T,Q2,K,d,lambda);

% grid search for best penalty parameter lambda
p=size(X,2);
lambda=repmat(lambda',1,p);
nl=size(lambda,1);

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

% 5-fold Cross-Validation;
Yi=Y(Ind);
Xi=X(Ind,:);
[n1,n2]=size(Yi);
np=size(X,2);
indp=randperm(n1);
sse=zeros(nl,1);
for(il=1:nl)
    for(ik=1:5)
        nk=floor(n1/5);
        if(ik<5)
            indk=indp(((ik-1)*nk+1):(ik*nk));
        end;
        if(ik==5)
            indk=indp(((ik-1)*nk+1):n1);
        end;
        indnk=setdiff(indp,indk);
        Yk=Yi(indk,:);
        Ynk=Yi(indnk,:);
        Xk=Xi(indk,:);
        Xnk=Xi(indnk,:);
        Bk=B(indk,:);
        Bnk=B(indnk,:);
        [betank,gamma,Ynkhat,ssenk]=...
            plsfit(Bnk,Q2,K,lambda(il,:),Xnk,Ynk); 
        betak=[];
        for (ip=1:np)
            gammap=gamma(:,ip);
            betap=Bk*gammap;
            betak=[betak betap];
        end;
        Ykhat=sum(Xk.*betak,2);
        sse(il)=sse(il)+sum((Yk-Ykhat).^2);
    end;
end;
%imagesc(Ykhat(1,:));
plot(sse);
[J,j]=min(sse); 
lambda_c=lambda(j,:);
[beta_hat,gamma_hat,Y_hat,sse_hat]=...
    plsfit(B,Q2,K,lambda_c,Xi,Yi); 
end
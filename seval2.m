function Z = seval2(V,T,c,sz,IndP,cnt,Lam)
%        Z = seval2(V,T,c,sz,IndP,cnt,Lam)
% This function evaluates the piecewise polynomials over the triangulation 
% [V,T] defined by the vector c of Bnet coeffs at the points with coords X 
% and Y .
% This matlab program is written by Ming-Jun Lai and Qianying Hong
% based on a version copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc.. 
k = size(T,1);
Z = nan*zeros(sz,1);
m = length(c)/k;
d=degree(m);
Ind1=zeros(1/6*d*(d+1)*(d+2),1);
Ind2=Ind1;
Ind3=Ind1;
cnt1=zeros(d+1,1);
for i=1:d
    cnt1(i+1)=cnt1(i)+(d-i+1)*(d-i+2)/2;
end
for i=1:d
    [I,J,K] = indices(d-i+1);
    [I1,J1,K1] = indices(d-i);
    Ind1(cnt1(i)+1:cnt1(i+1)) = locate(I1+1,J1,K1,I,J,K);
    Ind2(cnt1(i)+1:cnt1(i+1)) = locate(I1,J1+1,K1,I,J,K);
    Ind3(cnt1(i)+1:cnt1(i+1)) = locate(I1,J1,K1+1,I,J,K);
end

C=reshape(c,[m,k]);
B1=HQblkdiag(Lam(:,1),cnt);
B2=HQblkdiag(Lam(:,2),cnt);
B3=HQblkdiag(Lam(:,3),cnt);
B1=sparse(B1);
B2=sparse(B2);
B3=sparse(B3);
C=C(Ind1(cnt1(1)+1:cnt1(1+1)),:)*B1'+C(Ind2(cnt1(1)+1:cnt1(1+1)),:)*B2'+C(Ind3(cnt1(1)+1:cnt1(1+1)),:)*B3';


%initial C
% 
% for j=1:k
%     C=[C,repmat(c((j-1)*m+1:j*m),1,cnt2(j+1)-cnt2(j))];
% end

if d>=2
    n=length(IndP);
    for i=2:d
        C=C(Ind1(cnt1(i)+1:cnt1(i+1)),:)*spdiags(Lam(:,1),0,n,n)+C(Ind2(cnt1(i)+1:cnt1(i+1)),:)*spdiags(Lam(:,2),0,n,n)...
            +C(Ind3(cnt1(i)+1:cnt1(i+1)),:)*spdiags(Lam(:,3),0,n,n);
    end
end 
Z(IndP)=C;

function d = degree(m)
%        d = degree(m)
% d = the degree of the polynomial corresponding to the vector 
% with length m
d = (-3 + sqrt(8*m+1))/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: HQgetInd.m
% Function: [Ind,cnt,Lam]=HQgetInd(V,T,xx,yy)
%
% Input:
% V: vertices
% T: triangles
% xx & yy: latitude & longitude
% 
% Output:
% Ind:
% cnt:
% 
% This matlab program is copyrighted @2010 by Ming-Jun Lai and Qianying Hong.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ind,cnt,Lam]=HQgetInd(V,T,xx,yy)
xx=xx(:);
yy=yy(:);
Ind=[];
k=size(T,1); %number of triangles
cnt=zeros(k+1,1);
tol=100*eps; %some small number
L=HQbary(V,T,xx,yy); %the barycentric coordinates with respect to the triangle 
    % specified by [V,T] of the points (xx,yy);
Lam=[];
for j = 1:k
   I=find(L((j-1)*3+1,:)>=-tol & L((j-1)*3+2,:)>=-tol & L((j-1)*3+3,:)>=-tol);
   Lam=[Lam;[L((j-1)*3+1,I)',L((j-1)*3+2,I)',L((j-1)*3+3,I)']];
   Ind=[Ind;I'];
   cnt(j+1)=cnt(j)+length(I);
end;


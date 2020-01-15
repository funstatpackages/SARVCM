function [I,n]=inside(V,T,X,Y);
%This program determines which (X,Y) are inside of triangulation (V,T)
%and returns the index I and the number of points (X,Y) inside of (V,T);
m=size(T,1);
I=[];
for i=1:m
   V1=V(T(i,1),:);V2=V(T(i,2),:);V3=V(T(i,3),:);
   [lam1,lam2,lam3] = bary(V1,V2,V3,X,Y);
   J=find(lam1>0 & lam2>0 &lam3>0);
   I=[I;J];
end
n=size(I,1);

function Ind=insideVT(V,T,xx,yy)
%   Ind=insideVT(V,T,xx,yy)
%This function finds points (xx,yy) inside of triangulation V,T;
n=size(xx(:),1);
T=mkcc(V,T);
m=size(T,1); Ind=zeros(n,1);
for i=1:m
    v1=V(T(i,1),:); v2=V(T(i,2),:);v3=V(T(i,3),:);
    [lam1,lam2,lam3]=bary(v1,v2,v3,xx(:),yy(:));
    index=find(lam1>=0 & lam2>=0 & lam3>=0);
    Ind(index)=ones(size(index));
end
end
    
function [lam1,lam2,lam3] = bary(V1,V2,V3,X,Y);
%        [lam1,lam2,lam3] = bary(V1,V2,V3,X,Y);
% This function returns the barycentric coordinates with respect to the 
% triangle [V1,V2,V3] of the points (X,Y);
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
One = ones(size(X(:)'));
A = [1 1 1 ;V1(1),V2(1),V3(1);V1(2),V2(2),V3(2)];
lam = A\[One;X(:)';Y(:)'];
lam1 = reshape(lam(1,:),size(X));
lam2 = reshape(lam(2,:),size(X));
lam3 = reshape(lam(3,:),size(X));
end

function T = mkcc(V,T);
%        T = mkcc(V,T)
% If necessary, this function re-orders in counter clockwise order
% the vertices of each triangle in the given triangulation.
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
    m = size(T,1);
    for i = 1:m
      if triarea(V(T(i,1),:),V(T(i,2),:),V(T(i,3),:)) < 0
        T(i,:) = [T(i,1),T(i,3),T(i,2)];
      end;
    end; 
end
    
function  A = triarea(V1,V2,V3)
%         A = triarea(V1,V2,V3)
%This is a matlab program to compute the signed area of a triangle with vertices 
% V1, V2, V3. 
x = V1(:,1); y = V1(:,2);
a = V2(1);  b = V2(2);
c = V3(1);  d = V3(2);
     A=(a-x).*(d-y)-(c-x).*(b-y);
     A=A/2;
end
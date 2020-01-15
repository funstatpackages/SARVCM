%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: HQbary.m
% Function: Lam = HQbary(V,T,X,Y)
% This function returns the barycentric coordinates with respect to the 
% triangle specified by [V,T] of the points (X,Y);
%
% Input:
% V: vertices
% T: triangles
% X & Y: latitude & longitude
% 
% Output:
% Lam: the output Lam has dimension 3*T-points
% 
% This matlab program is copyrighted @2010 by Ming-Jun Lai.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Lam= HQbary(V,T,X,Y);
k=size(T,1);
A=T';
B=V(A(:,:),:);
One=ones(size(B,1),1);
C=[One,B];
cnt=3*(0:k)';

C=HQblkdiag(C,cnt);
C=sparse(C');

One=ones(size(X(:)'));
D=[One;X(:)';Y(:)'];
D=repmat(D,k,1);

Lam=C\D;

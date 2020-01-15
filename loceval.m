function y = loceval(lam1,lam2,lam3,Bcoeff);
%        y = loceval(lam1,lam2,lam3,Bcoeff);
% This function evaluates the poly with Bnet coefficients in the vector Bcoeff 
% at the points with barycentric coordinates lam1, lam2, lam3.
m = length(Bcoeff);
d = degree(m);
for j = 1:d
   Bcoeff = de_cast_step(lam1,lam2,lam3,Bcoeff);
end;
y = Bcoeff';

loceval(lam1(2601),lam2(2601),lam3(2601),c(((j-1)*m+1):(j*m)))

function Bout = de_cast_step(lam1,lam2,lam3,Bin);
%        Bout = de_cast_step(lam1,lam2,lam3,Bin);
% This function does one step of the de Casteljau algorithm
% Note: This function assumes that the input arguments are all column vectors
m = size(Bin,1);
d = degree(m);
n = length(lam1);
[I,J,K] = indices(d);
[I1,J1,K1] = indices(d-1);
indx1 = locate(I1+1,J1,K1,I,J,K);
indx2 = locate(I1,J1+1,K1,I,J,K);
indx3 = locate(I1,J1,K1+1,I,J,K);
if size(Bin,2) == 1
   Bout = Bin(indx1)*lam1' + Bin(indx2)*lam2' + Bin(indx3)*lam3';
else
   Bout = Bin(indx1,:)*spdiags(lam1,0,n,n) + Bin(indx2,:)*spdiags(lam2,0,n,n) + ...
      Bin(indx3,:)*spdiags(lam3,0,n,n);
end;

function d = degree(m);
%        d = degree(m)
% d = the degree of the polynomial corresponding to the vector 
% with length m
d = (-3 + sqrt(8*m+1))/2;

function [I,J,K] = indices(d)
%        [I,J,K] = indices(d)
% This function computes the index vectors I,J,K associated with our linear ordering
% for the Bnet coefficients of a polynomial of degree d over an arbitrary triangle.
m = (d+1)*(d+2)/2;
I = zeros(m,1);
J = I;
K = I;
Mark = 1;
for j = d:(-1):0
   I(Mark:(Mark+j)) =[j:(-1):0]';
   J(Mark:(Mark+j)) =[0:j]';
   K(Mark:(Mark+j)) =(d-j)*ones(j+1,1);
   Mark = Mark+j+1;
end;

function Index = locate(I1,J1,K1,I,J,K);
%        Index = locate(I1,J1,K1,I,J,K);
% Index = the location of the index [I1,J1,K1] in the list 
% [I,J,K]
Index = ismember([I,J,K],[I1,J1,K1],'rows');
Index = find(Index);


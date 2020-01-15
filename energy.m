function K = energy(V,T,d);
%        K = energy(V,T,d);
% This function computes the energy  matrix associated with the minimal energy norm.
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
n = size(T,1);
m = (d+1)*(d+2)/2;
msqr = m*m;
Mat = build(d-2);
Indx1 = zeros(n*msqr,1);
Indx2 = Indx1;
S = Indx1;
place = 1; 
for k = 1:n
   LocK = locEng(V(T(k,1),:),V(T(k,2),:),V(T(k,3),:),Mat,d);
   [i,j,s] = find(LocK);
   L = length(i);
   Indx1(place:(place + L-1)) = (k-1)*m + i;
   Indx2(place:(place + L-1)) = (k-1)*m + j;
   S(place:(place + L-1)) = s;
   place = place + L;
end;
K = sparse(Indx1(1:(place-1)),Indx2(1:(place-1)),S(1:(place-1)),n*m,n*m);

function K3 = locEng(V1,V2,V3,Mat,d)
%        K3 = locEng(V1,V2,V3,Mat,d)
% This function computes the local energy matrix.
m = (d+1)*(d+2)/2;
Id = eye(m);
vx = [1,0];
vy = [0,1];
[lam1x,lam2x,lam3x] = tcord(V1,V2,V3,vx);
[lam1y,lam2y,lam3y] = tcord(V1,V2,V3,vy);
Dx = dirder(Id,lam1x,lam2x,lam3x);
Dxx= dirder(Dx,lam1x,lam2x,lam3x);
Dxy= dirder(Dx,lam1y,lam2y,lam3y);
Dy = dirder(Id,lam1y,lam2y,lam3y);
Dyy = dirder(Dy,lam1y,lam2y,lam3y);
K3 = abs(triarea(V1,V2,V3))*(Dxx'*Mat*Dxx + 2*Dxy'*Mat*Dxy+ Dyy'*Mat*Dyy);

function [lam1,lam2,lam3] = tcord(V1,V2,V3,v);
%        [lam1,lam2,lam3] = tcord(V1,V2,V3,v)
% This function returns the T cords with respect to the 
% triangle with vertices V1,V2,V3.
[lam1,lam2,lam3] = bary(V1,V2,V3,v(1),v(2));
[a1,a2,a3] = bary(V1,V2,V3,0,0);
lam1 = lam1-a1;
lam2 = lam2-a2;
lam3 = lam3-a3;

function Mat = build(d)
%        Mat = build(d)
% This function builds the matrix needed to compute the inner product of two
% polynomials of degree d over an arbitrary triangle.
[I,J,K] = indices(d);
m = (d+1)*(d+2)/2;
Mat = zeros(m,m);
for j = 1:m
   for i = 1:m
      Mat(i,j) = choose(I(j),I(j)+I(i))*choose(J(j),J(j)+J(i))*choose(K(j),K(j)+K(i));
   end;
end;
Mat = Mat/(choose(d,2*d)*choose(2,2*d+2));

function y = choose(n,m)
%        y = choose(n,m)
% This function returns the number of ways of choosing n elements from a
% set of m elements. This function assumes that both n and m are nonnegative
% integers.
y = prod(1:m)/(prod(1:n)*prod(1:(m-n)));

function DerBcoeff = dirder(Bcoeff,lam1,lam2,lam3);
%        DerBcoeff = dirder(Bcoeff,lam1,lam2,lam3);
% Bcoeff = the Bnet coeffs over some triangle of a polynomial of 
% degree d. DerBcoeff = the Bnet coeffs of the directional derivative
% in the direction given by the vector having tcoords [lam1,lam2,lam3].
m = size(Bcoeff,1);
d = degree(m);
DerBcoeff = d*de_cast_step(lam1,lam2,lam3,Bcoeff);

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


function  A = triarea(V1,V2,V3)
%         A = triarea(V1,V2,V3)
%This is a matlab program to compute the signed area of a triangle with vertices 
% V1, V2, V3. 
x = V1(:,1); y = V1(:,2);
a = V2(1);  b = V2(2);
c = V3(1);  d = V3(2);
     A=(a-x).*(d-y)-(c-x).*(b-y);
     A=A/2;

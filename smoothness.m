%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: smoothness.m: function H = smoothness(V,T,d,r)
% This function returns the coefficient matrix for the system of eqs 
% that must be satisfied if the spline of degree d is to be C^r across 
% edges.
%
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul 
% Wenston through University of Georgia Research Foundation, Inc..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = smoothness(V,T,d,r)
[E,TE,TV,EV,B] = tdata(V,T);
int = find(sum(TE)>1);
n = length(int);
N = 0;
Neq = 0;
[I1,I2] = crarrays(d,r);
for j = 0:r
    N = N + ((j+1)*(j+2)/2+1)*(d+1-j);
    Neq = Neq + d + 1 -j;
    [LI,LJ,LK] = indices(j);
    I{j+1} = LI;
    J{j+1} = LJ;
    K{j+1} = LK;   
end;
% Neq = the number of equations generated in each triangle
% N = the number of nonzero coefficients in the equations generated in 
% each triangle.

m = (d+1)*(d+2)/2;
Index1 = zeros(N*n,1);
Index2 = zeros(N*n,1);
Values = Index1;
A=[1,2;2,3;3,1;2,1;3,2;1,3];
for j = 1:n
	k = int(j);
	v1 = E(k,1);
	v2 = E(k,2);
	AdjT = find(TE(:,k));
    t1 = AdjT(1);
    t2 = AdjT(2);
    T1 = T(AdjT(1),:);
    T2 = T(AdjT(2),:);
    i1 = find(T1 == v1);
    i2 = find(T1 == v2);
    j1 = find(T2 == v1);
    j2 = find(T2 == v2);
    e1 = find(ismember(A,[i1,i2],'rows'));
    e2 = find(ismember(A,[j1,j2],'rows'));
    if e1 > 3
        e1 = e1 - 3;
        Temp = T1;
        T1 = T2;
        T2 = Temp;
        Temp = e1;
        e1 = e2;
        e2 = Temp;
        Temp = t1;
        t1 = t2;
        t2 = Temp;
    else
        e2 = e2-3;
    end;
    v4 = setdiff(T2,[v1,v2]);
    V4 = V(v4,:);
    [lam1,lam2,lam3] = bary(V(T1(1),:),V(T1(2),:),V(T1(3),:),V4(1),V4(2));
    lambda = [lam1;lam2;lam3];
    switch e1
    case 2
        temp = lambda(1);
        lambda(1) = lambda(2);
        lambda(2) = lambda(3);
        lambda(3) = temp;
    case 3
        temp = lambda(2);
        lambda(2) = lambda(3);
        lambda(3) = temp;
    end;
    VarCT = 0;
    EqCt = 0;
    for i = 0:r
        Lambda = factorial(i)./(gamma(I{i+1}+1).*gamma(J{i+1}+1).*gamma(K{i+1}+1)).*...
               lambda(1).^I{i+1}.*lambda(2).^J{i+1}.*lambda(3).^K{i+1};
        T1mat = I1{i+1}{e1};
        T2vector = I2{i+1}{e2};
        numeq = size(T1mat,1);
        numVar = (size(T1mat,2)+1)*numeq;
        T1Values = ones(size(T1mat))*diag(Lambda);
        eqnums = (j-1)*Neq + EqCt + diag([1:numeq])*ones(numeq,size(T1mat,2)+1);
        Index1(((j-1)*N+VarCT+1):((j-1)*N+VarCT+numVar)) = eqnums;
        Index2(((j-1)*N+VarCT+1):((j-1)*N+VarCT+numVar)) = ...
                        [((t1-1)*m+T1mat),((t2-1)*m+T2vector)];
        Values(((j-1)*N+VarCT+1):((j-1)*N+VarCT+numVar)) = [T1Values(:);-ones(size(T2vector))];
        VarCT = VarCT + numVar;
        EqCt = EqCt + numeq;
    end;
end;
H = sparse(Index1,Index2,Values,n*Neq,size(T,1)*m);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: smoothness.m: function [I1,I2] = Crarrays(d,r);
% This function cell arrays needed to impose the C^r conditions 
%
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul 
% Wenston through University of Georgia Research Foundation, Inc..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I1,I2] = crarrays(d,r);  
for j = 0:r
   [J1,J2] = crcellarrays(d,j);
   I1{j+1} = J1;
   I2{j+1} = J2;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: smoothness.m: function [I1,I2] = CrCellArrays(d,r);
% This function returns the cell arrays of indices needed to implement 
% the C^r across an edge condition. Note that both I1 and I2 will be 
% 3 X 1 cell arrays.
%
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul 
% Wenston through University of Georgia Research Foundation, Inc..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I1,I2] = crcellarrays(d,r);   
[J1,J2] = cr_indices(r,d);
D1 = zeros(d+1,1);
D2 = zeros(d+1,1);
s1 = d+1;
s2 = 1;
for j = 1:(d+1)
	D1(j) = s1;
    s1 = s1 + d+1-j;
    D2(j) = s2;
    s2 = s2 + d + 2 - j;
end;
I2{1} = flipud(J2);
Temp = D1 - r;
I2{2} = flipud(Temp(1:(d+1-r)));
Temp = D2 + r;
I2{3} = Temp(1:(d+1-r));
I1{1} = J1;
Temp = zeros(size(J1));
for j = 0:r
    Temp(:,j+1)=D1((j+1):(d+1-r+j));
end;
loc = r+2;
back = r+1;
for j = 1:r
    for k = 0:(r-j)
        Temp(:,loc) = Temp(:,loc-back)-1;
        loc = loc+1;
    end;
    back = back-1;
end;
I1{2}=Temp;
Temp = zeros(size(J1));
for j = 0:r
    Temp(:,j+1)=D2((j+1):(d+1-r+j));
end;
loc = r+2;
back = r+1;
for j = 1:r
    for k = 0:(r-j)
        Temp(:,loc) = Temp(:,loc-back)+1;
        loc = loc+1;
    end;
    back = back-1;
end;
I1{3} = flipud(Temp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: smoothness.m: function [I1,I2] = Cr_indices(r,d);
% This function returns indices needed to formulate the equations between
% Bnet coefficients that are a consequence of the C^r condition
%
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul 
% Wenston through University of Georgia Research Foundation, Inc..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I1,I2] = cr_indices(r,d);
I1 = [];
Start = 1;
D = d + 1;
for j = 0:r
	for k = 0:(r-j)
        new_col = [(Start+k):(Start+k+d-r)]';
        I1 = [I1,new_col];
    end;
    Start = Start + D;
    D = D-1;
end;
I2 = [(-r*r/2+r*(d+3/2)+1):(-r*r/2+r*(d+1/2)+d+1)]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: smoothness.m: function [I,J,K] = indices(d)
% This function computes the index vectors I,J,K associated with our 
% linear ordering for the Bnet coefficients of a polynomial of degree 
% d over an arbitrary triangle.
%
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul 
% Wenston through University of Georgia Research Foundation, Inc..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,J,K] = indices(d)
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

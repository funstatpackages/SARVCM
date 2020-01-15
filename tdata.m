%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: tdata.m: function: [E,TE,TV,EV,B] = tdata(V,T)
% Input: 
% V: vertices
% T: triangles
%
% Output:
% Edges: a list of edges
% TE: Num(triangles) X Num(edges) matrix whose (i,j)th entry = 1 if 
% edge #j is an edge for triangle #i.
% TV: Num(triangles) X Num(vertices) matrix whose (i,j)th entry = 1 
% if vertex #j is a vertex for triangle #i.
% EV: Num(edges) X Num(vertices) matrix whose (i,j)th entry = 1 if 
% edge #i has vertex #j as an endpoint.
% Bdr: a list of the bdr vertice in cc order on the boundary
% BdrEdges: a list of the bdr edges in cc order on the boundary
%
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul 
% Wenston through University of Georgia Research Foundation, Inc..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Edges,TE,TV,EV,Bdr] = tdata(V,T);
[m,n] = size(T);
Edges = [];
TE = [];
numEdges = 0;
for i = 1:m
    Ti = T(i,:);
    for j = 1:3
        edge = [min([Ti(j) Ti(rem(j,3)+1)]),...
        max([Ti(j) Ti(rem(j,3)+1)])];
        if ~isempty(Edges)
            edgenum = find((edge(1) == Edges(:,1))&(edge(2) == Edges(:,2))); 
        else
            edgenum = [];
        end;
        if isempty(edgenum) 
            Edges = [Edges;edge];
            numEdges = numEdges + 1; 
            TE = [TE sparse(m,1)];
            edgenum = numEdges;
        end;
        TE(i,edgenum) = 1;
    end;
end;
numV = size(V);
numV = numV(1); 
TV = sparse(m,numV);
for i = 1:m
    TV(i,T(i,:)) = ones(1,3);
end; 
EV = sparse(numEdges,numV);
for i = 1:numEdges
    EV(i,Edges(i,:)) = ones(1,2);
end; 
Bdr = findbdt(T,V,Edges,TE,EV);

    
function B = findbdt(T,V,E,TE,EV);
%        B = findbdr(T,V,E,TE,EV);
% This function finds the boundary vertices for the above triangulation
% and lists them so that any two consecutive vertices are adjacent on
% the boundary.
 BE = find(sum(TE)==1);
 B = [];
%
% BE = a list of bdr edges, but not yet in cc order.
%
while ~isempty(BE)
    enum = BE(1);
    BEE = E(BE,:);
    e = E(enum,:);
    t = find(TE(:,enum));
    Tri = T(t,:);
    V1 = V(Tri(1),:);
    V2 = V(Tri(2),:);
    V3 = V(Tri(3),:);
    P = (V1 + V2 + V3)/3;
    VL = V(e(1),:);
    VR = V(e(2),:);
    if triarea(VL,VR,P) > 0
       B1 = e';
    else
       B1 = [e(2);e(1)]; 
    end;
    I1 = [1];
    loop = 0;
    i = 1;
    while ~loop
      v1 = B1(i);
      v2 = B1(i+1);
%
% [v1,v2] is the last bdr edge, with the cc direction from v1 to v2 
% 
      next_e = find(((v2==BEE(:,1))&(v1~=BEE(:,2)))|((v2==BEE(:,2))&(v1~=BEE(:,1))));
      if BEE(next_e,1) == v2 
        v3 = BEE(next_e,2);
      else
        v3 = BEE(next_e,1);
      end;
      if v3 ~= B1(1);
         B1 = [B1;v3];
         I1 = [I1;next_e];
         i = i+1;
      else
         I1 = [I1;next_e];
         loop = 1;
      end;
    end;   
    BE = del(BE,BE(I1));
    B = newcol(B,B1);
  end 
 
function A = newcol(B,c);
%        A = newcol(B,c)
% This function sets A = [B,c] with B or c padded with zeros as needed
  [n,k] = size(B);
  m = length(c);
  if n > m
    A = [B,[c;zeros(n-m,1)]];
  else
    A = [[B;zeros(m-n,k)],c];
  end;
  
 function [lnew,I] = del(l1,l2);
%        [lnew,I] = del(l1,l2)
% This function deletes the set of elements l2 from l1. Note: I = the 
% indices of those elements in l1 which are not in l2
  m = length(l2);
  n = length(l1);
  lnew = l1(:)'; 
  I = 1:n;
  for i = 1:m
    loc = find(l2(i) == lnew);
    if ~isempty(loc)
      k = length(loc);
      for j = 1:k
        J = find(l1 == l2(i));
        I(J) = zeros(size(J));
        lnew = [lnew(1:(loc(j)-1)),lnew((loc(j)+1):length(lnew))];
      end;
    end;
  end;
  I = I(find(I));
  n = size(l1);
  if n(1) > 1
    lnew = lnew';
  end;
     
function  A = triarea(V1,V2,V3)
%         A = triarea(V1,V2,V3)
%This is a matlab program to compute the signed area of a triangle with vertices 
% V1, V2, V3. 
x = V1(:,1); y = V1(:,2);
a = V2(1);  b = V2(2);
c = V3(1);  d = V3(2);
     A=(a-x).*(d-y)-(c-x).*(b-y);
     A=A/2;

 

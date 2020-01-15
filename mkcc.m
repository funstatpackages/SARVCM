%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: mkcc.m: function T = mkcc(V,T);
% This function re-orders in counter clockwise order the vertices of 
% each triangle in the given triangulation.
%
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul 
% Wenston through University of Georgia Research Foundation, Inc..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = mkcc(V,T);
    m = size(T,1);
    for i = 1:m
        if triarea(V(T(i,1),:),V(T(i,2),:),V(T(i,3),:)) < 0
            T(i,:) = [T(i,1),T(i,3),T(i,2)];
        end;
    end; 

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: mkcc.m: function  A = triarea(V1,V2,V3);
% This is a matlab program to compute the signed area of a triangle 
% with vertices V1, V2, V3. 
%
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul 
% Wenston through University of Georgia Research Foundation, Inc..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  A = triarea(V1,V2,V3)
    x = V1(:,1); y = V1(:,2);
    a = V2(1);  b = V2(2);
    c = V3(1);  d = V3(2);
    A=(a-x).*(d-y)-(c-x).*(b-y);
    A=A/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: refine.m: function [V1,T1] = refine(V0,T0);
% This function refines the input triangulation by adding nodes to 
% the midpoints of the sides of the input triangulation. 
%
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul 
% Wenston through University of Georgia Research Foundation, Inc..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V1,T1] = refine(V0,T0);
    n = size(V0,1);
    m = size(T0,1);
    V1 = V0;
    T1 = [];
    for i = 1:m
        J = T0(i,:);
        x1 = V0(J(1),1);  y1 = V0(J(1),2);
        x2 = V0(J(2),1);  y2 = V0(J(2),2);
        x3 = V0(J(3),1);  y3 = V0(J(3),2);
        x12 = (x1 + x2)/2;  y12 = (y1 + y2)/2;  
        x13 = (x1 + x3)/2;  y13 = (y1 + y3)/2;  
        x23 = (x2 + x3)/2;  y23 = (y2 + y3)/2;  

        % placeij = the location in V1 of (xij,yij)
        place12 = find((V1(:,1) == x12) & (V1(:,2) == y12));
        place13 = find((V1(:,1) == x13) & (V1(:,2) == y13));
        place23 = find((V1(:,1) == x23) & (V1(:,2) == y23));
        if isempty(place12)
            n = n+1;
            place12 =n;
            V1 = [V1;x12 y12];
        end;
        if isempty(place13)
            n = n+1;
            place13 =n;
            V1 = [V1;x13 y13];
        end;
        if isempty(place23)
            n = n+1;
            place23 =n;
            V1 = [V1;x23 y23];
        end;
        T1 = [T1;J(1) place12 place13];
        T1 = [T1;place12 place13 place23];
        T1 = [T1;J(2) place12 place23];
        T1 = [T1;J(3) place13 place23];
    end; 
T1 = mkcc(V1,T1);


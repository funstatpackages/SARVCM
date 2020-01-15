function Index = locate(I1,J1,K1,I,J,K);
%        Index = locate(I1,J1,K1,I,J,K);
% Index = the location of the index [I1,J1,K1] in the list 
% [I,J,K]
Index = ismember([I,J,K],[I1,J1,K1],'rows');
Index = find(Index);
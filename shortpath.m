%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SS: total points location for triangulation
% T_all: K by 3 matrix, denote triangles
% D: N by N matrix, show the paired distance between two points
% DG: N by N matrix, only show the distance between two points if they
%     share an edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist,dist_full,M]=shortpath(Vnew,Tnew,Z)
n=size(Z,1);
nV=size(Vnew,1);
nT=size(Tnew,1);

% find vertice in the same triangule;
M=zeros(nV,nV);
indT=1:nT;
for(i=1:nV)
	loc1=indT(sum(Tnew==i,2)>0);
    iV=Tnew(loc1,:);
    iV=iV(:);
    M(i,iV)=1;    
end;
M=M-diag(diag(M));

% calculate the distance;
dist_d=reshape(sqrt(sum((kron(Vnew',ones(1,nV))...
    -repmat(Vnew',1,nV)).^2)),nV,nV);

INF=1000*max(max(dist_d))*nV;

dist_d(M==0)=INF;
dist_d=dist_d-diag(diag(dist_d));
dist_d=min(dist_d,dist_d');
dist_full=dist_d;

% find the shortest path
for k=1:nV
    dist_d=min(dist_d,repmat(dist_d(:,k),1,nV)+...
        repmat(dist_d(k,:),nV,1));
end;
    
dist=dist_d(1:n,1:n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dg Wg]=GWeights2(V,T,bb,Z,~,~);
% find centroid of each triangle

K=size(T,1); % the number of triangles 
n=size(Z,1); % the number of locations
cent=[]; % K*2 matrix denote the centroid of triangle

for i=1:K
XY=V(T(i,:),:);
cent=[cent; mean(XY(:,1)) mean(XY(:,2))]; 
end; % find the centroids of each triangle

% calculate centroid distance matrix;
[cent_d Wg]=GWeights(V,T,bb,cent);

% find distance within the same triangle
dist_d=reshape(sqrt(sum((kron(Z',ones(1,n))...
    -repmat(Z',1,n)).^2)),n,n);

% find index of each point
[Ind,cnt,Lam]=HQgetInd(V,T,Z(:,1),Z(:,2));
ind_z=zeros(n,1);
for i=1:K
    st=cnt(i)+1;
    ed=cnt(i+1);
    ind_z(Ind(st:ed))=i;
end;

% calculate dist between point and centroid
cent_z=zeros(n,1);
for i=1:n
    cent_z(i)=sqrt(sum((Z(i,:)-cent(ind_z(i),:)).^2));
end;

% build up a matrix n*n
% (i,j)=(i,c1)+(j,c2)
% i-th point in triangle c1, j-th point in triangle c2
ptct=zeros(n,n);
for i=1:n
    for j=(i+1):n
        ptct(i,j)=cent_z(i)+cent_z(j);
    end;
end;
ptct=ptct'+ptct;

dist_z=dist_d;
n_all=1:n;

for i=1:K
     ind1=n_all(ind_z==i); % find points in the i-th triangle T_i
 for j=(i+1:K)
     ind2=n_all(ind_z==j); % find points in the j-th triangle T_j (i\neq j)
     dist_temp=ptct(ind1,ind2); 
     dist_z(ind1,ind2)=cent_d(i,j)+dist_temp; 
     dist_z(ind2,ind1)=cent_d(i,j)+dist_temp';
 end;
end;
Dg=dist_z;
Wg=weights(dist_z,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
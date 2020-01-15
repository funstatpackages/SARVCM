%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dg=dist0(V,T,bb,Z);
% find centroid of each triangle

K=size(T,1);
n=size(Z,1);

dist_d=reshape(sqrt(sum((kron(Z',ones(1,n))...
    -repmat(Z',1,n)).^2)),n,n);
INF=max(dist_d(:))*1000;
% find index of each point
[Ind,cnt,Lam]=HQgetInd(V,T,Z(:,1),Z(:,2));
ind_z=zeros(n,1);
for i=1:K
    st=cnt(i)+1;
    ed=cnt(i+1);
    ind_z(Ind(st:ed))=i;
end;


dist_z=dist_d;
n_all=1:n;

for i=1:K
     ind1=n_all(ind_z==i);
 for j=(i+1:K)
     ind2=n_all(ind_z==j);
     %------------------------
     [out1 out2]=dist_true(Z,ind1,ind2,bb); % length(ind1) * length(ind2)
     %------------------------
     temp=out1.*dist_z(ind1,ind2)+out2.*INF;
     dist_z(ind1,ind2)=temp;
     dist_z(ind2,ind1)=temp';
 end;
end;
%Dg=dist_z'+dist_z;
Dg=dist_z;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
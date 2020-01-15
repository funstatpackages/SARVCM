%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Spatial Autoregressive Partially Linear Varying Coefficient
% Model (SAR-CVM)
% Author: GuanNan Wang, JingRu Mu & Lily Wang
% Function: GWeights -- calculate geodesic distance and weight matrix
% Date: 12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dg Wg]=GWeights(V,T,bb,Z,~,~)
if nargin==4
    a=10;
    threshold=1e-3;
end

n=size(Z,1);
if (n+500<=1000)
    Vtarget=1000;
else
    Vtarget=n+500;
end
nadd=50;
if(n<=3000)
    ntrial=50;
elseif(n<10000)
    ntrial=15;
else
    ntrial=3;
end

% Find the rectangle which covers the entire domain;
bb_max=max(bb); bb_min=min(bb);
box=bb_max-bb_min;
box_h=box(1); box_v=box(2);
u_tmp=bb_min(1)+(bb_max(1)-bb_min(1))*rand(1000,1);
v_tmp=bb_min(2)+(bb_max(2)-bb_min(2))*rand(1000,1);
ind_tmp=sum(insideVT(V,T,u_tmp,v_tmp));

% calculate geo-distance;
if ind_tmp==1000
    n0=size(Z,1);
    Dg=reshape(sqrt(sum((kron(Z',ones(1,n0))-...
        repmat(Z',1,n0)).^2)),n0,n0);
    Wg=weights(Dg,1);
end
if ind_tmp<1000
    [Vnew,Tnew,~,threshold]=AddPoints(Vtarget,nadd,V,T,bb,Z);  
    [Dg0,D,M]=shortpath(Vnew,Tnew,Z);
    Wg=weights(Dg0,1);
    flag=0;
    for i=1:ntrial
        [Vnew,Tnew,step,threshold]=AddPoints(Vtarget,nadd,V,T,bb,Z);
        [Dg,D,M]=shortpath(Vnew,Tnew,Z);
        Dg_temp=[Dg0(:) Dg(:)]; 
        Dg=min(Dg_temp');
        Dg_old=Dg0;
        Dg0=reshape(Dg,n,n);  
        Wg_old=Wg;
        Wg=weights(Dg0,1);
        w_diff=sum(((Wg(:)-Wg_old(:)).^2))/(n^2);
        if(sqrt(w_diff)<1e-5)
            flag=1
            break;
        end    
    end
    Dg=Dg0;
end
end
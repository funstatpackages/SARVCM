%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Spatial Autoregressive Partially Linear Varying Coefficient
% Model (SAR-CVM)
% Example 2. Horseshoe Domain
% Author: GuanNan Wang, JingRu Mu & Lily Wang
% Function: DataGenerator -- generate sample points
% Date: 12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=500;
% sigma=1;
% alpha=0.5;
function [dat,Wg]=DataGenerator2(pop,n,sigma,alpha,bb,V,T,iter)
if n==200
    load('eg2_sam_ind_n200.mat');
elseif n==500
    load('eg2_sam_ind_n500.mat');
elseif n==1000
    load('eg2_sam_ind_n1000.mat');
end
ind=sort(sam_ind(1:n,iter));
sam=pop(ind,:);
Ystar=sam(1:n,1);
beta1=sam(1:n,2);
beta2=sam(1:n,3);
beta3=sam(1:n,4);
x1=sam(1:n,5);
x2=sam(1:n,6);
x3=sam(1:n,7);
u=sam(1:n,8);
v=sam(1:n,9);
Z=sam(1:n,8:9);

% Calculate the pairwise Euclidean distances and Weights
[Dg Wg]=GWeights(V,T,bb,Z);
A=diag(ones(n,1))-alpha*Wg;
Y=A\Ystar;
dat=[Y beta1 beta2 beta3 x1 x2 x3 u v];
end





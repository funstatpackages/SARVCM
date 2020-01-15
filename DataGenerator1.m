%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Spatial Autoregressive Partially Linear Varying Coefficient
% Model (SAR-CVM)
% Example 1. Rectangular Domain
% Author: GuanNan Wang, JingRu Mu & Lily Wang
% Function: DataGenerator -- generate sample points
% Date: 12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=500;
% sigma=1;
% alpha=0.5;
function [dat,D1]=DataGenerator1(n,sigma,alpha)
u=rand(n,1); v=rand(n,1); Z=[u v];
s=u.^2+v.^2;
beta1=sin(pi.*s);
beta2=cos(pi.*s);
beta3=exp(s);
beta4=ones(n,1)*(-1);
beta5=ones(n,1)*(1);
beta=[beta1 beta2 beta3 beta4 beta5];
X=normrnd(0,1,n,5);
eps=normrnd(0,sigma,n,1);
Ystar=sum(X.*beta,2)+eps;

% Calculate the true pairwise distances and Weights;
% (in this example, since the domain is rectangular, the tru distance
% should be Euclidean distance);
D1=reshape(sqrt(sum((kron(Z',ones(1,n))-repmat(Z',1,n)).^2)),n,n);
W=weights(D1,1);
A=diag(ones(n,1))-alpha*W;
Y=A\Ystar;
dat=[Y beta X u v];
end





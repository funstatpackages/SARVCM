%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Spatial Autoregressive Partially Linear Varying Coefficient
% Model (SAR-CVM)
% Example 2. Horseshoe Domain
% Author: GuanNan Wang, JingRu Mu & Lily Wang
% Function: Main Function
% Date: 12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
warning off;
format long g;

% Population
n1=180; n2=80;
load('eg2_pop.mat');
Ystar_pop=pop(:,1);
beta1_pop=pop(:,2);
beta2_pop=pop(:,3);
beta3_pop=pop(:,4);
beta_pop=[beta1_pop beta2_pop beta3_pop];
Xpop=pop(:,5:7);
Zpop=pop(:,8:9);

% Boundary and Triangulation (V1 T1, V2 T2 and V3 T3);
load('bnd_n60.mat');
load('eg2_V2.mat'); % Vertices;
load('eg2_T2.mat'); % Triangulation;
figure;
triplot(T,V(:,1),V(:,2),'k');

tmp=insideVT(V,T,Zpop(:,1),Zpop(:,2));
npop=size(pop,1);
ind_pop=1:npop;
ind_pop=ind_pop(tmp==1 & ~isnan(Ystar_pop));

% Bivariate Spline;
d=2; %d>r is the degree of spline functions over (V,T).
r=1; %r is the smoothness of the spline functions over (V,T).
H=smoothness(V,T,d,r); %H is the smoothness condition matrix.
Q2=qrH(H); %QR decomposition of the transpose of matrix H.
K=energy(V,T,d); %energy matrix;

% Simulation Parameters;
n=1000; % Sample size
np=size(beta_pop,2); % # Covariates;
alpha=0.5; 
sigma=1; % White noise;
ind_nl_true=1:2; ind_l_true=3; % Varying and non-varying coefficients;
nSIM=200; % # Replications;
penalty=true; % Penalized spline (true/false);

% For results;
alpha_mse_all=[];
beta_mse_all=[];
beta_mise_all=[];
ind_nl_all=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Main function;
tic
for iter=1:nSIM 
    %iter=1;
    rng(iter);
    iter
    [dat Wg]=DataGenerator2(pop,n,sigma,alpha,bb,V,T,iter);
    Y=dat(:,1);
    beta=dat(:,2:(np+1)); % True varying coefficient function;
    X=dat(:,(np+2):(2*np+1)); % Covariates;
    Z=dat(:,(2*np+2):(2*np+3)); % Location info;

    % estimation;
    ind_start=-6;
    ind_end=3;
    stepSize=0.5;
    index=ind_start:stepSize:ind_end;
    lambda=10.^index;

    [B B0 Ind]=Basis(Z,V,T,d);

    X_l=X(:,ind_l_true); X_nl=X(:,ind_nl_true);
    [theta_hat gamma_hat beta_hat...
        alpha_hat Ystar_hat Yhat mse mle aic bic]=...
        fitSVCAM(X_l,X_nl,Y,Wg,B,Q2,K,lambda,penalty);
    
    alpha_mse=(alpha-alpha_hat)^2;
    alpha_mse_all=[alpha_mse_all alpha_mse];

    beta_hat_all=zeros(n,np);
    beta_hat_all(:,ind_l_true)=repmat(theta_hat',n,1);
    beta_hat_all(:,ind_nl_true)=beta_hat;
    beta_mse=mean((beta-beta_hat_all).^2,1);
    beta_mse_all=[beta_mse_all; beta_mse];

    beta_pop_hat=beta_pred(Zpop,V,T,d,gamma_hat,theta_hat,ind_nl_true,...
        ind_l_true);
    beta_mise=nanmean((beta_pop_hat-beta_pop).^2,1);
    beta_mise_all=[beta_mise_all; beta_mise];

    display(['alpha_mse = ',num2str(alpha_mse)]);
    display(['beta_mise = ',num2str(beta_mise)]);
end
toc
alpha_mse=round(mean(alpha_mse_all),4)
beta_mise=round(mean(beta_mise_all,1),4)


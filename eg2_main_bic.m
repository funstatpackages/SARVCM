%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Spatial Autoregressive Partially Linear Varying Coefficient
% Model (SAR-CVM)
% Example 2. Horseshoe Domain
% Author: GuanNan Wang, JingRu Mu & Lily Wang
% Function: Main Function with Model Selection via BIC
% Date: 12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
warning off;
format long g;

path(path,'eg2_4Var');

% Population
load('pop_4Var.mat');
Ystar_pop=pop(:,1);
beta1_pop=pop(:,2);
beta2_pop=pop(:,3);
beta3_pop=pop(:,4);
beta4_pop=pop(:,5);
beta_pop=[beta1_pop beta2_pop beta3_pop beta4_pop];
Xpop=pop(:,6:9);
Zpop=pop(:,10:11);
n1=180;
n2=80;
np=size(beta_pop,2);
load('bnd_n60.mat');

% Boundary and Triangulation (V1 T1/V2 T2/V3 T3);
load('bnd_n60.mat');
load('eg2_V1.mat'); % Vertices;
load('eg2_T1.mat'); % Triangulation;
% figure;
% triplot(T,V(:,1),V(:,2),'k');

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
n=1000; % Sample size (500, 1000);
np=size(beta_pop,2); % # Covariates;
alpha=0.5; 
sigma=1; % White noise;
ind_nl_true=1:2; ind_l_true=3:4; % Varying and non-varying coefficients;
nSIM=200; % # Replications;
penalty=true; % Penalized spline (true/false);
criterion='bic';
alg=1; % Algorithm used to calculate the geodesic distance (1/2);

% For results;
alpha_mse_all=[];
beta_mse_all=[];
beta_mise_all=[];
ind_nl_all=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Main function;
for iter=1:nSIM
	% iter=80;
    rng(iter);
    iter
    % [dat Wg]=DataGenerator2(pop,n,sigma,alpha,bb,V,T,iter,alg);
    load(['dat_1000_' num2str(iter) '.mat']);
    load(['W_1000_' num2str(iter) '.mat']);
    Wg=W;
    Y=dat(:,1);
    beta=dat(:,2:(np+1)); % True varying coefficient function;
    X=dat(:,(np+2):(2*np+1)); % Covariates;
    Z=dat(:,(2*np+2):(2*np+3)); % Location info;
    
    tic
    % estimation;
    ind_start=-6;
    ind_end=3;
    stepSize=0.5;
    index=[ind_start:stepSize:ind_end]; 
    lambda=[10.^index];
    
    [B B0 Ind]=Basis(Z,V,T,d);

    [theta_hat gamma_hat beta_hat alpha_hat Ystar_hat Yhat ind_l ind_nl...
        mse mle aic bic]=SVCAM_sel(X,Y,Wg,B,Q2,K,lambda,penalty,criterion);

    indnl_all=zeros(np,1);
    indnl_all(ind_nl)=1;
    ind_nl_all=[ind_nl_all; indnl_all'];
    
    alpha_mse=(alpha-alpha_hat)^2;
    alpha_mse_all=[alpha_mse_all alpha_mse];

    beta_hat_all=zeros(n,np);
    beta_hat_all(:,ind_l)=repmat(theta_hat',n,1);
    beta_hat_all(:,ind_nl)=beta_hat;
    beta_mse=mean((beta-beta_hat_all).^2,1);
    beta_mse_all=[beta_mse_all; beta_mse];

    beta_pop_hat=beta_pred(Zpop,V,T,d,gamma_hat,theta_hat,ind_nl,ind_l);
    beta_mise=nanmean((beta_pop_hat-beta_pop).^2,1);
    beta_mise_all=[beta_mise_all; beta_mise];

    display(['ind_l = ',num2str(ind_l),' and ind_nl = ',num2str(ind_nl)]);
    display(['alpha_mse = ',num2str(alpha_mse)]);
    display(['beta_mise = ',num2str(beta_mise)]);
    toc
end

alpha_mse=round(mean(alpha_mse_all),4)
beta_mise=round(mean(beta_mise_all,1),4)
[unique_rows,~,ind]=unique(ind_nl_all,'rows');
unique_rows
counts=histc(ind,unique(ind))/nSIM


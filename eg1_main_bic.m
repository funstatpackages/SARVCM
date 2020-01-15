%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Spatially Varying Coefficient Models (SVCM)
% Example 1. Horse-Shoe Domain
% Author: GuanNan Wang, JingRu Mu & Lily Wang
% Date: 02012017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
warning off;
format long g;

% Population
n1=101; n2=101;
[beta_pop Xpop Zpop]=PopGenerator(n1,n2);

% Triangulation;
bb=[-0.01,-0.01;1.01,-0.01;1.01,1.01;-0.01,1.01;-0.01,-0.01];
load('eg1_V2.mat'); % Vertices; (V1/V2);
load('eg1_T2.mat'); % Triangulation; (V1/V2)

% Bivariate Spline;
d=2; %d>r is the degree of spline functions over (V,T).
r=1; %r is the smoothness of the spline functions over (V,T).
H=smoothness(V,T,d,r); %H is the smoothness condition matrix.
Q2=qrH(H); %QR decomposition of the transpose of matrix H.
K=energy(V,T,d); %energy matrix;

% Simulation Parameters;
n=1000; % Sample size (500/1000);
np=5; % # Covariates;
alpha=0.5;
sigma=1; % White noise;
ind_nl_true=1:3; ind_l_true=4:5; % Varying and non-varying coefficients;
nSIM=200; % # Replications;
penalty=true; % Penalized spline (true/false);
criterion='bic';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_mse_all=[];
beta_mse_all=[];
beta_mise_all=[];
ind_nl_all=[];

tic
for iter=1:nSIM
    %iter=1;
    rng(iter);
    iter
    % Generate sample data;
    [dat, W]=DataGenerator1(n,sigma,alpha);
    Y=dat(:,1); % Response;
    beta=dat(:,2:(np+1)); % True varying coefficient function;
    X=dat(:,(np+2):(2*np+1)); % Covariates;
    Z=dat(:,(2*np+2):(2*np+3)); % Location info;
    
    % Calculate geodesic distance and weights;
    [Dg, Wg]=GWeights(V,T,bb,Z);
    
    % estimation;
    ind_start=-6;
    ind_end=3;
    stepSize=0.5;
    index=ind_start:stepSize:ind_end; 
    lambda=10.^index;
    
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
end
toc

alpha_mse=round(mean(alpha_mse_all),4)
beta_mise=round(mean(beta_mise_all,1),4)
[unique_rows,~,ind]=unique(ind_nl_all,'rows');
unique_rows
counts=histc(ind,unique(ind))/nSIM

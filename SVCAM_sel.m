%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Spatial Autoregressive Partially Linear Varying Coefficient
% Model (SAR-CVM)
% Author: GuanNan Wang, JingRu Mu & Lily Wang
% Function: Model Selection for SARVCM
% Date: 12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta_hat gamma_hat beta_hat alpha_hat...
    Ystar_hat Yhat ind_l ind_nl mse mle aic bic]=...
    SVCAM_sel(X,Y,Wg,B,Q2,K,lambda,penalty,criterion)
n=length(Y);
np=size(X,2);
ind_l=[];
ind_nl=setdiff(1:np,ind_l);
X_l=X(:,ind_l);
X_nl=X(:,ind_nl);
[theta_hat gamma_hat beta_hat alpha_hat Ystar_hat Yhat...
    mse mle aic bic]=fitSVCAM(X_l,X_nl,Y,Wg,B,Q2,K,lambda,penalty);
if(strcmp(criterion,'aic'))
    abic_old=aic;
    abic_new=aic;
end;
if(strcmp(criterion,'bic'))
    abic_old=bic;
    abic_new=bic;
end;

jnl=0;
while(abic_new<=abic_old & length(ind_nl)>1)
    abic_old=abic_new;
    if(jnl~=0)
        ind_l=[ind_l ind_nl(jnl)];
        ind_nl=[ind_nl(1:(jnl-1)) ind_nl((jnl+1):length(ind_nl))];
    end;
    nnl=length(ind_nl);
    abic_all=[];
    for(inl=1:nnl)
        tmp_l=[ind_l ind_nl(inl)];
        tmp_nl=[ind_nl(1:(inl-1)) ind_nl((inl+1):nnl)];
        X_l=X(:,tmp_l);
        X_nl=X(:,tmp_nl);
        [theta_tmp gamma_tmp beta_tmp alpha_tmp Ystar_tmp Yhat_tmp...
            mse mle aic bic]=fitSVCAM(X_l,X_nl,Y,Wg,B,Q2,K,lambda,penalty);
        if(strcmp(criterion,'aic'))
            abic_all=[abic_all aic];
        end;
        if(strcmp(criterion,'bic'))
            abic_all=[abic_all bic];
        end;
    end;
    [abic_new jnl]=min(abic_all);
end;
ind_l=sort(ind_l);
ind_nl=sort(ind_nl);
X_l=X(:,ind_l);
X_nl=X(:,ind_nl);
[theta_hat gamma_hat beta_hat alpha_hat Ystar_hat Yhat...
    mse mle aic bic]=fitSVCAM(X_l,X_nl,Y,Wg,B,Q2,K,lambda,penalty);
end
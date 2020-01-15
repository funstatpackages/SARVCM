function [theta_hat gamma_hat beta_hat alpha_hat...
    Ystar_hat Yhat mse mle aic bic]=...
    SVCAM_est(X_l,X_nl,Y,Wg,B,Q2,K,lambda)
n=length(Y);
alpha_init=(1:99)./100;
nalpha=length(alpha_init);
    
mle_alpha=[];
mse_alpha=[];
det_alpha=[];
lam_alpha=[];
for(ia=1:nalpha)
	Ai=eye(n)-alpha_init(ia)*Wg;
	Ystar=Ai*Y;
	[theta gamma beta Yhat sse gcv df lam_c]=...
        plsfitGCV(B,Q2,K,X_l,X_nl,Ystar,lambda);
	mse=sse/n;
    detA=log(det(Ai));
    mle=detA-n*log(sqrt(mse));
    lam_alpha=[lam_alpha mean(lam_c)];
    mse_alpha=[mse_alpha mle];
    det_alpha=[det_alpha detA];
	mle_alpha=[mle_alpha mle];
end;
[mle j]=max(mle_alpha);
alpha_hat=alpha_init(j);

Ahat=eye(n)-alpha_hat*Wg;
Ystar=Ahat*Y;
lam_c=lam_alpha(j);
[theta_hat gamma_hat beta_hat Ystar_hat sse gcv df lamc]=...
    plsfitGCV(B,Q2,K,X_l,X_nl,Ystar,lam_c);
mse_star=sse/n;
detA=log(det(Ai));
mle=detA-n*log(sqrt(mse_star));
aic=df-mle;
bic=df*log(n)/2-mle;
Yhat=inv(Ahat)*Ystar_hat;
mse=mean((Y-Yhat).^2);
end

function [beta_pop Xpop Zpop]=PopGenerator(n1,n2)
u_min=0; u_max=1;
v_min=0; v_max=1;
step1=(u_max-u_min)/(n1-1);
step2=(v_max-v_min)/(n2-1);
u1=u_min:step1:u_max;
v1=v_min:step2:v_max;
u_vec=repmat(u1,n2,1);
v_vec=repmat(v1,1,n1);
Zpop=[u_vec(:) v_vec'];

s_pop=Zpop(:,1).^2+Zpop(:,2).^2;
beta1_pop=sin(pi.*s_pop);
beta2_pop=cos(pi.*s_pop);
beta3_pop=exp(s_pop);
beta4_pop=ones(n1*n1,1)*(-1);
beta5_pop=ones(n1*n1,1)*(1);
beta_pop=[beta1_pop beta2_pop beta3_pop...
    beta4_pop beta5_pop];
np=size(beta_pop,2);
% for(ip=1:np)
%     beta_tmp=beta_pop(:,ip);
% 	figure;
%     [C,h]=contourf(reshape(u_vec,n1,n2),...
%         reshape(v_vec,n1,n2),...
%         reshape(beta_tmp,n1,n2));
%     clabel(C,h);
% end;
Xpop=normrnd(0,1,n1*n2,np);
end
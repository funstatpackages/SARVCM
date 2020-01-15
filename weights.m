%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Spatial Autoregressive Partially Linear Varying Coefficient
% Model (SAR-CVM)
% Author: GuanNan Wang, JingRu Mu & Lily Wang
% Function: weights -- calculate weight matrix
% Date: 12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w=weights(dist,func,~,~)
if nargin==2
    a=10;
    threshold=1e-3;
end

n=size(dist,1);
if(func==1)
    temp1=sum(exp(-dist*a),2)-1;
    temp2=repmat(temp1,1,n);
    w=exp(-dist*a)./temp2;
    w=w-diag(diag(w));
end

if(func==2)
    temp1=1./dist;
    temp1(isinf(temp1))=0;
    temp2=sum(temp1,2);
    temp2=repmat(temp2,1,n);
    w=temp1./temp2;
end

w(w(:)<threshold)=0;
temp3=sum(w,2);
temp4=repmat(temp3,1,n);
w=w./temp4;
w=sparse(w);
end
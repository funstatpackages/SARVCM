%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: qrH.m
% Function: Q2=qrH(H)
% This matlab program implements QR decomposition of the transpose of 
% matrix H.
% 
% Input:
% H: smoothness condition matrix
%
% Output:
% Q2: H^{T}=QR=[Q1 Q2][R1 R2]^{T}
%
% This matlab program is copyrighted @2015 by Lily Wang and Guannan Wang.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q2=qrH(H)
    [QH RH]=qr(H');
    [pH qH]=size(H);
    r2=0;
    for j=1:qH
        r2=r2+(sum(RH(j,:)==0)==pH);
    end
    r1=qH-r2;
    Q2=QH(:,((r1+1):qH));
end
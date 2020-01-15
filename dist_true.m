function [out1 out2]=dist_true(Z,ind1,ind2,bb)
XY2=[];
n_bnd=size(bb,1);
if (bb(1,1)~=bb(n_bnd,1) & bb(1,2)~=bb(n_bnd,2))
    bb=[bb; bb(1,:)];
    n_bnd=n_bnd+1;
end;

for(i=1:(n_bnd-1))
    j=i+1;
    temp1=[bb(i,:) bb(j,:)];
    XY2=[XY2; temp1];
end

n_1=length(ind1);
n_2=length(ind2);


ind1=reshape(ind1,n_1,1);
ind2=reshape(ind2,1,n_2);
ind1=repmat(ind1,n_2,1);
ind2=repmat(ind2,n_1,1);
ind2=ind2(:);


XY1=[Z(ind1,:) Z(ind2,:)];


%-----------------------------------------------------------
validateattributes(XY1,{'numeric'},{'2d','finite'});
validateattributes(XY2,{'numeric'},{'2d','finite'});
[n_rows_1,n_cols_1] = size(XY1);
[n_rows_2,n_cols_2] = size(XY2);
if n_cols_1 ~= 4 || n_cols_2 ~= 4
    error('Arguments must be a Nx4 matrices.');
end
%-------------------------------------------------------------------------------
X1 = repmat(XY1(:,1),1,n_rows_2);
X2 = repmat(XY1(:,3),1,n_rows_2);
XY2 = XY2';
X3 = repmat(XY2(1,:),n_rows_1,1);
X4 = repmat(XY2(3,:),n_rows_1,1);
X2_X1 = (X2-X1); 
X4_X3 = (X4-X3); 
X1_X3 = (X1-X3);
%X2_X1 = (X2-X1);
clear X1 X2 X3 X4;
Y1 = repmat(XY1(:,2),1,n_rows_2);
Y2 = repmat(XY1(:,4),1,n_rows_2);
clear XY1
Y3 = repmat(XY2(2,:),n_rows_1,1);
Y4 = repmat(XY2(4,:),n_rows_1,1);
clear XY2
Y1_Y3 = (Y1-Y3);
Y4_Y3 = (Y4-Y3);
Y2_Y1 = (Y2-Y1);
clear Y1 Y2 Y3 Y4;
numerator_a = X4_X3 .* Y1_Y3 - Y4_Y3 .* X1_X3;
numerator_b = X2_X1 .* Y1_Y3 - Y2_Y1 .* X1_X3;
denominator = Y4_Y3 .* X2_X1 - X4_X3 .* Y2_Y1;
clear X4_X3 Y1_Y3 Y4_Y3 X1_X3 X2_X1 Y2_Y1;
u_a = numerator_a ./ denominator;
clear numberator_a;
temp1=(u_a >= 0) & (u_a <= 1);
clear u_a;
u_b = numerator_b ./ denominator;
INT_B = temp1 & (u_b >= 0) & (u_b <= 1);
clear u_b;
clear denominator;

temp = sum(INT_B,2);
out1=zeros(length(temp),1);
out1(temp==0)=1;
out1=reshape(out1,n_1,n_2);

out2=ones(length(temp),1);
out2(temp==0)=0;
out2=reshape(out2,n_1,n_2);
end
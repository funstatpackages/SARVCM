%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Spatially Varying Coefficient Models (SVCM)
% Example 2. Rectangular Domain
% Author: GuanNan Wang, JingRu Mu & Lily Wang
% Date: 02162017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
warning off;
format long g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate sample points;
load('sam_horseshoe.mat');
n=size(Z,1);

% make the boundary finer
%filename='boundary_r.txt';
%bb=textread(filename,'%f');
%bb=reshape(bb.',2,[]).';
load('bnd_n60.mat');

% Triangulation;
n_bb=size(bb,1)-1;
bb_start=1:(n_bb-1);
bb_end=2:n_bb;
bb_edge=[bb_start' bb_end'];
bb_edge=[bb_edge;[n_bb 1]];

% all points (sample points + boundary);
bb_t=bb(1:n_bb,:);
V=Z;
C=bb_edge+size(V,1);
V=[V; bb_t];
t_re=delaunayTriangulation(V(:,1),V(:,2),C);
io=t_re.isInterior();
T=t_re(io==1,:);
V=t_re.Points;
triplot(T,V(:,1),V(:,2),'k');
hold on;
plot(Z(:,1),Z(:,2),'r*');
hold off;

% euclidian distance;
D1=reshape(sqrt(sum((kron(Z',ones(1,n))-repmat(Z',1,n)).^2)),n,n);

D_max=[];
D_pmax=[];
W_max=[];
W_pmax=[];
D_all=[];
for(iter=1:30)
iter
% add random points and create new triangulation;
ts=1; % parameters to determine threshold;
pn=50;
[Vnew,Tnew,step,threshold]=...
    AddPoints(1000,50,V,T,bb,Z,ts,pn);

% calculate geo-distance;
[D,DG,M]=shortpath(Vnew,Tnew,Z);
D_all=[D_all D(:)];

% compare geo-distance with distance;
Dmax=max(max(abs(D-D1)));
D_max=[D_max Dmax];
temp1=max(D1)';
temp2=repmat(temp1,1,n);
Dp=abs(D-D1)./temp2;
Dpmax=max(max(Dp));
D_pmax=[D_pmax Dpmax];

% calculate weights based on distance;
func=1; % weight function: (1) exp. (2) reciprocal;
Dw_e=weights(D1,func,n); 
Dw_g=weights(D,func,n);
Dwmax=max(max(abs(Dw_e-Dw_g)));
W_max=[W_max Dwmax];
temp1=max(Dw_e')';
temp2=repmat(temp1,1,n);
Dwp=abs(Dw_g-Dw_e)./temp2;
Dwpmax=max(max(Dwp));
W_pmax=[W_pmax Dwpmax];
end;
Dnew=min(D_all');
Dnew=reshape(Dnew,n,n);
Dnewmax=max(max(abs(Dnew-D1)));
temp1=max(D1)';
temp2=repmat(temp1,1,n);
Dnewp=abs(Dnew-D1)./temp2;
Dnewpmax=max(max(Dnewp));
[Dnewmax Dnewpmax]

% calculate weights based on distance;
func=1; % weight function: (1) exp. (2) reciprocal;
Dw_new=weights(Dnew,func,n);
Dwnewmax=max(max(abs(Dw_new-Dw_e)));
temp1=max(Dw_new')';
temp2=repmat(temp1,1,n);
Dwp_new=abs(Dw_new-Dw_e)./temp2;
Dwpnewmax=max(max(Dwp_new));
[Dwnewmax Dwpnewmax]
%round([quantile(D_max,[0,0.25,0.5,0.75,1]) mean(D_max) std(D_max)],4)
%round([quantile(D_pmax,[0,0.25,0.5,0.75,1]) mean(D_pmax) std(D_pmax)],4)
%round([quantile(W_max,[0,0.25,0.5,0.75,1]) mean(W_max) std(W_max)],4)
%round([quantile(W_pmax,[0,0.25,0.5,0.75,1]) mean(W_pmax) std(W_pmax)],4)

% e.g. check distance via plot between 36th and 4th point
% [dss,p1,p2]=graphshortestpath(sparse(DG),36,4);
% triplot(Tnew,Vnew(:,1),Vnew(:,2))
% hold on
% check=[36 4]
% plot(Vnew(p1,1),Vnew(p1,2),'linewidth',2)
% hold off
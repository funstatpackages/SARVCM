%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project: Spatial Autoregressive Partially Linear Varying Coefficient
% Model (SAR-CVM)
% Author: GuanNan Wang, JingRu Mu & Lily Wang
% Function: AddPoints -- add random points to create traingle meshes
% Input Arguments:
% (1) Vtarget: # points we used to create triangule meshes
% (2) nadd: # points added in each step
% (3) V: Vertice
% (4) T: Triangulation
% (5) bb: Boundary
% (6) Z: Location info
% (7) ts: 
% (8) pn: # points assumed to be in the box domain
% Date: 12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vnew,Tnew,step,threshold]=...
    AddPoints(Vtarget,nadd,V,T,bb,Z,~,~)

if nargin==6
    pn=100; ts=1;
end

% Find the rectangle which covers the entire domain;
bb_max=max(bb); bb_min=min(bb);
box=bb_max-bb_min;
box_h=box(1); box_v=box(2);
threshold=sqrt(box_h*box_v/(pn^2))/ts;

% edges for boundary;
n_bb=size(bb,1);
if sum(bb(1,:)==bb(n_bb,:))==2
    n_bb=n_bb-1;
end
bb_start=1:(n_bb-1);
bb_end=2:n_bb;
bb_edge=[bb_start' bb_end'];
bb_edge=[bb_edge;[n_bb 1]];

% all points (sample points + boundary);
bb_t=bb(1:n_bb,:);
Vnew=Z;
C=bb_edge+size(Vnew,1);
Vnew=[Vnew; bb_t];

% calculate the distance between the sample points;
n0=size(Z,1);
DZ0=reshape(sqrt(sum((kron(Z',ones(1,n0))-...
    repmat(Z',1,n0)).^2)),n0,n0);

zmin=min(DZ0(DZ0>0));
zmax=max(DZ0(DZ0>0));

if(threshold/zmax<1)
    if(threshold/zmin>1)
        threshold=zmin;
    end;
    step=0;
    while(size(Vnew,1)<Vtarget)
        step=step+1;
        Zr1=unifrnd(bb_min(1),bb_max(1),nadd,1);
        Zr2=unifrnd(bb_min(2),bb_max(2),nadd,1);
        Zr=[Zr1 Zr2];
        [ind_temp,~]=inside(V,T,Zr(:,1),Zr(:,2));
        Zr=Zr(ind_temp,:);
        Ztemp=[Vnew;Zr];
        nz0=size(Vnew,1);
        nzt=size(Ztemp,1);
    
        DZ=reshape(sqrt(sum((kron(Ztemp',ones(1,nzt))-...
            repmat(Ztemp',1,nzt)).^2)),nzt,nzt);
        
        ind_temp=(nz0+1):nzt;
        ind_re=ind_temp(sum(DZ(ind_temp,1:nz0)<threshold,2)>0);
        Ztemp(ind_re,:)=[];
        Vnew=Ztemp;
    end;
end;
% plot(bb(:,1),bb(:,2),'k');
% hold on;
% plot(Vnew(:,1),Vnew(:,2),'bo');
% plot(Z(:,1),Z(:,2),'r*');
% hold off;
%t_re=delaunayTriangulation(Vnew(:,1),Vnew(:,2));
%triplot(t_re);
t_re=delaunayTriangulation(Vnew(:,1),Vnew(:,2),C);
io=t_re.isInterior();
Tnew=t_re(io==1,:);
Vnew=t_re.Points;
% figure;
% triplot(Tnew,Vnew(:,1),Vnew(:,2))
% hold on
% plot(Z(:,1),Z(:,2),'ro')
% hold off
end
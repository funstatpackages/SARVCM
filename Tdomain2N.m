function [VT,T] = Tdomain2N(Pt,n,H);
% Empty hole
if nargin == 2
	H{1} = [];
end

[X,Y] = gridpoly(Pt,n,H);
T = delaunay(X,Y);
VT = [X,Y];
T = del_tri(Pt,VT,T,1);
% Delete triangles within the holes
if ~isempty(H)
   for i=1:numel(H)
   T = del_tri(H{i},VT,T,-1);
   end
end; 
% Look for degenrate triangles
tol=10^(-18);
area=[];
for i=1:size(T(:,1))
area=[area triarea(VT(T(i,1),:),VT(T(i,2),:),VT(T(i,3),:))];
end;
dT=find(area<tol);
T(dT,:)=[];

 
function NewT = del_tri(Pt,V,T,order);
%        NewT = del_tri(Pt,V,T,order)
% This function deletes from the triangulation (V,T) all triangles
% which are outside the polygon with vertices Pt
  n = size(T);
  n = n(1);
  i = 1;
  if nargin < 4
     order = polyord(Pt);
  end;
  while i <= n
     C = (V(T(i,1),:)+V(T(i,2),:)+V(T(i,3),:))/3;
    if (~inpoly(C,Pt)&order > 0) | (inpoly(C,Pt)& order < 0)
      T = [T(1:(i-1),:);T((i+1):n,:)];
      n = n - 1;      
    else
      i = i + 1;
    end;
  end;
  NewT = T;    
    
  
 function [X,Y] = gridpoly(Pt,n,H);
%        [X,Y,bdr] = gridpoly(Pt,n)
% This function returns grid points inside and on the boundary of
% the polygon with extreme points Pt and interior holes defined by H
  a = min(Pt(:,1));
  b = max(Pt(:,1));
  c = min(Pt(:,2));
  d = max(Pt(:,2));
  m = size(Pt,1);
  delta_x = (b-a)/n;
  delta_y = (d-c)/n;
  delta = norm([delta_x,delta_y]);
  x = a:delta_x:b;
  y = c:delta_y:d;

  delta=min(delta_x,delta_y);
  x=a:delta:b;          % make the inner point denser by times 0.5
  y=c:delta:d;
  [X,Y] = meshgrid(x,y);
  X = X(:);
  Y = Y(:);
%
% We now delete all points which are not inside the polygon
%
  k = 1; 
  while k <= length(X)
      % Check for points in each hole
      inHole = [];
      for i=1:numel(H)
            inHole = [inHole inpoly([X(k),Y(k)],H{i})];    
        end
        if sum(inHole) >= 1
             inHoleTest = 1;
        else
             inHoleTest = 0;
        end;   
%    if ~inpoly([X(k),Y(k)],Pt)| inpoly([X(k),Y(k)],H)     
     if ~inpoly([X(k),Y(k)],Pt)|inHoleTest
      X = [X(1:(k-1));X((k+1):length(X))];
      Y = [Y(1:(k-1));Y((k+1):length(Y))];
    else
      k = k + 1;
    end; 
 end;
%
% Now we must delete all points which are too close to the bdr
%
  m = length(X);
  p = size(Pt,1);
%  % q = size(H,1);
%  % q=0;
  v = [Pt([2:p,1],1)-Pt(1:p,1),Pt([2:p,1],2)-Pt(1:p,2)];
  L = sqrt(v(:,1).^2 + v(:,2).^2);
  v = [v(:,1)./L,v(:,2)./L];
  q = []; % q vector for holes
  % For each hole
  for i=1:numel(H)
    q = [q size(H{i},1)];
    if q(i) > 0
       w{i} = [H{i}([2:q(i),1],1)-H{i}(1:q(i),1),H{i}([2:q(i),1],2)-H{i}(1:q(i),2)];
       LL{i} = sqrt(w{i}(:,1).^2 + w{i}(:,2).^2);
       w{i} = [w{i}(:,1)./LL{i},w{i}(:,2)./LL{i}];
    end;
  end;
  newX = [];
  newY = [];
  for k = 1:m;
     IP = (X(k)-Pt(:,1)).*v(:,1) + (Y(k)-Pt(:,2)).*v(:,2);
     Q = [Pt(:,1)+IP.*v(:,1),Pt(:,2)+IP.*v(:,2)];
     Indx1 = find(IP < 0);
     Indx2 = find(IP > L);
     Q(Indx1,:) = Pt(Indx1,:);
     Q(Indx2,:) = Pt(rem(Indx2,p)+1,:);
     D = sqrt((X(k)-Q(:,1)).^2 + (Y(k)-Q(:,2)).^2);
  % For each hole   
   dH=[];
  for i=1:numel(H)
     if q(i) > 0
       IP2 = (X(k)-H{i}(:,1)).*w{i}(:,1) + (Y(k)-H{i}(:,2)).*w{i}(:,2);
       Q2 = [H{i}(:,1)+IP2.*w{i}(:,1),H{i}(:,2)+IP2.*w{i}(:,2)];
       Indx1 = find(IP2 < 0);
       Indx2 = find(IP2 > LL{i});
       Q2(Indx1,:) = H{i}(Indx1,:);
       Q2(Indx2,:) = H{i}(rem(Indx2,q(i))+1,:);
       D2 = sqrt((X(k)-Q2(:,1)).^2 + (Y(k)-Q2(:,2)).^2);
     else
       D2 = inf;
     end;
     dH=[dH;D2];
  end;
     if min([D;dH]) > delta*1/3 %delta/2 is the thresh hold
      newX = [newX;X(k)];
      newY = [newY;Y(k)];
    end;
  end;
  X = newX;
  Y = newY;

% Next, we add on bdr points
% delta is the largest distance allowed between two point on the boundary;
 newPts=[];
 k=1;
 numofpt=0; %length of newPts;
 while k<=p;
    j=rem(k,p)+1; 
    newPts=[newPts;Pt(k,:)]; %adds the current pt;
    numofpt=numofpt+1;
    D=norm(Pt(j,:)-Pt(k,:)); %distance between two adjacent pts;
    numpts=floor(D/(delta)-0.01);     %number of points should be added;
    if numpts>0
        for i=1:numpts
            newPt=(i*Pt(j,:)+(numpts+1-i)*Pt(k,:))/(numpts+1);
            newPts=[newPts;newPt];
            numofpt=numofpt+1;
        end
    end
    k=k+1;
 end
 X=[X;newPts(:,1)];
 Y=[Y;newPts(:,2)];
 % For each hole
 for nH=1:numel(H)
    if q(nH)>0
        newPts=[];
        k=1;
        numofpt=0; %length of newPts;
        while k<=q(nH);
        j=rem(k,q(nH))+1; 
        newPts=[newPts;H{nH}(k,:)];%adds the current pt;
        numofpt=numofpt+1;
        D=norm(H{nH}(j,:)-H{nH}(k,:));%distance between two adjacent pts;
        numpts=floor(D/(delta)-0.01);%number of points should be added;
            if numpts>0
                for i=1:numpts
                    newPt=(i*H{nH}(j,:)+(numpts+1-i)*H{nH}(k,:))/(numpts+1);
                    newPts=[newPts;newPt];
                    numofpt=numofpt+1;
                end
            end
        k=k+1;
        end
        X=[X;newPts(:,1)];
        Y=[Y;newPts(:,2)];
    end
 end   

     
 


  function boolean = inpoly(Point,Poly);
%        boolean = inpoly(Point,Poly)
% Poly is an N X 2 matrix whose ith row are the coordinates of the ith
% vertex of a possibly non-convex polygon. If Point is inside the Poly
% then boolean = 1 else boolean = 0.
% This matlab program is copyrighted @2001 by Ming-Jun Lai and Paul Wenston
% through University of Georgia Research Foundation, Inc..
  N = size(Poly);
  N = N(1);
  s = 0;
  for i = 1:N
   j = rem(i,N)+1;
   s = s + rotation(Point,Poly(i,:),Poly(j,:));
  end;
  if abs(s) > 1
    boolean = 1;
  else
    boolean = 0;
  end; 
  
  function theta = rotation(P1,P2,P3);
%        theta = rotation(P1,P2,P3)
% theta returns the angle between the rotation from the line joining P1
% and P2 to the line joining P1 and P3.
 v1 = [P2(1)-P1(1),P2(2)-P1(2)]';
 v2 = [P3(1)-P1(1),P3(2)-P1(2)]';
 if det([v1,v2]) > 0
   s = 1;
 else 
   s = -1;
 end;
 a = sqrt((P2(1)-P1(1))^2 + (P2(2)-P1(2))^2);
 b = sqrt((P3(1)-P1(1))^2 + (P3(2)-P1(2))^2);
 c = sqrt((P3(1)-P2(1))^2 + (P3(2)-P2(2))^2);
 theta = real(s*acos((a^2 + b^2 - c^2)/(2*a*b)));
 
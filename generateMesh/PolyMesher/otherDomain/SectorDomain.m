%------------------------------ PolyTop ----------------------------------%
% Ref: A Pereira, C Talischi, GH Paulino, IFM Menezes, MS Carvalho        %
%      "Fluid Flow Topology Optimization in PolyTop: Stability and        %
%      Computational Implementation", Struct Multidisc Optim,             %
%      DOI 10.1007/s00158-014-1182-z                                      %
%-------------------------------------------------------------------------%
function [x] = SectorDomain(Demand,Arg)
  r1=4; r2=8;
  angle = 45; % in degrees
  % First sector
  xcs1 = 0; ycs1 = 0; a1s1 = (90-angle)*pi/180;
  a2s1 = pi/2;
  % Second sector
  xcs2 = (r1+r2)*cos(a1s1); ycs2 = (r1+r2)*sin(a1s1);
  a1s2 = pi + a1s1; a2s2 = pi + a2s1;
  Coord = [xcs1+r1*cos(a1s1), ycs1+r1*sin(a1s1); ...
           xcs1+r2*cos(a1s1), ycs1+r2*sin(a1s1); ...
           xcs1+r1*cos(a2s1), ycs1+r1*sin(a2s1); ...
           xcs1+r2*cos(a2s1), ycs1+r2*sin(a2s1); ...
           xcs2+r1*cos(a1s2), ycs2+r1*sin(a1s2); ...
           xcs2+r2*cos(a1s2), ycs2+r2*sin(a1s2); ...
           xcs2+r1*cos(a2s2), ycs2+r1*sin(a2s2); ...
           xcs2+r2*cos(a2s2), ycs2+r2*sin(a2s2)];
  BdBox = [min(Coord(:,1)),max(Coord(:,1)), ...
           min(Coord(:,2)),max(Coord(:,2))];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox,xcs1,ycs1,r1,r2,a1s1,a2s1,xcs2,ycs2,a1s2,a2s2);
    case('BC');    x = BndryCnds(Arg{:},BdBox,xcs1,ycs1,r1,r2,a1s1,a2s1,xcs2,ycs2,a1s2,a2s2);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(r1,r2,a1s1);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox,xcs1,ycs1,r1,r2,a1s1,a2s1,xcs2,ycs2,a1s2,a2s2)
  d1 = DistFncSector1(P,xcs1,ycs1,r1,r2,a1s1,a2s1);
  d2 = DistFncSector1(P,xcs2,ycs2,r1,r2,a1s2,a2s2);
  Dist = dUnion(d1,d2);
  d3 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
  Dist = dDiff(d3,Dist);
function Dist = DistFncSector1(P,xc,yc,r1,r2,a1,a2)
  d1 = dCircle(P,xc,yc,r2);
  d2 = dCircle(P,xc,yc,r1);
  d3 = dLine(P,xc+r1*cos(a2),yc+r1*sin(a2),xc,yc);
  d4 = dLine(P,xc,yc,xc+r1*cos(a1),yc+r1*sin(a1));
  Dist = dIntersect(dIntersect(dDiff(d1,d2),d3),d4);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox,xcs1,ycs1,r1,r2,a1s1,a2s1,xcs2,ycs2,a1s2,a2s2)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  I = PolyBoundary(Node,Element);
  NNodeBC = size(I,1);
  NodeBC = zeros(2*NNodeBC,3); 
  NodeBC(:,1) = I(reshape([1:NNodeBC;1:NNodeBC],2*NNodeBC,1));
  NodeBC(:,2) = mod(2:2*NNodeBC+1,2)+1;
  tmp = (r1+r2)*sin(a1s1);
  for i=1:size(I,1)
    if ((abs(Node(I(i),1)-BdBox(1))<eps))
      if (Node(I(i),2)>((r1+1/6*(r2-r1))+eps)) && (Node(I(i),2)<((r1+2/6*(r2-r1))-eps))
        r = Node(I(i),2) - (r1+1/4*(r2-r1));
        NodeBC(2*i-1,:) = [I(i),1,1-r^2/((r2-r1)*1/12)^2];
      elseif (Node(I(i),2)>((r1+4/6*(r2-r1))+eps)) && (Node(I(i),2)<((r1+5/6*(r2-r1))-eps))
        r = Node(I(i),2) - (r1+3/4*(r2-r1));
        NodeBC(2*i-1,:) = [I(i),1,1-r^2/((r2-r1)*1/12)^2];
      end
    elseif ((abs(Node(I(i),1)-(r1+r2)*cos(a1s1))<eps))
      if (Node(I(i),2)<((tmp-(r1+1/6*(r2-r1)))-eps)) && (Node(I(i),2)>((tmp-(r1+2/6*(r2-r1)))+eps))
        r = Node(I(i),2) - (tmp-(r1+1/4*(r2-r1)));
        NodeBC(2*i-1,:) = [I(i),1,1-r^2/((r2-r1)*1/12)^2];
      elseif (Node(I(i),2)<((tmp-(r1+4/6*(r2-r1)))-eps)) && (Node(I(i),2)>((tmp-(r1+5/6*(r2-r1)))+eps))
        r = Node(I(i),2) - (tmp-(r1+3/4*(r2-r1)));
        NodeBC(2*i-1,:) = [I(i),1,1-r^2/((r2-r1)*1/12)^2];
      end
    end
  end
  x = {NodeBC,[]};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(r1,r2,a1)
  xc = (r1+r2)*cos(a1); yc = (r1+r2)*sin(a1);
  PFix = [0,r1+1/6*(r2-r1);
          0,r1+1/3*(r2-r1);
          0,r1+2/3*(r2-r1);
          0,r1+5/6*(r2-r1);
          xc,yc-(r1+1/6*(r2-r1));
          xc,yc-(r1+1/3*(r2-r1));
          xc,yc-(r1+2/3*(r2-r1));
          xc,yc-(r1+5/6*(r2-r1))];
%-------------------------------------------------------------------------%
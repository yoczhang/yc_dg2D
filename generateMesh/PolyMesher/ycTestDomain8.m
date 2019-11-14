%------------------------------ PolyTop ----------------------------------%
% Ref: A Pereira, C Talischi, GH Paulino, IFM Menezes, MS Carvalho        %
%      "Fluid Flow Topology Optimization in PolyTop: Stability and        %
%      Computational Implementation", Struct Multidisc Optim,             %
%      DOI 10.1007/s00158-014-1182-z                                      %
%-------------------------------------------------------------------------%
function [x] = ycTestDomain8(Demand,Arg) % from SectorDomain
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
    case('BC');    x = BndryCnds(Arg{:},BdBox);
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
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  LeftEdgeNodes = find(abs(Node(:,1)-BdBox(1))<eps);
  LeftUpperNode = find(abs(Node(:,1)-BdBox(1))<eps & ...
                       abs(Node(:,2)-BdBox(4))<eps);
  RigthBottomNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
                         abs(Node(:,2)-BdBox(3))<eps);
  FixedNodes = [LeftEdgeNodes; RigthBottomNode];
  Supp = zeros(length(FixedNodes),3);
  Supp(:,1)=FixedNodes; Supp(1:end-1,2)=1; Supp(end,3)=1;
  Load = [LeftUpperNode,0,-0.5];
  x = {Supp,Load};
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
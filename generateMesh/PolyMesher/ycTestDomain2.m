%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = ycTestDomain2(Demand,Arg) % from MbbDomain
  BdBox = [0, 2, 0, 1.25];
  BdBox2 = [0., 0.8, 0., 0.2];
  BdBox3 = [0.8, 1.3, 0., 0.5];
  BdBox4 = [1.3, 2., 0., 1.];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox,BdBox2,BdBox3,BdBox4);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox,BdBox2,BdBox3,BdBox4);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox,BdBox2,BdBox3,BdBox4)
  %--- part 1
  d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
  d2 = dRectangle(P,BdBox2(1),BdBox2(2),BdBox2(3),BdBox2(4));
  d3 = dRectangle(P,BdBox3(1),BdBox3(2),BdBox3(3),BdBox3(4));
  d4 = dRectangle(P,BdBox4(1),BdBox4(2),BdBox4(3),BdBox4(4));
  d_diff = dUnion(d4,dUnion(d2,d3));
  Dist = dDiff(d1,d_diff);

%   %--- part 2, --- failed
%   xy_1 = [0,0.2];
%   xy_2 = [0.8,0.2];
%   xy_3 = [0.8,0.5];
%   xy_4 = [1.3,0.5];
%   xy_5 = [1.3,1];
%   xy_6 = [2,1];
%   xy_7 = [2,1.25];
%   xy_8 = [0,1.25];
%   d1 = dLine(P,xy_1(1),xy_1(2),xy_2(1),xy_2(2));
%   d2 = dLine(P,xy_2(1),xy_2(2),xy_3(1),xy_3(2));
%   d3 = dLine(P,xy_3(1),xy_3(2),xy_4(1),xy_4(2));
%   d4 = dLine(P,xy_4(1),xy_4(2),xy_5(1),xy_5(2));
%   d5 = dLine(P,xy_5(1),xy_5(2),xy_6(1),xy_6(2));
%   d6 = dLine(P,xy_6(1),xy_6(2),xy_7(1),xy_7(2));
%   d7 = dLine(P,xy_7(1),xy_7(2),xy_8(1),xy_8(2));
%   d8 = dLine(P,xy_8(1),xy_8(2),xy_1(1),xy_1(2));
%   Dist = dIntersect(d8,dIntersect(d7,dIntersect(d6,dIntersect(d5,dIntersect(d4,dIntersect(d3,dIntersect(d2,d1)))))));
  
%   %--- part 3, --- succeed
%   xy_1 = [0,0];
%   xy_2 = [2,0];
%   xy_3 = [2,1.25];
%   xy_4 = [0,1.25];
%   d1 = dLine(P,xy_1(1),xy_1(2),xy_2(1),xy_2(2));
%   d2 = dLine(P,xy_2(1),xy_2(2),xy_3(1),xy_3(2));
%   d3 = dLine(P,xy_3(1),xy_3(2),xy_4(1),xy_4(2));
%   d4 = dLine(P,xy_4(1),xy_4(2),xy_1(1),xy_1(2));
%   Dist = dIntersect(d4,dIntersect(d3,dIntersect(d2,d1)));

%   %--- part 4,--- failed
%   xy_1 = [0,0];
%   xy_2 = [1.8,0];
%   xy_3 = [1.8,0.2];
%   xy_4 = [2,0.2];
%   xy_5 = [2,1.25];
%   xy_6 = [0,1.25];
%   d1 = dLine(P,xy_1(1),xy_1(2),xy_2(1),xy_2(2));
%   d2 = dLine(P,xy_2(1),xy_2(2),xy_3(1),xy_3(2));
%   d3 = dLine(P,xy_3(1),xy_3(2),xy_4(1),xy_4(2));
%   d4 = dLine(P,xy_4(1),xy_4(2),xy_5(1),xy_5(2));
%   d5 = dLine(P,xy_5(1),xy_5(2),xy_6(1),xy_6(2));
%   d6 = dLine(P,xy_6(1),xy_6(2),xy_1(1),xy_1(2));
%   Dist = dIntersect(d6,dIntersect(d5,dIntersect(d4,dIntersect(d3,dIntersect(d2,d1)))));

%   %--- part 5,--- succeed
%   xy_1 = [0,0];
%   xy_2 = [1.8,0];
%   xy_3 = [2,0.2];
%   xy_4 = [2,1.25];
%   xy_5 = [0,1.25];
%   d1 = dLine(P,xy_1(1),xy_1(2),xy_2(1),xy_2(2));
%   d2 = dLine(P,xy_2(1),xy_2(2),xy_3(1),xy_3(2));
%   d3 = dLine(P,xy_3(1),xy_3(2),xy_4(1),xy_4(2));
%   d4 = dLine(P,xy_4(1),xy_4(2),xy_5(1),xy_5(2));
%   d5 = dLine(P,xy_5(1),xy_5(2),xy_1(1),xy_1(2));
%   Dist = dIntersect(d5,dIntersect(d4,dIntersect(d3,dIntersect(d2,d1))));
  
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
function [PFix] = FixedPoints(BdBox,BdBox2,BdBox3,BdBox4)

  %--- part1, ---
  xx1 = (BdBox2(1): (BdBox2(2)-BdBox2(1))/10 : BdBox2(2))';
  yy1 = BdBox2(4) * ones(length(xx1),1);
  
  yy2 = (BdBox2(4): (BdBox3(4)-BdBox2(4))/10 : BdBox3(4))';
  xx2 = BdBox2(2) * ones(length(yy2),1);
  
  xx3 = (BdBox3(1): (BdBox3(2)-BdBox3(1))/10 : BdBox3(2))';
  yy3 = BdBox3(4) * ones(length(xx3),1);
  
  yy4 = (BdBox3(4): (BdBox4(4)-BdBox3(4))/10 : BdBox4(4))';
  xx4 = BdBox3(2) * ones(length(yy4),1);
  
  xx5 = (BdBox4(1): (BdBox4(2)-BdBox4(1))/10 : BdBox4(2))';
  yy5 = BdBox4(4) * ones(length(xx5),1);
  
  xx = [xx1;xx2;xx3;xx4;xx5];
  yy = [yy1;yy2;yy3;yy4;yy5];
  PFix = [xx, yy];
  PFix = unique(PFix,'rows','stable');
  
%-------------------------------------------------------------------------%
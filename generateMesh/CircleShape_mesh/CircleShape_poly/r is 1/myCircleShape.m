function [x] = myCircleShape(Demand,Arg)
  BdBox = [-1, 1, -1, 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dRectangle(P,-1, 1, -1, 1);
  d2 = dCircle(P,0,0,1);
  Dist = dIntersect(d1,d2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  Supp = find(abs(sqrt(Node(:,1).^2+Node(:,2).^2)-1)<1e-2);
  Load = [];
  x = {Supp, Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%
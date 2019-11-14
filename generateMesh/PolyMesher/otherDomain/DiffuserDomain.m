%------------------------------ PolyTop ----------------------------------%
% Ref: A Pereira, C Talischi, GH Paulino, IFM Menezes, MS Carvalho        %
%      "Fluid Flow Topology Optimization in PolyTop: Stability and        %
%      Computational Implementation", Struct Multidisc Optim,             %
%      DOI 10.1007/s00158-014-1182-z                                      %
%-------------------------------------------------------------------------%
function [x] = DiffuserDomain(Demand,Arg)
  r = sqrt((7/18)^2+(1/2)^2);
  BdBox = [0 1 1/2-r 1/2+r];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dCircle(P,7/18,1/2,sqrt((7/18)^2+(1/2)^2));
  d2 = dLine(P,0,0,0,1);
  d3 = dLine(P,1,1,1,0);
  Dist = dDiff(d1,dUnion(d2,d3));
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  [I] = PolyBoundary(Node,Element);
  NodeBC = zeros(2*size(I,1),3);  ElementBC = [];
  for i=1:size(I,1)
    NodeBC(2*i-1,:) = [I(i), 1, 0]; NodeBC(2*i,:) = [I(i), 2, 0];
    x = Node(I(i),1); y = Node(I(i),2);
    r = y - 1/2;
    if     ((abs(x-BdBox(1))<eps) && (y>0) && (y<1))
      NodeBC(2*i-1,:) = [I(i), 1, 1*(1-(r/(1/2))^2)];
    elseif ((abs(x-BdBox(2))<eps) && (y>1/3) && (y<2/3))
      NodeBC(2*i-1,:) = [I(i), 1, 3*(1-(r/(1/6))^2)];
    end
  end
  x = {NodeBC,ElementBC};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%
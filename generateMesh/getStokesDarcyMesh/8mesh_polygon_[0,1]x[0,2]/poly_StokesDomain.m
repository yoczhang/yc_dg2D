%------------------------------ PolyMesher -------------------------------%
%
%-------------------------------------------------------------------------%
function [x] = poly_StokesDomain(Demand,Arg)
  BdBox = [0 1 1 2];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(Arg,BdBox);
  end
end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
end
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
    x = cell(2,1); % No boundary conditions specified for this problem
end
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(h,BdBox)
%   xx = (BdBox(1):h:BdBox(1))';
%   yy = BdBox(3) * ones(length(xx),1);
%   PFix = [xx, yy];
  PFix = [];
end
%-------------------------------------------------------------------------%

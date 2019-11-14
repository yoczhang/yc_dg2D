% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DLevelSetMain(xxElem, ...
    yyElem, ffElem, NodesAct1, NodesAct2, CrackTip, nQ);

% Organize the computation of integration points. Note that almost-polar
% integration is used in the crack-tip element. 

if length(find(NodesAct1==0))==4 & length(find(NodesAct2==0))==4 
    % Standard finite element.
    %patch(xxElem, yyElem, 'b')
    [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DRefElemQuad(2);
elseif isempty(find(xxElem>CrackTip.xx))==0 & isempty(find(xxElem<CrackTip.xx))==0 & ...
    isempty(find(yyElem>CrackTip.yy))==0 & isempty(find(yyElem<CrackTip.yy))==0
    % Crack-tip element.
    %patch(xxElem, yyElem, 'r')
    [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DCrackTipElem(xxElem, yyElem, ffElem, 0, 0, 5);
    %[xxIntRef, yyIntRef, wwIntRef] = IntPoints2DLevelSet(ffElem, nQ);
else % All other elements.
    [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DLevelSet(ffElem, nQ);
end

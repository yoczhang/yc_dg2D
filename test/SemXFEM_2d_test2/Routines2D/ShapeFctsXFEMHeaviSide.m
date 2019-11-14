% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [N, dNdx, dNdy, M, dMdx, dMdy, xxInt, yyInt, wwInt, ffInt] = ...
    ShapeFctsXFEMHeaviSide(xxElem, yyElem, ffElem, NodesAct, xxIntRef, ...
    yyIntRef, wwIntRef, nQ);

% Compute shape functions and enrichment functions at integration 
% points in the real element.

% Get standard FE shape functions and level-set values at int.points.
[N, dNdx, dNdy, xxInt, yyInt, wwInt] = ShapeFctsStrdFEM(xxElem, yyElem, ...
    xxIntRef, yyIntRef, wwIntRef, nQ);
ffInt = N' * ffElem;

% Define the enrichment functions M, dMdx, dMdy for enr. element nodes.
M     = zeros(4, nQ);
dMdx  = zeros(4, nQ);
dMdy  = zeros(4, nQ);
if sum(NodesAct) == 0
    return % Go back if enr. is not active in this element.
end

for i = 1 : nQ 
    % Strd. FEM shape fcts. at int. point.
    NInt  = N(:, i);
    NxInt = dNdx(:, i);
    NyInt = dNdy(:, i);
    % HeaviSide(Levelset) at int. point. and nodes.
    PsiElem = [0; 0; 0; 0];
    PsiInt = 0.5*(1+sign(ffInt(i))); % This the Heaviside-function.
    dPsidxInt = 0;
    dPsidyInt = 0;
    % Enrichment function at int. point.
    M(:, i)    = NInt.*(PsiInt-PsiElem);
    dMdx(:, i) = NxInt.*(PsiInt-PsiElem) + NInt.*dPsidxInt;
    dMdy(:, i) = NyInt.*(PsiInt-PsiElem) + NInt.*dPsidyInt;
end

Pos = find(NodesAct == 0);
M(Pos, :)    = 0;
dMdx(Pos, :) = 0;
dMdy(Pos, :) = 0;

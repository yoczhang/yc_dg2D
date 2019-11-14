% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [N, dNdx, dNdy, M, dMdx, dMdy, xxInt, yyInt, wwInt, ffInt] = ...
    ShapeFctsXFEMAbsMoes(xxElem, yyElem, ffElem, NodesAct, xxIntRef, ...
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
if sum(NodesAct) ~= 4
    return % Go back if this element is not fully enriched.
end

for i = 1 : nQ 
    % Strd. FEM shape fcts. at int. point.
    NInt  = N(:, i);
    NxInt = dNdx(:, i);
    NyInt = dNdy(:, i);
    % Abs(Levelset) at int. point. and nodes.
    PsiElem = abs(ffElem);
    PsiInt = abs(ffInt(i));
    dPsidxInt = sign(ffInt(i)) * NxInt' * ffElem;
    dPsidyInt = sign(ffInt(i)) * NyInt' * ffElem;
    % Enrichment function at int. point.
    PsiMod = PsiElem'*NInt - PsiInt';
    dPsidxMod = PsiElem'*NxInt - dPsidxInt';
    dPsidyMod = PsiElem'*NyInt - dPsidyInt';
    M(:, i)    = NInt.*PsiMod;
    dMdx(:, i) = NxInt.*PsiMod + NInt.*dPsidxMod;
    dMdy(:, i) = NyInt.*PsiMod + NInt.*dPsidyMod;
    % plot3(xxInt(i), yyInt(i), PsiMod, 'k*')
end

Pos = find(NodesAct == 0);
M(Pos, :)    = 0;
dMdx(Pos, :) = 0;
dMdy(Pos, :) = 0;
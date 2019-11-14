% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [N, dNdx, dNdy, M, dMdx, dMdy, F1, dF1dx, dF1dy, F2, dF2dx, dF2dy, ...
    F3, dF3dx, dF3dy, F4, dF4dx, dF4dy, xxInt, yyInt, wwInt, ffInt] = ...
    ShapeFctsXFEMCrack(xxElem, yyElem, ffElem, NodesAct1, NodesAct2, ...
    xxIntRef, yyIntRef, wwIntRef, CrackTip, nQ);

% Compute shape functions, sign- and branch-enrichment functions at 
% integration points in the real element.
%
% Description of some output data:
% N, dNdx, dNdy    : Standard bi-linear FE shape functions.
% M, dMdx, dMdy    : Sign-enriched shape functions.
% Fi, dFidx, dFidy : 4 branch-enriched shape functions.

% Get standard FE shape functions and level-set values at int.points.
[N, dNdx, dNdy, xxInt, yyInt, wwInt] = ShapeFctsStrdFEM(xxElem, yyElem, ...
    xxIntRef, yyIntRef, wwIntRef, nQ);
ffInt = N' * ffElem;

% Initialization.
M  = zeros(4, nQ); dMdx  = zeros(4, nQ); dMdy  = zeros(4, nQ);
F1 = zeros(4, nQ); dF1dx = zeros(4, nQ); dF1dy = zeros(4, nQ);
F2 = zeros(4, nQ); dF2dx = zeros(4, nQ); dF2dy = zeros(4, nQ);
F3 = zeros(4, nQ); dF3dx = zeros(4, nQ); dF3dy = zeros(4, nQ);
F4 = zeros(4, nQ); dF4dx = zeros(4, nQ); dF4dy = zeros(4, nQ);
if sum(NodesAct1) == 0 & sum(NodesAct2) == 0
    return % Go back if enr. are not active in this element.
end

if sum(NodesAct1) ~= 0 % Sign-enrichment is activated for at least one node.
    PsiElem = sign(ffElem);
end
if sum(NodesAct2) ~= 0 % Branch-enrichment is activated for at least one node.
    for i = 1 : 4
        [f1Elem(i,1), df1dxElem(i,1), df1dyElem(i,1), ...
         f2Elem(i,1), df2dxElem(i,1), df2dyElem(i,1), ...
         f3Elem(i,1), df3dxElem(i,1), df3dyElem(i,1), ...
         f4Elem(i,1), df4dxElem(i,1), df4dyElem(i,1)] = ...
            EvalCrackTipEnr(CrackTip, xxElem(i), yyElem(i));
    end
end

% Define the enrichment functions at the integration points in the REAL element.
for i = 1 : nQ 
    % Strd. FEM shape functions at integration point.
    NInt  = N(:, i);
    NxInt = dNdx(:, i);
    NyInt = dNdy(:, i);

    % Sign-enrichment.
    if sum(NodesAct1) ~= 0 % Sign-enrichment is activated for at least one node.
        PsiInt = sign(ffInt(i));
        dPsidxInt = 0;
        dPsidyInt = 0;
        M(:, i)    = NInt.*(PsiInt-PsiElem);
        dMdx(:, i) = NxInt.*(PsiInt-PsiElem) + NInt.*dPsidxInt;
        dMdy(:, i) = NyInt.*(PsiInt-PsiElem) + NInt.*dPsidyInt;
    end

    % Branch-enrichment.
    if sum(NodesAct2) ~= 0 % Branch-enrichment is activated for at least one node.
        [f1Int, df1dxInt, df1dyInt, f2Int, df2dxInt, df2dyInt, ...
            f3Int, df3dxInt, df3dyInt, f4Int, df4dxInt, df4dyInt] = ...
            EvalCrackTipEnr(CrackTip, xxInt(i), yyInt(i));
        F1(:, i)    = NInt.* (f1Int-f1Elem);
        dF1dx(:, i) = NxInt.*(f1Int-f1Elem) + NInt.*df1dxInt;
        dF1dy(:, i) = NyInt.*(f1Int-f1Elem) + NInt.*df1dyInt;
        F2(:, i)    = NInt.* (f2Int-f2Elem);
        dF2dx(:, i) = NxInt.*(f2Int-f2Elem) + NInt.*df2dxInt;
        dF2dy(:, i) = NyInt.*(f2Int-f2Elem) + NInt.*df2dyInt;
        F3(:, i)    = NInt.* (f3Int-f3Elem);
        dF3dx(:, i) = NxInt.*(f3Int-f3Elem) + NInt.*df3dxInt;
        dF3dy(:, i) = NyInt.*(f3Int-f3Elem) + NInt.*df3dyInt;
        F4(:, i)    = NInt.* (f4Int-f4Elem);
        dF4dx(:, i) = NxInt.*(f4Int-f4Elem) + NInt.*df4dxInt;
        dF4dy(:, i) = NyInt.*(f4Int-f4Elem) + NInt.*df4dyInt;
    end
end

Pos = find(NodesAct1 == 0);
M(Pos, :) = 0; dMdx(Pos, :) = 0; dMdy(Pos, :) = 0;

Pos = find(NodesAct2 == 0);
F1(Pos, :) = 0; dF1dx(Pos, :) = 0; dF1dy(Pos, :) = 0;
F2(Pos, :) = 0; dF2dx(Pos, :) = 0; dF2dy(Pos, :) = 0;
F3(Pos, :) = 0; dF3dx(Pos, :) = 0; dF3dy(Pos, :) = 0;
F4(Pos, :) = 0; dF4dx(Pos, :) = 0; dF4dy(Pos, :) = 0;

% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function XFEM2dRodBiMat_AbsEnr(TestCase)

% 2D-XFEM with abs-enrichment for bi-material rod test case.

addpath Routines2D

nQ = 5;  % Number of Gauss integration points (from 2 to 21 possible).

% Get mesh.
nElem = 19;
[NodeNum, ElemNum, xx, yy, Mesh] = GetMesh(1, 1, nElem, nElem);
PlotMesh(Mesh, xx, yy, ElemNum)

% Define test case parameters.
EE1 = 100; EE2 = 10;
nu1 = 0.0; nu2 = 0.0;
mu1 = EE1/(2*(1+nu1));
mu2 = EE2/(2*(1+nu2));
lambda1 = EE1*nu1/((1+nu1)*(1-2*nu1)); % plane strain
lambda2 = EE2*nu2/((1+nu2)*(1-2*nu2)); % plane strain
%lambda1 = EE1*nu1/(1-nu1*nu1); % plane stress
%lambda2 = EE2*nu2/(1-nu2*nu2); % plane stress

if TestCase == 1 % Area force in domain.
    LineLoadRight = 0;
    fx = 1; fy = 0;
elseif TestCase == 2 % Line force on right boundary.
    LineLoadRight = 100;
    fx = 0; fy = 0;
else
    error('This value for TestCase is not defined!')
end

% Get level-set function.
aa=1; bb=-1; cc=-0.5; % Vertical disc. at x=0.5.
ff = GetLevelSet(aa, bb, cc, xx, yy);
PlotLevelSet(Mesh, xx, yy, ff, ElemNum)

% Define Dirichlet boundary conditions.
if TestCase == 1 % Area force in domain.
    uDirNodes  = find(xx==min(xx) | xx==max(xx));
elseif TestCase == 2 % Line force on right boundary.
    uDirNodes  = find(xx==min(xx));
end
uDirNum = length(uDirNodes);
uDirValues = zeros(uDirNum, 1);

vDirNodes  = find(yy==min(yy));
vDirNum = length(vDirNodes);
vDirValues = zeros(vDirNum, 1);

% Get enriched elements and nodes.
[ElemsEnriched, NodesEnriched] = GetEnrichedNodesElems(Mesh, ff);
PlotEnrichment(Mesh, xx, yy, ElemsEnriched, NodesEnriched)
pause

MAT = sparse(4*NodeNum, 4*NodeNum);
RHS = zeros(4*NodeNum, 1);

% Domain integration.
for CurrElem = 1 : ElemNum

    Nodes = Mesh(CurrElem, :);
    xxElem = xx(Nodes)';
    yyElem = yy(Nodes)';
    ffElem = ff(Nodes);

    % Activated nodes are enriched.
    NodesAct = [length(find(NodesEnriched==Nodes(1))) length(find(NodesEnriched==Nodes(2))) ...
        length(find(NodesEnriched==Nodes(3))) length(find(NodesEnriched==Nodes(4)))];

    % Set integration points in REFERENCE element.
    [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DLevelSet(ffElem, nQ);
    Curr_nQ = length(xxIntRef);

    % Get shape functions in element.
    [N, dNdx, dNdy, M, dMdx, dMdy, xxInt, yyInt, wwInt, ffInt] = ...
        ShapeFctsXFEMAbs(xxElem, yyElem, ffElem, NodesAct, xxIntRef, yyIntRef, wwIntRef, Curr_nQ);
    % plot(xxInt, yyInt, 'k*')

    % Integrate PDE.
    [MAT, RHS] = BuildMatRhs_Hooke(N, dNdx, dNdy, M, dMdx, dMdy, ...
        xxInt, yyInt, wwInt, ffInt, Nodes, lambda1, lambda2, ...
        mu1, mu2, fx, fy, Curr_nQ, NodeNum, MAT, RHS);

    % Store the shape function data for computing the L2-norm later.
    ShapeFctData(CurrElem).xxInt = xxInt;
    ShapeFctData(CurrElem).yyInt = yyInt;
    ShapeFctData(CurrElem).wwInt = wwInt;
    ShapeFctData(CurrElem).ffInt = ffInt;
    ShapeFctData(CurrElem).nQ = Curr_nQ;
    ShapeFctData(CurrElem).N = N;
    ShapeFctData(CurrElem).M = M; %enrichment functions.

    if (CurrElem/100) == round(CurrElem/100)
        disp([num2str(CurrElem) ' / ' num2str(ElemNum)])
    end
end

% Integration over boundary of domain (Neumann boundary conditions).
[ElemNum1D, BoundaryElems1D, CorrBoundaryElems2D] = DetectBoundaryElems(...
    NodeNum, ElemNum, xx, yy, Mesh);
for m = 1 : ElemNum1D
    Nodes1D = BoundaryElems1D(m, :); % Get nodes of the 1D boundary element.
    xxElem = xx(Nodes1D);
    yyElem = yy(Nodes1D);

    % Only consider elements on the right.
    if length(find( abs(xxElem-1) < 1.e-12 )) ~= 2
        continue
    end

    hh = max(yyElem) - min(yyElem);
    RHS(Nodes1D) = RHS(Nodes1D) + 0.5 * hh * LineLoadRight;
    % plot(xx(Nodes1D), yy(Nodes1D), 'bo')
end

% Insert Dirichlet BCs.
MAT(uDirNodes, :)         = 0;
MAT(uDirNodes, uDirNodes) = speye(uDirNum);
MAT(vDirNodes+NodeNum, :) = 0;
MAT(vDirNodes+NodeNum, vDirNodes+NodeNum) = speye(vDirNum);
RHS(uDirNodes)            = uDirValues;
RHS(vDirNodes+NodeNum)    = vDirValues;

% Reduce system of equations.
Pos = [[1:1:2*NodeNum]'; NodesEnriched+2*NodeNum; NodesEnriched+3*NodeNum];
MAT = MAT(Pos, Pos);
RHS = RHS(Pos);

%disp(sprintf('Condition number of final system         : %15.5e', condest(MAT)))

% Solve system of equations for solution.
Sol = zeros(4*NodeNum, 1);
Sol(Pos) = MAT \ RHS;

% Plot solution.
if TestCase == 1
    Scal = 50;
elseif TestCase == 2
    Scal = 1;
end
uu = Sol(1:NodeNum); vv = Sol(NodeNum+1:2*NodeNum);
PlotMesh(Mesh, xx+Scal*uu, yy+Scal*vv, ElemNum)
pause(0.5)

% Compute L2-norm.
[L2Norm] = GetL2Norm_RodBiMat(ElemNum, NodeNum, Mesh, ShapeFctData, Sol, EE1, EE2, ...
    LineLoadRight, TestCase, cc);
disp(sprintf('Resulting L2-error                       : %18.8e', L2Norm))
view(10,10)

rmpath Routines2D


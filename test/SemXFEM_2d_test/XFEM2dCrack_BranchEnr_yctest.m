% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function XFEM2dCrack_BranchEnr_yctest

% 2D-XFEM with sign-enrichment along crack and branch-enrichment 
% at crack-tip.

% addpath Routines2D
% addpath RoutinesCrack
clc
clearvars
close all


rr = 0.2; % Radius of branch-enrichment around crack-tip.
nQ = 5;  % Number of Gauss integration points (from 2 to 21 possible).

% Get mesh.
nElem = 3;
[NodeNum, ElemNum, xx, yy, Mesh] = GetMesh(1, 1, nElem, nElem);
xx = xx*2-1; yy = yy*2-1;
% %------------------------------------------
% yc_node = [xx,yy];
% yc_elem = Mesh;
% yc_meshInfo = polyMeshAuxStructure(yc_node, yc_elem);
% plotPolyMsh(yc_meshInfo)
% %-------------------------------------------
% figure(2)
PlotMesh(Mesh, xx, yy, ElemNum)
% divideMesh2Submesh_test(yc_node, yc_elem);


% Define test case parameters.
k1 = 1.0;
EE = 10000;
nu = 0.3;
fx = 0; fy = 0; 
mu = EE/(2*(1+nu));
%lambda = EE*nu/((1+nu)*(1-2*nu)); % plane strain
%kappa = 3-4*nu; % plane strain
lambda = EE*nu/(1-nu*nu); % plane stress
kappa  = (3-nu)/(1+nu); % plane stress

% Get level-set function.
ff = yy;
% ffy = yy;
% ffx = xx-0.52;
% PlotLevelSet_yctest(Mesh, xx, yy, ffx, ffy, ElemNum)
PlotLevelSet(Mesh, xx, yy, ff, ElemNum);

% Get exact solution at the nodes.
uuExact = zeros(NodeNum, 1); vvExact = zeros(NodeNum, 1);
for i = 1 : NodeNum
    [uuExact(i), vvExact(i), duudx, duudy, dvvdx, dvvdy, Eps11, Eps12, Eps22, Sigma11, Sigma12, Sigma22] = ...
        ExactSol_Mode1(xx(i), yy(i), k1, kappa, mu, lambda);
end

% Define Dirichlet boundary conditions.
Bound = find(xx==min(xx) | xx==max(xx) | yy==min(yy) | yy==max(yy));
uDirNodes = Bound;
vDirNodes = Bound;

uDirNum = length(uDirNodes);
vDirNum = length(vDirNodes);

uDirValues = uuExact(uDirNodes);
vDirValues = vvExact(vDirNodes);
% Scal = 1000;
% plot(xx(uDirNodes)+uDirValues*Scal, yy(uDirNodes)+vDirValues*Scal, 'k*')
% axis([-1.2 1.2 -1.2 1.2])

% Define position and orientation of crack tip (must fit to the 
% level-set function).
CrackTip.xx = 0;
CrackTip.yy = 0;
CrackTip.th = 0;

% Get enriched nodes and elements around crack-tip.
Dist = sqrt(xx.^2+yy.^2);
NodesEnriched2 = find(Dist-rr<0);
Count = 0;
for i = 1 : ElemNum
    Nodes = Mesh(i, :);
    if isempty(setdiff(Nodes, NodesEnriched2)) == 1
        Count = Count + 1;
        ElemsEnriched2(Count, 1) = i;
    end
end

% Get enriched nodes and elements along crack.
[ElemsEnriched1Dummy, NodesEnriched1] = GetEnrichedNodesElems(Mesh, ff);
Pos = find(xx(NodesEnriched1)<0);
NodesEnriched1 = NodesEnriched1(Pos);
NodesEnriched1 = setdiff(NodesEnriched1, NodesEnriched2);
Count = 0;
for i = 1 : ElemNum
    Nodes = Mesh(i, :);
    if isempty(setdiff(Nodes, NodesEnriched1)) == 1
        Count = Count + 1;
        ElemsEnriched1(Count, 1) = i;
    end
end

PlotEnrichmentCrack(Mesh, xx, yy, ElemsEnriched1, NodesEnriched1, ElemsEnriched2, NodesEnriched2)
pause

MAT = sparse(12*NodeNum, 12*NodeNum);
RHS = zeros(12*NodeNum, 1);

% Domain integration.
for CurrElem = 1 : ElemNum

    Nodes = Mesh(CurrElem, :);
    xxElem = xx(Nodes)';
    yyElem = yy(Nodes)';
    ffElem = ff(Nodes);

    % Activated nodes for sign-enrichment.
    NodesAct1 = [length(find(NodesEnriched1==Nodes(1))) length(find(NodesEnriched1==Nodes(2))) ...
        length(find(NodesEnriched1==Nodes(3))) length(find(NodesEnriched1==Nodes(4)))];

    % Activated nodes for branch-enrichment.
    NodesAct2 = [length(find(NodesEnriched2==Nodes(1))) length(find(NodesEnriched2==Nodes(2))) ...
        length(find(NodesEnriched2==Nodes(3))) length(find(NodesEnriched2==Nodes(4)))];

    % Set integration points in REFERENCE element.
    [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DLevelSetMain(xxElem, yyElem, ffElem, ...
        NodesAct1, NodesAct2, CrackTip, nQ);
    Curr_nQ = length(xxIntRef);

    % Get shape functions in element.
    [N, dNdx, dNdy, M, dMdx, dMdy, F1, dF1dx, dF1dy, F2, dF2dx, dF2dy, ...
        F3, dF3dx, dF3dy, F4, dF4dx, dF4dy, xxInt, yyInt, wwInt, ffInt] = ...
        ShapeFctsXFEMCrack(xxElem, yyElem, ffElem, NodesAct1, NodesAct2, ...
        xxIntRef, yyIntRef, wwIntRef, CrackTip, Curr_nQ);
    %plot(xxInt, yyInt, 'k*')

    % Integrate PDE.
    [MAT, RHS] = BuildMatRhs_HookeCrack(N, dNdx, dNdy, M, dMdx, dMdy, ...
        F1, dF1dx, dF1dy, F2, dF2dx, dF2dy, F3, dF3dx, dF3dy, ...
        F4, dF4dx, dF4dy, xxInt, yyInt, wwInt, Nodes, lambda, ...
        mu, fx, fy, Curr_nQ, NodeNum, MAT, RHS);

    % Store the shape function data for computing the L2-norm later.
    ShapeFctData(CurrElem).xxInt = xxInt;
    ShapeFctData(CurrElem).yyInt = yyInt;
    ShapeFctData(CurrElem).wwInt = wwInt;
    ShapeFctData(CurrElem).ffInt = ffInt;
    ShapeFctData(CurrElem).nQ = Curr_nQ;
    ShapeFctData(CurrElem).N  =  N; ShapeFctData(CurrElem).Nx  =  dNdx; ShapeFctData(CurrElem).Ny  =  dNdy;
    ShapeFctData(CurrElem).M  =  M; ShapeFctData(CurrElem).Mx  =  dMdx; ShapeFctData(CurrElem).My  =  dMdy;
    ShapeFctData(CurrElem).F1 = F1; ShapeFctData(CurrElem).F1x = dF1dx; ShapeFctData(CurrElem).F1y = dF1dy;
    ShapeFctData(CurrElem).F2 = F2; ShapeFctData(CurrElem).F2x = dF2dx; ShapeFctData(CurrElem).F2y = dF2dy;
    ShapeFctData(CurrElem).F3 = F3; ShapeFctData(CurrElem).F3x = dF3dx; ShapeFctData(CurrElem).F3y = dF3dy;
    ShapeFctData(CurrElem).F4 = F4; ShapeFctData(CurrElem).F4x = dF4dx; ShapeFctData(CurrElem).F4y = dF4dy;

    if (CurrElem/100) == round(CurrElem/100)
        disp([num2str(CurrElem) ' / ' num2str(ElemNum)])
    end
end

% Insert Dirichlet BCs.
MAT(uDirNodes, :)         = 0;
MAT(uDirNodes, uDirNodes) = speye(uDirNum);
MAT(vDirNodes+NodeNum, :) = 0;
MAT(vDirNodes+NodeNum, vDirNodes+NodeNum) = speye(vDirNum);
RHS(uDirNodes)            = uDirValues;
RHS(vDirNodes+NodeNum)    = vDirValues;

% Reduce system of equations.
Pos = [[1:1:2*NodeNum]'; ...
    NodesEnriched1+ 2*NodeNum; NodesEnriched1+ 3*NodeNum; ...
    NodesEnriched2+ 4*NodeNum; NodesEnriched2+ 5*NodeNum; NodesEnriched2+ 6*NodeNum; NodesEnriched2+ 7*NodeNum; ...
    NodesEnriched2+ 8*NodeNum; NodesEnriched2+ 9*NodeNum; NodesEnriched2+10*NodeNum; NodesEnriched2+11*NodeNum];
MAT = MAT(Pos, Pos);
RHS = RHS(Pos);

%disp(sprintf('Condition number of final system         : %15.5e', condest(MAT)))

% Solve system of equations for solution.
Sol = zeros(12*NodeNum, 1);
Sol(Pos) = MAT \ RHS;

uuTotal = Sol(1:NodeNum);
vvTotal = Sol(NodeNum+1:2*NodeNum);
disp(sprintf('Nodal norm of the solution               : %18.8e', norm(uuTotal-uuExact)+norm(vvTotal-vvExact)))

% Plot solution.
Scal = 1000;
PlotMesh(Mesh, xx+Scal*uuTotal, yy+Scal*vvTotal, ElemNum)

% Compute L2-norm.
[L2Norm, EnergyNorm] = GetL2Norm_Mode1Branch(ElemNum, NodeNum, Mesh, ShapeFctData, ...
    Sol, k1, kappa, mu, lambda);
disp(sprintf('Resulting L2-norm                        : %18.8e', L2Norm))
disp(sprintf('Resulting energy-norm                    : %18.8e', EnergyNorm))
disp(sprintf('   %1.15e %1.15e', L2Norm, EnergyNorm))

% Compute stress intensity factors.
[k1_1, k2_1] = StressIntFactorRadial(Mesh, xx, yy, Sol, ShapeFctData, ...
    0.2, k1, EE, nu, kappa, mu, lambda, NodeNum, ElemNum);
disp(sprintf('Stress intensity factor, rr=0.2 (exact)  : %17.7e  (%15.5e)', k1_1, k1))

[k1_2, k2_2] = StressIntFactorRadial(Mesh, xx, yy, Sol, ShapeFctData, ...
    0.4, k1, EE, nu, kappa, mu, lambda, NodeNum, ElemNum);
disp(sprintf('Stress intensity factor, rr=0.4 (exact)  : %17.7e  (%15.5e)', k1_2, k1))

[k1_3, k2_3] = StressIntFactorRadial(Mesh, xx, yy, Sol, ShapeFctData, ...
    0.6, k1, EE, nu, kappa, mu, lambda, NodeNum, ElemNum);
disp(sprintf('Stress intensity factor, rr=0.6 (exact)  : %17.7e  (%15.5e)', k1_3, k1))
disp(sprintf('%17.7e %14.7e %14.7e', k1_1, k1_2, k1_3))

rmpath Routines2D
rmpath RoutinesCrack

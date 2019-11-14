% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function XFEM2dVoid_HeaviEnr

% 2D-XFEM with sign-enrichment for plate test case with a circular
% hole. The continuity is enforced by 1d Lagrange multipliers.

addpath Routines2D

nQ = 5;  % Number of Gauss integration points (from 2 to 21 possible).

% Get mesh.
nElem = 20;
[NodeNum, ElemNum, xx, yy, Mesh] = GetMesh(1, 1, nElem, nElem);
xx = xx*2-1; yy = yy*2-1;
PlotMesh(Mesh, xx, yy, ElemNum)

% Define test case parameters.
aa = 0.401;
EE1 = 10000; EE2 = 0; 
nu1 = 0.3; nu2 = 0; 
fx = 0; fy = 0; 
mu1 = EE1/(2*(1+nu1));
mu2 = EE2/(2*(1+nu2));
lambda1 = EE1*nu1/((1+nu1)*(1-2*nu1)); % plane strain
lambda2 = EE2*nu2/((1+nu2)*(1-2*nu2)); % plane strain
kappa1  = 3-4*nu1; % plane strain
kappa2  = 3-4*nu2; % plane strain
% lambda1 = EE1*nu1/(1-nu1*nu1); % plane stress
% lambda2 = EE2*nu2/(1-nu2*nu2); % plane stress
% kappa1  = (3-nu1)/(1+nu1); % plane stress
% kappa2  = (3-nu2)/(1+nu2); % plane stress

% Get level-set function.
Dist = sqrt( xx.^2 + yy.^2 );
ff = Dist - aa;
PlotLevelSet(Mesh, xx, yy, ff, ElemNum)

% Get exact solution at the nodes.
uuExact = zeros(NodeNum, 1); vvExact = zeros(NodeNum, 1);
for i = 1 : NodeNum
    if ff(i) > 0
        [uuExact(i), vvExact(i), Eps11, Eps12, Eps22, Sigma11, Sigma12, Sigma22] = ...
            ExactSol_Hole(xx(i), yy(i), aa, kappa1, lambda1, mu1);
    end
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

% Get enriched elements and nodes.
[ElemsEnriched, NodesEnriched] = GetEnrichedNodesElems(Mesh, ff);
PlotEnrichment(Mesh, xx, yy, ElemsEnriched, NodesEnriched)
pause(0.1)

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
        ShapeFctsXFEMHeaviSide(xxElem, yyElem, ffElem, NodesAct, xxIntRef, yyIntRef, wwIntRef, Curr_nQ);
    %plot(xxInt, yyInt, 'k*')

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
    ShapeFctData(CurrElem).N = N; ShapeFctData(CurrElem).Nx = dNdx; ShapeFctData(CurrElem).Ny = dNdy;
    ShapeFctData(CurrElem).M = M; ShapeFctData(CurrElem).Mx = dMdx; ShapeFctData(CurrElem).My = dMdy;

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

RemoveNodes = find(ff < 0);
for i = 1 : length(ElemsEnriched)
    Nodes = Mesh(ElemsEnriched(i), :);
    RemoveNodes = setdiff(RemoveNodes, Nodes);
end
SolveNodes = setdiff([1:1:NodeNum], RemoveNodes)';
InnerSolveNodes = setdiff(SolveNodes, find(ff>0));

% plot(xx(RemoveNodes), yy(RemoveNodes), 'b*')
% plot(xx(SolveNodes) , yy(SolveNodes) , 'bo')
% plot(xx(InnerSolveNodes) , yy(InnerSolveNodes) , 'bs')

% Set Dirichlet BCs for enriched elements: prescribe u1=u2=u3=u4=0 but leave
% enrichments a1, a2, a3, a4 free; same for v_i and b_i.
MAT(NodesEnriched+0*NodeNum, :) = 0;
MAT(NodesEnriched+0*NodeNum, NodesEnriched+0*NodeNum) = speye(length(NodesEnriched));
MAT(NodesEnriched+1*NodeNum, :) = 0;
MAT(NodesEnriched+1*NodeNum, NodesEnriched+1*NodeNum) = speye(length(NodesEnriched));

% Reduce system of equations.
Pos = [SolveNodes; SolveNodes+NodeNum; NodesEnriched+2*NodeNum; NodesEnriched+3*NodeNum];
MAT = MAT(Pos, Pos);
RHS = RHS(Pos);

%disp(sprintf('Condition number of final system         : %15.5e', condest(MAT)))

% Solve system of equations for solution.
Sol = zeros(4*NodeNum, 1);
Sol(Pos) = MAT \ RHS;
uuTotal = Sol(1:NodeNum);
vvTotal = Sol(NodeNum+1:2*NodeNum);

% Plot solution.
Scal = 1000;
PlotMesh(Mesh, xx+Scal*uuTotal, yy+Scal*vvTotal, ElemNum)

% Compute L2-norm.
pause(0.5)
[L2Norm, EnergyNorm] = GetL2Norm_Hole(ElemNum, NodeNum, Mesh, ShapeFctData, Sol, ...
    aa, mu1, lambda1, EE1, nu1, kappa1);
disp(sprintf('Resulting L2-error                       : %18.8e', L2Norm))
disp(sprintf('Resulting energy-norm                    : %18.8e', EnergyNorm))
disp(sprintf('   %1.5e %1.5e', L2Norm, EnergyNorm))

rmpath Routines2D

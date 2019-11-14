% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function FEM2dVoid

% 2D standard FEM for plate test case with a circular hole. Integration
% points inside the hole are neglected.

addpath Routines2D

nQ = 5;  % Number of Gauss integration points (from 2 to 21 possible).

% Get mesh.
nElem = 20;
[NodeNum, ElemNum, xx, yy, Mesh] = GetMesh(1, 1, nElem, nElem);
xx = xx*2-1; yy = yy*2-1;
PlotMesh(Mesh, xx, yy, ElemNum)

% Define test case parameters.
aa = 0.401;
EE = 10000;
nu = 0.3;
fx = 0; fy = 0; 
mu = EE/(2*(1+nu));
lambda = EE*nu/((1+nu)*(1-2*nu)); % plane strain
kappa  = 3-4*nu; % plane strain
% lambda = EE*nu/(1-nu*nu); % plane stress
% kappa  = (3-nu)/(1+nu); % plane stress

% Get level-set function.
Dist = sqrt( xx.^2 + yy.^2 );
ff = Dist - aa;
PlotLevelSet(Mesh, xx, yy, ff, ElemNum)
pause

% Get exact solution at the nodes.
uuExact = zeros(NodeNum, 1); vvExact = zeros(NodeNum, 1);
for i = 1 : NodeNum
    if ff(i) > 0
        [uuExact(i), vvExact(i), Eps11, Eps12, Eps22, Sigma11, Sigma12, Sigma22] = ...
            ExactSol_Hole(xx(i), yy(i), aa, kappa, lambda, mu);
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

MAT = sparse(2*NodeNum, 2*NodeNum);
RHS = zeros(2*NodeNum, 1);

% Domain integration.
for CurrElem = 1 : ElemNum

    Nodes = Mesh(CurrElem, :);
    xxElem = xx(Nodes)';
    yyElem = yy(Nodes)';
    ffElem = ff(Nodes);

    % Set integration points in REFERENCE element.
    [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DLevelSet(ffElem, nQ);
    Curr_nQ = length(xxIntRef);

    % Get shape functions in element.
    [N, dNdx, dNdy, xxInt, yyInt, wwInt] = ShapeFctsStrdFEM(...
        xxElem, yyElem, xxIntRef, yyIntRef, wwIntRef, Curr_nQ);
    ffInt = N' * ffElem;
    % plot(xxInt, yyInt, 'k*')

    % Set integration weights of nodes in the hole to zero.
    Pos = find(ffInt < 0);
    wwInt(Pos) = 0;

    % Integrate PDE.
    [MAT, RHS] = BuildMatRhs_HookeStrdFEM(N, dNdx, dNdy, xxInt, yyInt, ...
        wwInt, Nodes, lambda, mu, fx, fy, Curr_nQ, NodeNum, MAT, RHS);

    % Store the shape function data for computing the L2-norm later.
    ShapeFctData(CurrElem).xxInt = xxInt;
    ShapeFctData(CurrElem).yyInt = yyInt;
    ShapeFctData(CurrElem).wwInt = wwInt;
    ShapeFctData(CurrElem).ffInt = ffInt;
    ShapeFctData(CurrElem).nQ = Curr_nQ;
    ShapeFctData(CurrElem).N = N; ShapeFctData(CurrElem).Nx = dNdx; ShapeFctData(CurrElem).Ny = dNdy;
    ShapeFctData(CurrElem).M  = zeros(4, Curr_nQ);
    ShapeFctData(CurrElem).Mx = zeros(4, Curr_nQ);
    ShapeFctData(CurrElem).My = zeros(4, Curr_nQ);

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
for i = 1 : ElemNum 
    Nodes = Mesh(i, :);
    if min(sign(ff(Nodes))) ~= max(sign(ff(Nodes)))
        RemoveNodes = setdiff(RemoveNodes, Nodes);
    end
end
SolveNodes = setdiff([1:1:NodeNum], RemoveNodes)';
plot(xx(SolveNodes), yy(SolveNodes), 'bo')
pause

% Reduce system of equations.
Pos = [SolveNodes; SolveNodes+NodeNum];
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
    aa, mu, lambda, EE, nu, kappa);
disp(sprintf('Resulting L2-error                       : %18.8e', L2Norm))
disp(sprintf('Resulting energy-norm                    : %18.8e', EnergyNorm))
disp(sprintf('   %1.5e %1.5e', L2Norm, EnergyNorm))

rmpath Routines2D

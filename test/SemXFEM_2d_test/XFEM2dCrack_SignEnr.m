% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function XFEM2dCrack_SignEnr

% 2D-XFEM with sign-enrichment along crack (NO branch-enrichment 
% at crack-tip).

addpath Routines2D
addpath RoutinesCrack

nQ = 5;  % Number of Gauss integration points (from 2 to 21 possible).

% Get mesh.
% Note: The element number in x-direction must be EVEN, because due to
% the sign-enrichment, the crack is elongated onto the element edge.
% The element number in y-direction must be ODD, so that the crack does
% not align with the element edges.
nElem = 19; % Element number in y-direction.
[NodeNum, ElemNum, xx, yy, Mesh] = GetMesh(1, 1, nElem-1, nElem);
xx = xx*2-1; yy = yy*2-1;
PlotMesh(Mesh, xx, yy, ElemNum)

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
PlotLevelSet(Mesh, xx, yy, ff, ElemNum)

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

% Get enriched elements and nodes.
[ElemsEnriched, NodesEnriched] = GetEnrichedNodesElems(Mesh, ff);

% Modify ElemsEnriched and NodesEnriched, so that the discontinuity
% goes into domain center only.
Count = 0;
MarkerVect = zeros(length(ElemsEnriched), 1);
for i = 1 : length(ElemsEnriched)
    Nodes = Mesh(ElemsEnriched(i), :);
    xxElem = xx(Nodes);
    Pos = find(xxElem<0);
    if length(Pos)==4
        MarkerVect(i) = 1;
    end
end
ElemsEnriched = ElemsEnriched(find(MarkerVect == 1));

Pos = find(xx(NodesEnriched)<0);
NodesEnriched = NodesEnriched(Pos);

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
        ShapeFctsXFEMSign(xxElem, yyElem, ffElem, NodesAct, xxIntRef, yyIntRef, wwIntRef, Curr_nQ);
    %plot(xxInt, yyInt, 'k*')

    % Integrate PDE.
    [MAT, RHS] = BuildMatRhs_Hooke(N, dNdx, dNdy, M, dMdx, dMdy, ...
        xxInt, yyInt, wwInt, ffInt, Nodes, lambda, lambda, ...
        mu, mu, fx, fy, Curr_nQ, NodeNum, MAT, RHS);

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

uuTotal = Sol(1:NodeNum);
vvTotal = Sol(NodeNum+1:2*NodeNum);
disp(sprintf('Nodal norm of the solution               : %18.8e', norm(uuTotal-uuExact)+norm(vvTotal-vvExact)))

% Plot solution.
Scal = 1000;
PlotMesh(Mesh, xx+Scal*uuTotal, yy+Scal*vvTotal, ElemNum)

% Compute L2-norm.
[L2Norm] = GetL2Norm_Mode1Sign(ElemNum, NodeNum, Mesh, ShapeFctData, Sol, k1, ...
    kappa, mu, lambda);
disp(sprintf('Resulting L2-error                       : %18.8e', L2Norm))

rmpath Routines2D
rmpath RoutinesCrack


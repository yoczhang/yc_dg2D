% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function XFEM2dHoleBiMat_SignEnr_1dLM

% 2D-XFEM with sign-enrichment for bi-material test-case with a circular
% inclusion. The continuity is enforced by 1d Lagrange multipliers.

addpath Routines2D

nQ = 5;  % Number of Gauss integration points (from 2 to 21 possible).

% Get mesh.
nElem = 20;
[NodeNum, ElemNum, xx, yy, Mesh] = GetMesh(1, 1, nElem, nElem);
xx = xx*2-1; yy = yy*2-1;
PlotMesh(Mesh, xx, yy, ElemNum)

% Define test case parameters.
aa = 0.401;
bb = 2.0;
EE1  = 1   ; EE2  = 10 ; 
nu1 = 0.25; nu2 = 0.3; 
fx = 0; fy = 0; 
mu1 = EE1/(2*(1+nu1));
mu2 = EE2/(2*(1+nu2));
lambda1 = EE1*nu1/((1+nu1)*(1-2*nu1)); % plane strain
lambda2 = EE2*nu2/((1+nu2)*(1-2*nu2)); % plane strain
%lambda1 = EE1*nu1/(1-nu1*nu1); % plane stress
%lambda2 = EE2*nu2/(1-nu2*nu2); % plane stress
delta = (lambda1+mu1+mu2)*bb*bb / ((lambda2+mu2)*aa*aa+(lambda1+mu1)*(bb*bb-aa*aa)+mu2*bb*bb);

% Get level-set function.
Dist = sqrt( xx.^2 + yy.^2 );
ff = -(Dist - aa);
PlotLevelSet(Mesh, xx, yy, ff, ElemNum)

% Get exact solution at the nodes.
uuExact = zeros(NodeNum, 1); vvExact = zeros(NodeNum, 1);
for i = 1 : NodeNum
    [uuExact(i), vvExact(i), Eps11, Eps12, Eps22, Sigma11, Sigma12, Sigma22] = ...
        ExactSol_BiMat(xx(i), yy(i), aa, bb, delta, mu1, mu2, lambda1, lambda2);
end

% Define Dirichlet boundary conditions.
Bound = find(xx==min(xx) | xx==max(xx) | yy==min(yy) | yy==max(yy));
uDirNodes = Bound;
vDirNodes = Bound;

uDirNum = length(uDirNodes);
vDirNum = length(vDirNodes);

uDirValues = uuExact(uDirNodes);
vDirValues = vvExact(vDirNodes);
% Scal = 0.1;
% plot(xx(uDirNodes)+uDirValues*Scal, yy(uDirNodes)+vDirValues*Scal, 'k*')
% axis([-1.2 1.2 -1.2 1.2])

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
    [ffElem] = ModifyLevelSetElem(ffElem);

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

% Enforce continuity with 1d Lagrange multipliers.
% Place unknowns at intersection points of the discontinuity with element
% edges. Use 1d shape functions for Lagrange multipliers along discontinuity.
xxLagrMult = []; yyLagrMult = []; Count = 0;
[zz, ww] = IntPoints1DGauss(nQ);
for i = 1 : length(ElemsEnriched)  % Loop over cut elements.
    
    CurrElem = ElemsEnriched(i);
    Nodes = Mesh(CurrElem, :);
    xxElem = xx(Nodes)';
    yyElem = yy(Nodes)';
    ffElem = ff(Nodes);
    [ffElem] = ModifyLevelSetElem(ffElem); % THIS IS CRUCIAL HERE.

    % Activated nodes are enriched.
    NodesAct = [length(find(NodesEnriched==Nodes(1))) length(find(NodesEnriched==Nodes(2))) ...
        length(find(NodesEnriched==Nodes(3))) length(find(NodesEnriched==Nodes(4)))];

    % Set integration points in REFERENCE element along discontinuity.
    [xxIntRef, yyIntRef, wwInt, nxVect, nyVect, xx1, yy1, xx2, yy2, nxRef, nyRef] = ...
        IntPoints2DAlongLevelSet(xxElem, yyElem, ffElem, nQ);
    Curr_nQ = length(xxIntRef);
    plot(xx1, yy1, 'g*')
    plot(xx2, yy2, 'g*')

    % Evaluate shape functions on one side of the disc.
    xxIntRef1 = xxIntRef + 1.e-6*nxRef;
    yyIntRef1 = yyIntRef + 1.e-6*nyRef;
    if isempty(find(xxIntRef1>=1 | xxIntRef1<=-1 | yyIntRef1>=1 | yyIntRef1<=-1)) == 0
        error('In Lagrange mult. loop: Reference integration point out of reference element!')
    end
    [N, dNdx, dNdy, M1, dM1dx, dM1dy, xxInt, yyInt, wwIntDummy, ffInt1] = ...
        ShapeFctsXFEMSign(xxElem, yyElem, ffElem, NodesAct, xxIntRef1, yyIntRef1, zeros(Curr_nQ, 1), Curr_nQ);
    plot(xxInt, yyInt, 'k*')

    % Evaluate shape functions on the other side of the disc.
    xxIntRef2 = xxIntRef - 1.e-6*nxRef;
    yyIntRef2 = yyIntRef - 1.e-6*nyRef;
    if isempty(find(xxIntRef2>=1 | xxIntRef2<=-1 | yyIntRef2>=1 | yyIntRef2<=-1)) == 0
        error('In Lagrange mult. loop: Reference integration point out of reference element!')
    end
    [N, dNdx, dNdy, M2, dM2dx, dM2dy, xxInt, yyInt, wwIntDummy, ffInt2] = ...
        ShapeFctsXFEMSign(xxElem, yyElem, ffElem, NodesAct, xxIntRef2, yyIntRef2, zeros(Curr_nQ, 1), Curr_nQ);
    plot(xxInt, yyInt, 'b*')

    % Ensure that the level-set functions have the same sign on each side, but opposite 
    % between both sides.
    if sum(abs(sign(ffInt1)))~=nQ | sum(abs(sign(ffInt2)))~=nQ | ...
        min(sign(ffInt2))~=max(sign(ffInt2))
        [ffInt1 ffInt2]
        error('Level-set data or integration points are wrong!')
    end

    % Construct new Lagrange multiplier if not yet existing.
    if isempty(find(abs(xxLagrMult-xx1)<1.e-12 & abs(yyLagrMult-yy1)<1.e-12)) == 1
        Count = Count + 1;
        xxLagrMult(Count, 1) = xx1;
        yyLagrMult(Count, 1) = yy1;
        MAT(4*NodeNum+2*Count-1, 4*NodeNum+2*Count-1) = 0; % Extend system matrix.
        MAT(4*NodeNum+2*Count  , 4*NodeNum+2*Count  ) = 0; % Extend system matrix.
    end
    if isempty(find(abs(xxLagrMult-xx2)<1.e-12 & abs(yyLagrMult-yy2)<1.e-12)) == 1
        Count = Count + 1;
        xxLagrMult(Count, 1) = xx2;
        yyLagrMult(Count, 1) = yy2;
        MAT(4*NodeNum+2*Count-1, 4*NodeNum+2*Count-1) = 0; % Extend system matrix.
        MAT(4*NodeNum+2*Count  , 4*NodeNum+2*Count  ) = 0; % Extend system matrix.
    end
    
    % Get positions of the current two Lagrange multipliers.
    Pos1 = find(abs(xxLagrMult-xx1)<1.e-12 & abs(yyLagrMult-yy1)<1.e-12);
    Pos2 = find(abs(xxLagrMult-xx2)<1.e-12 & abs(yyLagrMult-yy2)<1.e-12);
    
    if length(Pos1)~=1 | length(Pos2)~=1
        error('Internal error.')
    end
    
    ElemMAT = zeros(2,4);
    MM1 = 0.5*(1-zz);
    MM2 = 0.5*(1+zz);
    Sum = 0;
    for j = 1 : Curr_nQ
        ElemMAT(1,:) = ElemMAT(1,:) + wwInt(j) * MM1(j)*(M1(:, j)'-M2(:, j)');
        ElemMAT(2,:) = ElemMAT(2,:) + wwInt(j) * MM2(j)*(M1(:, j)'-M2(:, j)');
        Sum = Sum + wwInt(j);
    end
    %Sum

    % Enforce continuity on u.
    MAT(4*NodeNum+2*[Pos1 Pos2]-1, Nodes+2*NodeNum) = MAT(4*NodeNum+2*[Pos1 Pos2]-1, Nodes+2*NodeNum) + ElemMAT;
    MAT(Nodes+2*NodeNum, 4*NodeNum+2*[Pos1 Pos2]-1) = MAT(Nodes+2*NodeNum, 4*NodeNum+2*[Pos1 Pos2]-1) + ElemMAT';

    % Enforce continuity on v.
    MAT(4*NodeNum+2*[Pos1 Pos2]  , Nodes+3*NodeNum) = MAT(4*NodeNum+2*[Pos1 Pos2]  , Nodes+3*NodeNum) + ElemMAT;
    MAT(Nodes+3*NodeNum, 4*NodeNum+2*[Pos1 Pos2]  ) = MAT(Nodes+3*NodeNum, 4*NodeNum+2*[Pos1 Pos2]  ) + ElemMAT';
    
    % Draw the shape function for the qq-th Lagrange multiplier.
%     qq = 2;
%     if isempty(find(Pos1==qq))==0
%         plot3(xxInt, yyInt, MM1, 'k*')
%     elseif isempty(find(Pos2==qq))==0
%         plot3(xxInt, yyInt, MM2, 'k*')
%     end
end
LagrMultNum = 2*Count;
RHS(4*NodeNum+[1:LagrMultNum]) = 0; %Set all Lagrange multipliers to zero.

% Insert Dirichlet BCs.
MAT(uDirNodes, :)         = 0;
MAT(uDirNodes, uDirNodes) = speye(uDirNum);
MAT(vDirNodes+NodeNum, :) = 0;
MAT(vDirNodes+NodeNum, vDirNodes+NodeNum) = speye(vDirNum);
RHS(uDirNodes)            = uDirValues;
RHS(vDirNodes+NodeNum)    = vDirValues;

% Reduce system of equations.
Pos = [[1:1:2*NodeNum]'; NodesEnriched+2*NodeNum; NodesEnriched+3*NodeNum; [1:1:LagrMultNum]'+4*NodeNum];
MAT = MAT(Pos, Pos);
RHS = RHS(Pos);

%disp(sprintf('Condition number of final system         : %15.5e', condest(MAT)))

% Solve system of equations for solution.
Sol = zeros(4*NodeNum+LagrMultNum, 1);
Sol(Pos) = MAT \ RHS;

% For nodal norm: Remove node at (0,0) because uuExact = vvExact = NaN.
uuTotal = Sol(1:NodeNum); vvTotal = Sol(NodeNum+1:2*NodeNum); 
PosCenter = find(xx==0 & yy==0);
Pos = setdiff([1:NodeNum], PosCenter);
disp(sprintf('Nodal norm of the solution               : %18.8e', ...
    norm(uuTotal(Pos)-uuExact(Pos))+norm(vvTotal(Pos)-vvExact(Pos))))

% Plot solution.
Scal = 1;
PlotMesh(Mesh, xx+Scal*uuTotal, yy+Scal*vvTotal, ElemNum)

% disp(sprintf('Solution of the %d Lagrange multipliers:', LagrMultNum))
% Sol(4*NodeNum+1:4*NodeNum+LagrMultNum)

% Compute L2-norm.
[L2Norm] = GetL2Norm_HoleBiMat(ElemNum, NodeNum, Mesh, ShapeFctData, Sol, ...
    aa, bb, delta, mu1, mu2, lambda1, lambda2, EE1, EE2, nu1, nu2);
disp(sprintf('Resulting L2-error                       : %18.8e', L2Norm))

rmpath Routines2D

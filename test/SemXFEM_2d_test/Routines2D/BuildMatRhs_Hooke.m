% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [MAT, RHS] = BuildMatRhs_Hooke(NMat, dNdxMat, dNdyMat, ...
    MMat, dMdxMat, dMdyMat, xxInt, yyInt, wwInt, ffInt, Nodes, lambda1, ...
    lambda2, mu1, mu2, fx, fy, nQ, NodeNum, MAT, RHS)

% Integrate the weak form in the domain, corresponding to a Hooke solid.

% Initialization.
ElemMAT11 = zeros(8, 8);
ElemMAT12 = zeros(8, 8);
ElemMAT21 = zeros(8, 8);
ElemMAT22 = zeros(8, 8);
ElemRHS1  = zeros(8, 1);
ElemRHS2  = zeros(8, 1);

% Loop over integration points.
for i = 1 : nQ

    if ffInt(i) > 0
        lambda = lambda1;
        mu = mu1;
    elseif ffInt(i) < 0
        lambda = lambda2;
        mu = mu2;
    end

    N  = NMat(:, i);
    Nx = dNdxMat(:, i);
    Ny = dNdyMat(:, i);

    M  = MMat(:, i);
    Mx = dMdxMat(:, i);
    My = dMdyMat(:, i);

    NxNxT = [Nx; Mx] * [Nx; Mx]';
    NxNyT = [Nx; Mx] * [Ny; My]';
    NyNxT = [Ny; My] * [Nx; Mx]';
    NyNyT = [Ny; My] * [Ny; My]';

    % Compute element matrices.
    ElemMAT11 = ElemMAT11 + wwInt(i) * ( (lambda+2*mu) * NxNxT + mu * NyNyT );
    ElemMAT12 = ElemMAT12 + wwInt(i) * ( lambda * NxNyT + mu  * NyNxT );
    ElemMAT21 = ElemMAT21 + wwInt(i) * ( lambda * NyNxT + mu  * NxNyT );
    ElemMAT22 = ElemMAT22 + wwInt(i) * ( (lambda+2*mu) * NyNyT + mu * NxNxT );

    % Compute right hand side.
    ElemRHS1 = ElemRHS1 + wwInt(i) * ( [N; M] * fx );
    ElemRHS2 = ElemRHS2 + wwInt(i) * ( [N; M] * fy );

end

% Add element contribution to global matrix.
uuNodes = [Nodes          Nodes+2*NodeNum];
vvNodes = [Nodes+NodeNum  Nodes+3*NodeNum];

MAT(uuNodes, uuNodes) = MAT(uuNodes, uuNodes) + ElemMAT11;
MAT(uuNodes, vvNodes) = MAT(uuNodes, vvNodes) + ElemMAT12;
MAT(vvNodes, uuNodes) = MAT(vvNodes, uuNodes) + ElemMAT21;
MAT(vvNodes, vvNodes) = MAT(vvNodes, vvNodes) + ElemMAT22;

RHS(uuNodes) = RHS(uuNodes) + ElemRHS1;
RHS(vvNodes) = RHS(vvNodes) + ElemRHS2;

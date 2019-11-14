% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [N, dNdx, dNdy, xxInt, yyInt, wwInt] = ShapeFctsStrdFEM(...
    xxElem, yyElem, xxIntRef, yyIntRef, wwIntRef, nQ)

% Compute standard bi-linear shape functions at integration points in 
% the real element.

% Initialization.
N     = zeros(4, nQ);
dNdx  = zeros(4, nQ);
dNdy  = zeros(4, nQ);

xxInt = zeros(1, nQ); 
yyInt = zeros(1, nQ); 
wwInt = zeros(1, nQ); 

% Define N, dNdr, dNds at the integration points in the REFERENCE element.
N1 = 0.25*(1-xxIntRef).*(1-yyIntRef);
N2 = 0.25*(1+xxIntRef).*(1-yyIntRef);
N3 = 0.25*(1+xxIntRef).*(1+yyIntRef);
N4 = 0.25*(1-xxIntRef).*(1+yyIntRef);

dN1dr = -0.25*(1-yyIntRef);
dN2dr = +0.25*(1-yyIntRef);
dN3dr = +0.25*(1+yyIntRef);
dN4dr = -0.25*(1+yyIntRef);

dN1ds = -0.25*(1-xxIntRef);
dN2ds = -0.25*(1+xxIntRef);
dN3ds = +0.25*(1+xxIntRef);
dN4ds = +0.25*(1-xxIntRef);

NN    = [ N1  ;  N2  ;  N3  ;  N4  ];
dNNdr = [dN1dr; dN2dr; dN3dr; dN4dr];
dNNds = [dN1ds; dN2ds; dN3ds; dN4ds];

% Define N, dNdx, dNdy at the integration points in the REAL element.
N = NN;

for i = 1 : nQ
    dxdr = xxElem*dNNdr(:, i); % J11
    dydr = yyElem*dNNdr(:, i); % J12
    dxds = xxElem*dNNds(:, i); % J21
    dyds = yyElem*dNNds(:, i); % J22

    % det(J) = J11*J22 - J21*J12
    detJ = dxdr*dyds - dxds*dydr;
    wwInt(i) = wwIntRef(i) * detJ;

    dNdx(:, i) = ( dyds*dNNdr(:, i) - dydr*dNNds(:, i)) / detJ;
    dNdy(:, i) = (-dxds*dNNdr(:, i) + dxdr*dNNds(:, i)) / detJ;
end

% Map integration points from reference to real element.
xxInt = xxElem * NN;
yyInt = yyElem * NN;

% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function globalK = updateStiffness(globalK,omega,DOF,enrElem,pHDOF)
% This function updates only the non-constant portions of the global
% stiffness matrix assuming that the geometries represented by the
% inclusion and/or voids do not change as a function of time.

global CHI CONNEC CRACK DOMAIN MAT NODES PHI PSI XYZ ZETA

nXElem = DOMAIN(1);                                                         % Number of elements in the x-direction
nYElem = DOMAIN(2);                                                         % Number of elements in the y-direction
lXElem = DOMAIN(3);                                                         % Length of elements in the x-direction
Em     = MAT(1);                                                            % Young's modulus for the matrix
vm     = MAT(2);                                                            % Poisson's ratio for the matrix
Ef     = MAT(3);                                                            % Young's modulus for the fiber
vf     = MAT(4);                                                            % Poisson's ratio for the fiber
plane  = MAT(5);                                                            % Plane stress or plane strain
nCT    = size(PHI,2);                                                       % Number of crack tips 

% Update the original stiffness matrix to remove enriched variables
tradDOF = 2*(nXElem+1)*(nYElem+1);                                          % Number of traditional degrees of freedom
incDOF  = 2*max(NODES(:,30));                                               % Maximum inclusion degree of freedom
consDOF = max([tradDOF,incDOF,pHDOF]);                                      % Maximum constant degree of freedom

% Set nonconstant stiffness components equal to zero
globalK((consDOF+1):DOF,:) = 0;
globalK(:,(consDOF+1):DOF) = 0;

% Initialize updated stiffness matrix
enrichK = sparse(size(globalK,1),size(globalK,1));

% Create elastic constant matrix
if plane == 1                                                               % Plane stress
    h = MAT(6);                                                             % Plane stress thickness
    
    C1 = Em/(1-vm^2);                                                       % Constant for elastic constant matrix
    C2 = Em*vm/(1-vm^2);                                                    % Constant for elastic constant matrix
    C3 = Em/2/(1+vm);                                                       % Constant for elastic constant matrix
    Cm = h*[C1 C2  0;...
            C2 C1  0;...
             0  0 C3];
    C1 = Ef/(1-vf^2);                                                       % Constant for elastic constant matrix
    C2 = Ef*vf/(1-vf^2);                                                    % Constant for elastic constant matrix
    C3 = Ef/2/(1+vf);                                                       % Constant for elastic constant matrix
    Cf = h*[C1 C2  0;...
            C2 C1  0;...
             0  0 C3];
elseif plane == 2                                                           % Plane strain
    C1 = Em*(1-vm)/(1+vm)/(1-2*vm);                                         % Constant for elastic constant matrix
    C2 = Em*vm/(1+vm)/(1-2*vm);                                             % Constant for elastic constant matrix
    C3 = Em/2/(1+vm);                                                       % Constant for elastic constant matrix
    Cm  = [C1 C2  0;...
           C2 C1  0;...
            0  0 C3];
    C1 = Ef*(1-vf)/(1+vf)/(1-2*vf);                                         % Constant for elastic constant matrix
    C2 = Ef*vf/(1+vf)/(1-2*vf);                                             % Constant for elastic constant matrix
    C3 = Ef/2/(1+vf);                                                       % Constant for elastic constant matrix
    Cf  = [C1 C2  0;...
           C2 C1  0;...
            0  0 C3];              
end

m = size(CRACK,1);                                                          % Determine number of data points defining crack
if m > 0
    if nCT == 1
        xCT = CRACK(m,1);                                                   % X-coordinate of crack tip
        yCT = CRACK(m,2);                                                   % Y-coordinate of crack tip
    elseif nCT == 2
        xCT = [CRACK(1,1) CRACK(m,1)];                                      % X-coordinates of crack tips
        yCT = [CRACK(1,2) CRACK(m,2)];                                      % Y-coordinates of crack tips
    end
end

% Update the enriched elements
for iElem = 1:length(enrElem)
    nElem = enrElem(iElem);
    
    N1  = CONNEC(nElem,2);                                                  % Node 1 for current element
    N2  = CONNEC(nElem,3);                                                  % Node 2 for current element
    N3  = CONNEC(nElem,4);                                                  % Node 3 for current element
    N4  = CONNEC(nElem,5);                                                  % Node 4 for current element
    NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element
    CN  = [CHI(N1) CHI(N2) CHI(N3) CHI(N4)];                                % Nodal chi level set values            
    CTN = nnz(NN(:,4));                                                     % Number of nodes with crack tip enrichment    
    HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
    IEN = nnz(NN(:,30));                                                    % Number of inclusion nodes    
    NEN = HEN+CTN+IEN;                                                      % Number of enriched nodes
    
    localK = 0;                                                             % Initialize stiffness for current element
    local  = [N1*2-1 N1*2 N2*2-1 N2*2 N3*2-1 N3*2 N4*2-1 N4*2 NaN(1,96)];   % Traditional index locations
    iLoc   = 9;                                                             % Next index location
    
    X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2);         % Nodal x-coordinates
    Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3);         % Nodal y-coordinates
    xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                        % Nodal coordinate matrix
        
    if NEN == 4                                                             % Fully enriched element
        PN = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];                         % Nodal crack level set values
        ZN = [ZETA(N1) ZETA(N2) ZETA(N3) ZETA(N4)];                         % Nodal inclusion level set values
            
        if HEN == 4                                                         % Full Heaviside enrichment
            [gp,gw,J] = subDomain(3,CN,PN,ZN,xyz,0,0,[]);
        elseif CTN == 4                                                     % Full crack tip enrichmnet
            [gp,gw,J] = subDomain(7,CN,PN,ZN,xyz,1,0,[]);
        else                                                                % Full heaviside/crack tip enrichment
            [gp,gw,J] = subDomain(7,CN,PN,ZN,xyz,0,0,[]);
        end
        
        if min(CN) < 0
            xi      = gp(:,1);
            eta     = gp(:,2);
            N       = 1/4*[(1-xi).*(1-eta) (1+xi).*(1-eta) (1+xi).*(1+eta) (1-xi).*(1+eta)];
            chi     = N(:,1)*CN(1)+N(:,2)*CN(2)+N(:,3)*CN(3)+N(:,4)*CN(4);
            R       = find(chi <= 0);
            J(R,:)  = [];
            gp(R,:) = [];
            gw(R,:) = [];
        end
    else                                                                    % Partially enriched element
        if min(CN) < 0
            ZN        = [ZETA(N1) ZETA(N2) ZETA(N3) ZETA(N4)];              % Nodal inclusion level set values
            [gp,gw,J] = subDomain(3,CN,[],ZN,xyz,0,0,[]);
            xi        = gp(:,1);
            eta       = gp(:,2);
            N         = 1/4*[(1-xi).*(1-eta) (1+xi).*(1-eta) (1+xi).*(1+eta) (1-xi).*(1+eta)];
            chi       = N(:,1)*CN(1)+N(:,2)*CN(2)+N(:,3)*CN(3)+N(:,4)*CN(4);
            R         = find(chi <= 0);
            J(R,:)    = [];
            gp(R,:)   = [];
            gw(R,:)   = [];
        else
            [gp,gw] = gauss(6,'QUAD');
            J = [];
        end
    end
    
    for i = 1:length(gp)
        xi = gp(i,1); eta = gp(i,2);                                        % Gauss points
        W  = gw(i);                                                         % Gauss weights
        if isempty(J) == 0
            Ji   = [J(i,1) J(i,2);J(i,3) J(i,4)];                           % Jacobian of subdomain
            detJ = det(Ji);                                                 % Determinant of the Jacobian
        else
            detJ  = lXElem/2*lXElem/2;                                      % Determinant of the Jacobian
        end
        
        N  = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);...                         % Shape functions
                  (1+xi)*(1+eta);(1-xi)*(1+eta)];
        Nx = 2/lXElem*1/4*[-(1-eta);1-eta;1+eta;-(1+eta)];                  % Derivative of shape functions with respect to x
        Ny = 2/lXElem*1/4*[-(1-xi);-(1+xi);1+xi;1-xi];                      % Derivative of shape functions with respect to y
        
        X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2);     % Nodal x-coordinates
        Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3);     % Nodal y-coordinates
        
        Xgp = N(1)*X1+N(2)*X2+N(3)*X3+N(4)*X4;                              % The global X for the current gauss point
        Ygp = N(1)*Y1+N(2)*Y2+N(3)*Y3+N(4)*Y4;                              % The global Y for the current gauss point
        Zgp = N(1)*ZETA(N1)+N(2)*ZETA(N2)+N(3)*ZETA(N3)+N(4)*ZETA(N4);      % Material level set at current gauss point
        
        Benr = [];
        Bu = [Nx(1)   0   Nx(2)   0   Nx(3)   0   Nx(4)   0;...
                0   Ny(1)   0   Ny(2)   0   Ny(3)   0   Ny(4);...
              Ny(1) Nx(1) Ny(2) Nx(2) Ny(3) Nx(3) Ny(4) Nx(4)];
        
        index = 1;
        for iN = 1:4
            if NN(iN,2) ~= 0
                psi1 = PSI(N1);                                             % Psi level set value at node 1
                psi2 = PSI(N2);                                             % Psi level set value at node 2
                psi3 = PSI(N3);                                             % Psi level set value at node 3
                psi4 = PSI(N4);                                             % Psi level set value at node 4
                psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;             % Psi level set value at current gauss point
                
                Hgp = sign(psi);                                            % Heaviside value at current gauss point
                Hi  = NN(iN,3);                                             % Nodal Heaviside value
                H   = Hgp-Hi;                                               % Shifted Heaviside value
                
                Ba = [Nx(iN)*H     0;
                          0    Ny(iN)*H;
                      Ny(iN)*H Nx(iN)*H];
                Benr(:,index:(index+1)) = Ba;
                index = index+2;
                
                if (i == length(gp))
                    local(iLoc:(iLoc+1)) = [2*NN(iN,2)-1 2*NN(iN,2)];
                    iLoc = iLoc+2;
                end
            elseif NN(iN,4) ~= 0
                if nCT == 1
                    X     = Xgp-xCT;                                        % Horizontal distance from crack tip to gauss point
                    Y     = Ygp-yCT;                                        % Vertical distance from crack tip to gauss point
                    CCS   = [cos(omega) sin(omega);-sin(omega) cos(omega)];
                    XYloc = CCS*[X Y]';                                     % Change to crack tip coordinates
                    r     = sqrt(XYloc(1)^2+XYloc(2)^2);                    % Radius from crack tip to current gauss point
                    if r < 0.001*lXElem; r = 0.05*lXElem; end
                    theta = atan2(XYloc(2),XYloc(1));                       % Angle from crack tip to current gauss point
                elseif nCT == 2
                    X1  = Xgp-xCT(1);
                    Y1  = Ygp-yCT(1);
                    X2  = Xgp-xCT(2);
                    Y2  = Ygp-yCT(2);
                    CCS = [cos(omega(1)) sin(omega(1));-sin(omega(1)) cos(omega(1))];
                    XY1 = CCS*[X1 Y1]';
                    CCS = [cos(omega(2)) sin(omega(2));-sin(omega(2)) cos(omega(2))];
                    XY2 = CCS*[X2 Y2]';
                    r1  = sqrt(XY1(1)^2+XY1(2)^2);                          % Radius from crack tip to current gauss point
                    r2  = sqrt(XY2(1)^2+XY2(2)^2);
                    if r1 > r2
                        r = r2; theta = atan2(XY2(2),XY2(1));
                        CCS = [cos(omega(2)) sin(omega(2));-sin(omega(2)) cos(omega(2))];
                    elseif r2 > r1
                        r = r1; theta = atan2(XY1(2),XY1(1));
                        CCS = [cos(omega(1)) sin(omega(1));-sin(omega(1)) cos(omega(1))];
                    end
                    if r < 0.001*lXElem; r = 0.05*lXElem; end
                end
                
                c = 1/2/sqrt(r); ct = CCS(1,1); st = CCS(1,2);              % Constants
                
                if NN(iN,12) == 0                                           % Crack tip enrichment
                    a1gp = sqrt(r)*sin(theta/2);                            % Node 1 crack tip enrichment value
                    a2gp = sqrt(r)*cos(theta/2);                            % Node 2 crack tip enrichment value
                    a3gp = sqrt(r)*sin(theta)*sin(theta/2);                 % Node 3 crack tip enrichment value
                    a4gp = sqrt(r)*sin(theta)*cos(theta/2);                 % Node 4 crack tip enrichment value
                    
                    a1 = a1gp-NN(iN,5);                                     % Shifted alpha 1 enrichment value
                    a2 = a2gp-NN(iN,7);                                     % Shifted alpha 2 enrichment value
                    a3 = a3gp-NN(iN,9);                                     % Shifted alpha 3 enrichment value
                    a4 = a4gp-NN(iN,11);                                    % Shifted alpha 4 enrichment value
                    
                    % Derivative of crack tip enrichment functions with respect to X
                    Px = c*[-sin(theta/2)*ct              + cos(theta/2)*-st;...
                             cos(theta/2)*ct              + sin(theta/2)*-st;...
                            -sin(3*theta/2)*sin(theta)*ct + (sin(theta/2)+sin(3*theta/2)*cos(theta))*-st;...
                            -cos(3*theta/2)*sin(theta)*ct + (cos(theta/2)+cos(3*theta/2)*cos(theta))*-st];
                    
                    % Derivative of crack tip enrichment functions with respect to Y
                    Py = c*[-sin(theta/2)*st              + cos(theta/2)*ct;...
                             cos(theta/2)*st              + sin(theta/2)*ct;...
                            -sin(3*theta/2)*sin(theta)*st + (sin(theta/2)+sin(3*theta/2)*cos(theta))*ct;...
                            -cos(3*theta/2)*sin(theta)*st + (cos(theta/2)+cos(3*theta/2)*cos(theta))*ct];
                    
                    B1 = [Nx(iN)*a1+N(iN)*Px(1)          0;
                                   0            Ny(iN)*a1+N(iN)*Py(1);
                          Ny(iN)*a1+N(iN)*Py(1) Nx(iN)*a1+N(iN)*Px(1)];
                    
                    B2 = [Nx(iN)*a2+N(iN)*Px(2)          0;
                                   0            Ny(iN)*a2+N(iN)*Py(2);
                          Ny(iN)*a2+N(iN)*Py(2) Nx(iN)*a2+N(iN)*Px(2)];
                    
                    B3 = [Nx(iN)*a3+N(iN)*Px(3)          0;
                                   0            Ny(iN)*a3+N(iN)*Py(3);
                          Ny(iN)*a3+N(iN)*Py(3) Nx(iN)*a3+N(iN)*Px(3)];
                    
                    B4 = [Nx(iN)*a4+N(iN)*Px(4)          0;
                                   0            Ny(iN)*a4+N(iN)*Py(4);
                          Ny(iN)*a4+N(iN)*Py(4) Nx(iN)*a4+N(iN)*Px(4)];
                    
                    Bb = [B1 B2 B3 B4];
                    Benr(:,index:(index+7)) = Bb;
                    index = index+8;
                    
                    if (i == length(gp))
                        local(iLoc:(iLoc+7)) = [2*NN(iN,4)-1 2*NN(iN,4) 2*NN(iN,6)-1  2*NN(iN,6)...
                                                2*NN(iN,8)-1 2*NN(iN,8) 2*NN(iN,10)-1 2*NN(iN,10)];
                        iLoc = iLoc+8;
                    end
                else                                                        % Bimaterial crack tip enrichment
                    G1 = Em/2/(1+vm);                                       % Shear modulus for the matrix
                    G2 = Ef/2/(1+vf);                                       % Shear modulus for the fiber
                    
                    % Kosolov constant
                    if plane == 1                                           % Plane stress
                        k1 = (3-vm)/(1+vm);
                        k2 = (3-vf)/(1+vf);
                    elseif plane == 2                                       % Plane strain
                        k1 = 3-4*vm;
                        k2 = 3-4*vf;
                    end
                    
                    b = (G1*(k2-1)-G2*(k1-1))/(G1*(k2+1)+G2*(k1+1));        % Second Dundur's parameter
                    e = 1/(2*pi)*log((1-b)/(1+b));                          % Material constant
                    
                    % Common variables
                    sr  = sqrt(r);
                    st  = sin(theta);
                    st2 = sin(theta/2);
                    ct2 = cos(theta/2);
                    
                    a1gp  = sr*cos(e*log(r))*exp(-e*theta)*st2;             % Alpha 1 crack tip enrichment value
                    a2gp  = sr*cos(e*log(r))*exp(-e*theta)*ct2;             % Alpha 2 crack tip enrichment value
                    a3gp  = sr*cos(e*log(r))*exp(e*theta)*st2;              % Alpha 3 crack tip enrichment value
                    a4gp  = sr*cos(e*log(r))*exp(e*theta)*ct2;              % Alpha 4 crack tip enrichment value
                    a5gp  = sr*cos(e*log(r))*exp(e*theta)*st2*st;           % Alpha 5 crack tip enrichment value
                    a6gp  = sr*cos(e*log(r))*exp(e*theta)*ct2*st;           % Alpha 6 crack tip enrichment value
                    a7gp  = sr*sin(e*log(r))*exp(-e*theta)*st2;             % Alpha 7 crack tip enrichment value
                    a8gp  = sr*sin(e*log(r))*exp(-e*theta)*ct2;             % Alpha 8 crack tip enrichment value
                    a9gp  = sr*sin(e*log(r))*exp(e*theta)*st2;              % Alpha 9 crack tip enrichment value
                    a10gp = sr*sin(e*log(r))*exp(e*theta)*ct2;              % Alpha 10 crack tip enrichment value
                    a11gp = sr*sin(e*log(r))*exp(e*theta)*st2*st;           % Alpha 11 crack tip enrichment value
                    a12gp = sr*sin(e*log(r))*exp(e*theta)*ct2*st;           % Alpha 12 crack tip enrichment value
                    
                    a = [a1gp-NN(iN,5);...                                  % Shifted alpha 1 enrichment value
                         a2gp-NN(iN,7);...                                  % Shifted alpha 2 enrichment value
                         a3gp-NN(iN,9);...                                  % Shifted alpha 3 enrichment value
                         a4gp-NN(iN,11);...                                 % Shifted alpha 4 enrichment value
                         a5gp-NN(iN,13);...                                 % Shifted alpha 5 enrichment value
                         a6gp-NN(iN,15);...                                 % Shifted alpha 6 enrichment value
                         a7gp-NN(iN,17);...                                 % Shifted alpha 7 enrichment value
                         a8gp-NN(iN,19);...                                 % Shifted alpha 8 enrichment value
                         a9gp-NN(iN,21);...                                 % Shifted alpha 9 enrichment value
                         a10gp-NN(iN,23);...                                % Shifted alpha 10 enrichment value
                         a11gp-NN(iN,25);...                                % Shifted alpha 11 enrichment value
                         a12gp-NN(iN,27)];                                  % Shifted alpha 12 enrichment value
                    
                    % Derivative of bimaterial crack tip enrichment functions with respect to x1 (crack tip coordinate system)
                    px = c*[-exp(-e*theta)*sin(theta/2)*(cos(e*log(r))+2*e*sin(e*log(r)-theta));...
                             exp(-e*theta)*cos(theta/2)*(cos(e*log(r))-2*e*sin(e*log(r)-theta));...
                            -exp(e*theta)*sin(theta/2)*(cos(e*log(r))+2*e*sin(e*log(r)+theta));...
                             exp(e*theta)*cos(theta/2)*(cos(e*log(r))-2*e*sin(e*log(r)+theta));...
                            -exp(e*theta)*sin(theta)*(cos(e*log(r))*sin(3*theta/2)+2*e*sin(e*log(r)+theta)*sin(theta/2));...
                            -exp(e*theta)*sin(theta)*(cos(e*log(r))*cos(3*theta/2)+2*e*sin(e*log(r)+theta)*cos(theta/2));...
                             exp(-e*theta)*sin(theta/2)*(-sin(e*log(r))+2*e*cos(e*log(r)-theta));...
                             exp(-e*theta)*cos(theta/2)*(sin(e*log(r))+2*e*cos(e*log(r)-theta));...
                             exp(e*theta)*sin(theta/2)*(-sin(e*log(r))+2*e*cos(e*log(r)+theta));...
                             exp(e*theta)*cos(theta/2)*(sin(e*log(r))+2*e*cos(e*log(r)+theta));...
                             exp(e*theta)*sin(theta)*(-sin(e*log(r))*sin(3*theta/2)+2*e*cos(e*log(r)+theta)*sin(theta/2));...
                             exp(e*theta)*sin(theta)*(-sin(e*log(r))*cos(3*theta/2)+2*e*cos(e*log(r)+theta)*cos(theta/2))];
                    
                    % Derivative of bimaterial crack tip enrichment functions with respect to x2 (crack tip coordinate system)
                    py = c*[exp(-e*theta)*(cos(e*log(r))*cos(theta/2)-2*e*cos(e*log(r)-theta)*sin(theta/2));...
                            exp(-e*theta)*(cos(e*log(r))*sin(theta/2)-2*e*cos(e*log(r)-theta)*cos(theta/2));...
                            exp(e*theta)*(cos(e*log(r))*cos(theta/2)+2*e*cos(e*log(r)+theta)*sin(theta/2));...
                            exp(e*theta)*(cos(e*log(r))*sin(theta/2)+2*e*cos(e*log(r)+theta)*cos(theta/2));...
                            exp(e*theta)*(cos(e*log(r))*(sin(theta/2)+sin(3*theta/2)*cos(theta))+2*e*cos(e*log(r)+theta)*sin(theta/2)*sin(theta));...
                            exp(e*theta)*(cos(e*log(r))*(cos(theta/2)+cos(3*theta/2)*cos(theta))+2*e*cos(e*log(r)+theta)*cos(theta/2)*sin(theta));...
                            exp(-e*theta)*(sin(e*log(r))*cos(theta/2)-2*e*sin(e*log(r)-theta)*sin(theta/2));...
                            exp(-e*theta)*(sin(e*log(r))*sin(theta/2)-2*e*sin(e*log(r)-theta)*cos(theta/2));...
                            exp(e*theta)*(sin(e*log(r))*cos(theta/2)+2*e*sin(e*log(r)+theta)*sin(theta/2));...
                            exp(e*theta)*(sin(e*log(r))*sin(theta/2)+2*e*sin(e*log(r)+theta)*cos(theta/2));...
                            exp(e*theta)*(sin(e*log(r))*(sin(theta/2)+sin(3*theta/2)*cos(theta))+2*e*sin(e*log(r)+theta)*sin(theta/2)*sin(theta));...
                            exp(e*theta)*(sin(e*log(r))*(cos(theta/2)+cos(3*theta/2)*cos(theta))+2*e*sin(e*log(r)+theta)*cos(theta/2)*sin(theta))];
                    
                    % Derivative of bimaterial crack tip enrichment functions with respect to X
                    Px = [px(1)*ct+py(1)*-st;...
                          px(2)*ct+py(2)*-st;...
                          px(3)*ct+py(3)*-st;...
                          px(4)*ct+py(4)*-st;...
                          px(5)*ct+py(5)*-st;...
                          px(6)*ct+py(6)*-st;...
                          px(7)*ct+py(7)*-st;...
                          px(8)*ct+py(8)*-st;...
                          px(9)*ct+py(9)*-st;...
                          px(10)*ct+py(10)*-st;...
                          px(11)*ct+py(11)*-st;...
                          px(12)*ct+py(12)*-st];
                    
                    % Derivative of bimaterial crack tip enrichment functions with respect to Y
                    Py = [px(1)*st+py(1)*ct;...
                          px(2)*st+py(2)*ct;...
                          px(3)*st+py(3)*ct;...
                          px(4)*st+py(4)*ct;...
                          px(5)*st+py(5)*ct;...
                          px(6)*st+py(6)*ct;...
                          px(7)*st+py(7)*ct;...
                          px(8)*st+py(8)*ct;...
                          px(9)*st+py(9)*ct;...
                          px(10)*st+py(10)*ct;...
                          px(11)*st+py(11)*ct;...
                          px(12)*st+py(12)*ct];
                    
                    Bb = zeros(3,24);
                    for iB = 1:12
                        Balpha = [Nx(iN)*a(iB)+N(iN)*Px(iB)             0;...
                                              0             Nx(iN)*a(iB)+N(iN)*Py(iB);...
                                  Nx(iN)*a(iB)+N(iN)*Py(iB) Nx(iN)*a(iB)+N(iN)*Px(iB)];
                        Bb(:,(2*iB-1):2*iB) = Balpha;
                    end
                    Benr(:,index:(index+23)) = Bb;
                    index = index+24;

                    if i == length(gp)
                        for iB = 4:2:26
                            local(iLoc:(iLoc+1)) = [2*NN(iN,iB)-1 2*NN(iN,iB)];
                            iLoc = iLoc+2;
                        end
                    end
                end
            end
        end
        
        if Zgp > 0, C = Cm; else C = Cf; end
        B = [Bu Benr];
        localK = localK + W*B'*C*B*detJ;     
    end
    
    local(:,isnan(local(1,:))) = [];
    temp = find(local <= pHDOF);                                            % Find the constant degrees of freedom
    localK(temp,temp) = 0;                                                  % Set stiffness values associated with constant DOF to zero
    enrichK(local,local) = enrichK(local,local)+localK;
end

globalK = globalK+enrichK;                                                  % Add constant and updated stiffness together
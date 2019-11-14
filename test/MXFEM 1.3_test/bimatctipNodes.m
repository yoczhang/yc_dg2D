% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function bimatctipNodes
% This function assigns nodes enriched with the bimaterial crack tip 
% enrichment function values associated with the bimaterial near tip 
% asymptotic displacement field.

global CRACK DOMAIN MAT NODES PHI XYZ

lElem = DOMAIN(3);
Em    = MAT(1);                                                             % Young's modulus for the matrix
vm    = MAT(2);                                                             % Poisson's ratio for the matrix
Ef    = MAT(3);                                                             % Young's modulus for the fiber
vf    = MAT(4);                                                             % Poisson's ratio for the fiber
plane = MAT(5);                                                             % Plane stress or plane strain
G1    = Em/2/(1+vm);                                                        % Shear modulus for the matrix
G2    = Ef/2/(1+vf);                                                        % Shear modulus for the fiber

% Kosolov constant
if plane == 1                                                               % Plane stress
    k1 = (3-vm)/(1+vm);
    k2 = (3-vf)/(1+vf);
elseif plane == 2                                                           % Plane strain
    k1 = 3-4*vm;
    k2 = 3-4*vf;
end

b = (G1*(k2-1)-G2*(k1-1))/(G1*(k2+1)+G2*(k1+1));                            % Second Dundur's parameter
e = 1/(2*pi)*log((1-b)/(1+b));                                              % Material constant

iSeg = size(CRACK,1);                                                       % Number of crack segments
nCT  = size(PHI,2);                                                         % Number of crack tips

% Define coordinates of crack tip(s)
if nCT == 1
    xCT = CRACK(iSeg,1);                                                    % X-coordinate of crack tip
    yCT = CRACK(iSeg,2);                                                    % Y-coordinate of crack tip
elseif nCT == 2
    xCT = [CRACK(iSeg,1) CRACK(1,1)];                                       % X-coordinate of crack tip
    yCT = [CRACK(iSeg,2) CRACK(1,2)];                                       % Y-coordinate of crack tip    
end

for iNode = 1:size(NODES,1)
    if NODES(iNode,4) ~= 0
        XN    = XYZ(iNode,2);                                               % Nodal x-coordinate
        YN    = XYZ(iNode,3);                                               % Nodal y-coordinate
        X     = XN-xCT;                                                     % Horizontal distance from crack tip to current node
        Y     = YN-yCT;                                                     % Vertical distance from crack tip to current node
        for i = 1:length(X)
            rt = sqrt(X(i)^2+Y(i)^2);
            tt = atan2(Y(i),X(i));
            if i == 1; r = rt; theta = tt; end
            if i > 1; if rt < r, r = rt; theta = tt; end, end
        end
        
        if r < 0.001*lElem; r = 0.05*lElem; end                           
        
        % Common variables
        sr  = sqrt(r); 
        st  = sin(theta);
        st2 = sin(theta/2);
        ct2 = cos(theta/2);
        
        NODES(iNode,5)  = sr*cos(e*log(r))*exp(-e*theta)*st2;               % Alpha 1 crack tip enrichment value
        NODES(iNode,7)  = sr*cos(e*log(r))*exp(-e*theta)*ct2;               % Alpha 2 crack tip enrichment value
        NODES(iNode,9)  = sr*cos(e*log(r))*exp(e*theta)*st2;                % Alpha 3 crack tip enrichment value
        NODES(iNode,11) = sr*cos(e*log(r))*exp(e*theta)*ct2;                % Alpha 4 crack tip enrichment value
        NODES(iNode,13) = sr*cos(e*log(r))*exp(e*theta)*st2*st;             % Alpha 5 crack tip enrichment value
        NODES(iNode,15) = sr*cos(e*log(r))*exp(e*theta)*ct2*st;             % Alpha 6 crack tip enrichment value
        NODES(iNode,17) = sr*sin(e*log(r))*exp(-e*theta)*st2;               % Alpha 7 crack tip enrichment value
        NODES(iNode,19) = sr*sin(e*log(r))*exp(-e*theta)*ct2;               % Alpha 8 crack tip enrichment value
        NODES(iNode,21) = sr*sin(e*log(r))*exp(e*theta)*st2;                % Alpha 9 crack tip enrichment value
        NODES(iNode,23) = sr*sin(e*log(r))*exp(e*theta)*ct2;                % Alpha 10 crack tip enrichment value
        NODES(iNode,25) = sr*sin(e*log(r))*exp(e*theta)*st2*st;             % Alpha 11 crack tip enrichment value
        NODES(iNode,27) = sr*sin(e*log(r))*exp(e*theta)*ct2*st;             % Alpha 12 crack tip enrichment value
    end
end
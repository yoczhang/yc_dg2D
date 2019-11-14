% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function connectivity
% This function defines connectivity information about the problem.  This
% includes defining node numbers for the problem, defining the XYZ 
% coordinates of the node numbers, creating a connectivity matrix, and
% creating matricies storing the X and Y coordinates for each element.

global CONNEC DOMAIN NODES XYZ

nXElem = DOMAIN(1);                                                         % Number of elements in the x-direction
nYElem = DOMAIN(2);                                                         % Number of elements in the y-direction
lXElem = DOMAIN(3);                                                         % Elemental length in the x-direction
lYElem = DOMAIN(4);                                                         % Elemental length in the y-direction
nNode  = (nXElem+1)*(nYElem+1);                                             % The number of nodes in the domain
NODES  = zeros(nNode,31);                                                   % Initialize the matrix storing nodal information
XYZ    = zeros(nNode,3);                                                    % Initialize the matrix storing nodal position

% Create global node numbering (NN), global xyz coordinates (XYZ), and index of enriched nodes (NODES)
% XYZ = [NodeNumber,X-Coordinate,Y-Coordinate,Z-Coordinate;...]
% NODES = [NodeNumber,EnrichedNodeNumber,HeavisideValue;...]
nNode = 1;
NN = zeros(nXElem+1,nYElem+1);
for iYNode = 1:(nYElem+1)
    for iXNode = 1:(nXElem+1)
        NN(iXNode,iYNode) = nNode;
        XYZ(nNode,:) = [nNode (iXNode-1)*lXElem (iYNode-1)*lYElem];
        NODES(nNode,1) = nNode;
        nNode = nNode+1;
    end
end

% Create a global connectivity matrix (CONNEC) and elemental coordinate matricies (Ex, Ey)
% CONNEC  = [ElementNumber,LocalNode1,LocalNode2,LocalNode3,LocalNode4,InclusionElementNum]
nElem = 1;
CONNEC = zeros(nXElem*nYElem,5);
for iYElem = 1:nYElem
    for iXElem = 1:nXElem
        N1 = NN(iXElem,iYElem);
        N2 = NN(iXElem+1,iYElem);
        N3 = NN(iXElem+1,iYElem+1);
        N4 = NN(iXElem,iYElem+1);
        CONNEC(nElem,1:5) = [nElem N1 N2 N3 N4];
        nElem = nElem+1;
    end
end
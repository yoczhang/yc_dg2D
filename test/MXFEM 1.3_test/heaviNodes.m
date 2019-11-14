% Written By: Matthew Jon Pais, University of Florida (2009)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function heaviNodes
% This function assigns nodes enriched with the Heaviside function as
% either above (+1) or below (-1) the crack.

global NODES PSI

for iNode = 1:size(NODES,1)
    if NODES(iNode,2) ~= 0
        NODES(iNode,3) = sign(PSI(iNode));
    end
end
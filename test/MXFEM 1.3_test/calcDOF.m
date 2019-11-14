% Written By: Matthew Jon Pais, University of Florida (2009)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [DOF,DISP] = calcDOF
% This function calculates the total degrees of freedom consisting of
% classical and enriched DOF.

global NODES

heavDOF = 2*nnz(NODES(:,2));                                                % Define the number of Heavi DOF
ctipDOF = 8*nnz(NODES(:,4));                                                % Define the number of crack tip DOF
bmctDOF = 24*nnz(NODES(:,12));                                              % Define the number of bimaterial crack tip DOF
if heavDOF > 0, heaviNodes;     end                                         % Define heaviside nodal enrichment values
if ctipDOF > 0, ctipNodes;      end                                         % Define crack tip nodal enrichment values
if bmctDOF > 0, bimatctipNodes; end                                         % Define bimaterial crack tip nodal enrichment values
DOF     = 2*max(max(NODES));                                                % Total number of degrees of freedom
DISP    = sparse(DOF,1);                                                    % Initialize displacement vector
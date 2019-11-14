function [G_fault, originalFracture] = getFaultsMeshFromMRST_exam1()
%
%
%   ComplexFaultNetwork
%   MRST 
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Just copy form ycComplexFaultsWells.m
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   YcZhang Apr.25.2019
%   Last modified Apr.25.2019
%

addpath(genpath([ROOTDIR,'modules/upr/'])); 

% faults = {
%     [0,1/2; 1/2,1/2;];
%     [1/2,1/2; 1, 1/2];
%     };
faults = {
    [0.008, 0.5; 1, 0.5;];
    };
faults = faults'; % notice that this cell-type varible, must be a ROW-vector.

% % [faults, ~, ~] = splitAtInt(faults, {});

originalFracture = faults;

%%
if (1)
    %% Set pebiGrid Parameters:
    %We now use the other wrapper function to create a PEBI-grid using distmesh:
    gS  = 1/16; % The grid size
    wGf = 0.8;  % The relative size of well cells compared to gS
    fGf = 0.3;  % The relative size of fault cells compared to gS
    eps = 1/(20); % This parameter defines the refinement around the wells. The            
    % cell size increase exponentially from the wells with            
    % exp(-(distance from well)/eps);

    %% Generate grid
        %> distmesh will most likely not converge in the maximum number 
        %> of iterations. This is usually not a problem, since the grid most 
        %> likely is good before the convergence requirement is met.
    %------------------------

    %--- only fault
    G_fault = pebiGrid(gS, [1, 1], 'faultlines', faults, 'faultGridFactor',fGf*0.4, ...
        'faultRefinement',true, 'faultEps',eps);

    %---
    %- firstly, we remove the short edges.
    G_fault = removeShortEdges(G_fault, 8e-8);
    %- geometry
    G_fault = computeGeometry(G_fault);
    %- save 
%     save G_faultbackup G_fault
    
    % Plot pebiGrid
    figure(); hold on
    plotGrid(G_fault,'facecolor','none')
    centF = G_fault.faces.centroids(G_fault.faces.tag,:);
    plot(centF(:,1), centF(:,2),'.','markersize',10)
    axis equal tight
    title('pebiGrid(...)')
end % if (0)


end % function
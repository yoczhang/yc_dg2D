function [G_fault, originalFracture] = getFaultsMeshFromMRST_case2()
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

faults = {
    [0.60, 0.72; 0.68, 0.74; 0.75, 0.76; 0.83, 0.76; 0.92, 0.75; 0.95, 0.78;] - [0.05, 0.1];
    [0.85, 0.60; 0.88, 0.65; 0.89, 0.68; 0.90, 0.70; 0.91, 0.73; 0.92, 0.75;] - [0.05, 0.1];
    [0.55,0.55; 0.60,0.60; 0.65,0.62; 0.7,0.65; 0.85,0.85; 0.9,0.87; ] - [0.05, 0.1];
    [0.13,0.24; 0.16,0.26; 0.18,0.28; 0.20,0.3; 0.22,0.32; 0.24,0.38; 0.26,0.6;] - [0, 0.13];
    [0.2,0.6; 0.22,0.56; 0.24,0.52; 0.26,0.5; 0.28,0.46; 0.32,0.43; 0.36,0.34; 0.38,0.28;] - [0, 0.13];
    };
faults = faults'; % notice that this cell-type varible, must be a ROW-vector.

% % [faults, ~, ~] = splitAtInt(faults, {});

originalFracture = faults;

%%
if (1)
    %% Set pebiGrid Parameters:
    %We now use the other wrapper function to create a PEBI-grid using distmesh:
    gS  = 1/32; % The grid size
    wGf = 0.5;  % The relative size of well cells compared to gS
    fGf = 0.5;  % The relative size of fault cells compared to gS
    eps = 1/(10); % This parameter defines the refinement around the wells. The            
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
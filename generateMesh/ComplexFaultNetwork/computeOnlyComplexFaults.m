%
%   ComplexFaultNetwork
%   MRST 
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Just copy form ycComplexFaultsWells.m
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   YcZhang Feb.21.2019
%   Last modified Feb.21.2019
%

close all
clearvars
clc
addpath(genpath([ROOTDIR,'modules/upr/'])); 
%% load Geometry.
% % pth = fullfile(mrstPath('upr'), 'datasets', 'gridBranets.mat');
% % load(pth);

% % aa = 1:2:length(fault);
% % ycFault = cell(length(aa),1);
% % for n = 1:length(ycFault)
% %     data_fault = fault{aa(n)};
% % %     max_x = max(abs(data_fault(:,1)));
% % %     min_x = min(abs(data_fault(:,1)));
% % %     max_y = max(abs(data_fault(:,2)));
% % %     min_y = min(abs(data_fault(:,2)));
% % %     
% % %     % normalize the data to [0,1]x[0,1]
% % %     data_fault(:,1) = (abs(data_fault(:,1)) - min_x)/(max_x-min_x);
% % %     data_fault(:,2) = (abs(data_fault(:,2)+max_y) - min_y)/(max_y-min_y);
% % %     ycFault{n} = data_fault;
% %     test_fault = {data_fault};
% %     color = get(gca,'ColorOrder');
% %     plotLinePath(test_fault,'color',color(2,:));
% %     hold on
% % end 
scale1 = 10;
faults = {%[0.1,0.1; 0.2, 0.25; 0.25, 0.25; 0.3, 0.27; 0.325, 0.32; 0.33, 0.37; 0.37, 0.45];
    %[0.07, 0.4; 0.13, 0.37; 0.16, 0.29; 0.2, 0.25];
    %[0.35, 0.795; 0.4, 0.79; 0.5, 0.77; 0.55, 0.75; 0.6, 0.741; 0.68, 0.74; 0.75, 0.76; 0.83, 0.76; 0.92, 0.75; 0.95, 0.78; 0.99, 0.8;];
    %[0.68, 0.74; 0.70, 0.738; 0.72, 0.737; 0.73, 0.735; 0.74, 0.73; 0.75, 0.725; 0.76, 0.72; 
    %0.77, 0.714; 0.78, 0.705; 0.795, 0.705; 0.85, 0.7; 0.89, 0.68;];
    [0.80, 0.60; 0.88, 0.65; 0.89, 0.68; 0.90, 0.70; 0.91, 0.73; 0.92, 0.75;];
    [0.50, 0.711; 0.54, 0.716; 0.60, 0.72; 0.68, 0.74;];
    %[0.2, 0.40; 0.23, 0.47; 0.235, 0.48; 0.25, 0.5; 0.28, 0.54; 0.3, 0.57; 0.34, 0.6; 0.36, 0.7; 0.415, 0.80; 0.42, 0.83; 0.43, 0.85];
    %[0.36, 0.7; 0.415, 0.80; 0.42, 0.83; 0.43, 0.85];
    [0.13,0.24; 0.16,0.26; 0.18,0.28; 0.20,0.3; 0.22,0.32; .24,0.38; 0.26,0.6;];
    [0.2,0.6; 0.22,0.56; 0.24,0.52; 0.26,0.5; 0.28,0.46; 0.32,0.40; 0.36,0.34; 0.38,0.28; ...
    0.4,0.22; 0.405,0.21; 0.5,0.1; 0.53,0.06];
    %[0.5,0.53; 0.55,0.55; 0.60,0.60; 0.65,0.62; 0.7,0.65; 0.85,0.85; ...
    %0.9,0.87; 0.92,0.86; ];
    %[0.4,0.1; 0.55,0.15; 0.65,0.26; 0.75,0.3;];
    %[0.6,0.45; 0.62,0.35; 0.66,0.25; 0.68,0.15; 0.7,0.1];
    }; % [x-coord,y-coord]

% % faults = {
% %     [0.80, 0.60; 0.88, 0.65; 0.89, 0.68; 0.90, 0.70; 0.91, 0.73; 0.92, 0.75;];
% %     };
faults = faults'; % notice that this cell-type varible, must be a ROW-vector.



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
    save G_onefault G_fault
    
    % Plot pebiGrid
    figure(); hold on
    plotGrid(G_fault,'facecolor','none')
    centF = G_fault.faces.centroids(G_fault.faces.tag,:);
    plot(centF(:,1), centF(:,2),'.','markersize',10)
    axis equal tight
    title('pebiGrid(...)')
end % if (0)

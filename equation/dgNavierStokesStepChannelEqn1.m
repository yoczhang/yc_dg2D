function [sysErr,sysTime] = dgNavierStokesStepChannelEqn1(pde,option,varargin)
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       Using the Saclar-Stokes, 
%       the dgNavierStokes2 using the Tensor-Stokes.
%   %---------------------------------------------
%
%   dgNavierStokesEqn solve the stokes equation by DG methods
%
%
%	YcZhang 17/10/2017
%
%   Last modified 21/108/2017
%

%% Default setting of mesh and pde data
if ~exist('option','var') 
    option = []; 
end
if ~exist('pde','var')
    pde = navierstokesData(0,4); % default data
end


%% Parameters
option = dgOption(option);
maxIt = 1;


%% Initialize err
sysErr = zeros(maxIt,1); % not used
sysTime = zeros(maxIt,1); % not used

h = zeros(maxIt,1);
% Nnodes = zeros(maxIt,1);
% Ndof = zeros(maxIt,1);
dgNavierStokesStepChannel_output = cell(maxIt,1);

%% set path and file name
load('setpath_pwd.mat')
dgfunc_name = 'dgNavierStokesStepChannel_';

%>>>>>>>>>  creat log file >>>>>>>>
date = datestr(now,31); 
    %> capture the now time, 31 stands for the scheme: 2000-03-01 15:45:17
logFilename=[setpath_pwd,'/logs/',dgfunc_name,date,'_log.txt'];
diary(logFilename);
diary on; % begin diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<

disp('********************  dg NavierStokesStepChannel  ********************')
disp('------------------------------------------------------')
disp('NavierStokesStepChannel Equation:')
disp(['   k_11 = ',func2str( pde.k_11)])
disp(['   k_22 = ',func2str( pde.k_22)])
disp(['   u1 = ',func2str(pde.u1)])
disp(['   u2 = ',func2str(pde.u2)])
disp('DG parameters: ')
disp(['   basesType_u: ', option.basesType_u])
disp(['   basesType_p: ', option.basesType_p])
disp(['   penalty pars: p_epsilon: ', num2str(option.p_epsilon)])
disp(['   penalty pars: p_sigma: ', num2str(option.p_sigma)])
disp(['   penalty pars: p_beta: ', num2str(option.p_beta)])
disp(['   maxIt: ', num2str(maxIt)])

%% DG methods
disp('------------------------------------------------------')
disp('Solving dgNavierStokesStepChannel: ')
disp('------------------------------------------------------')
% format long e
for n = 1:maxIt
    n_str = num2str(n);

    %---- mesh 5
    %----------------- NS Step Channel ------------------
    file=fopen('NS_StepChannel_nodes.dat');points=fscanf(file,'%f');fclose(file);
    nppp=size(points);npp=nppp(1)/3;node=zeros(npp,2);
    for i=1:npp
        node(i,1)=points(3*i-1);node(i,2)=points(3*i);
    end
    file=fopen('NS_StepChannel_elems.dat');connecticy=fscanf(file,'%f');fclose(file);
    neee=size(connecticy);nee=neee(1)/9;elem=zeros(nee,3);
    for i=1:nee
        elem(i,1)=connecticy(9*i-6)+1;elem(i,2)=connecticy(9*i-3)+1;elem(i,3)=connecticy(9*i)+1;
    end
    %-----------------------------------------------------------
    
    
    disp('')
    n_th_cycle = ['the ',n_str,'-th cycle(i.e. the ', n_str, '-th mesh)'];
    disp(n_th_cycle)
    
    
    %% get mesh information
    meshInfo = polyMeshAuxStructure(node, elem);
	%plotPolyMsh(meshInfo)
    patchPlotMesh(node, elem);

    %% solve equations
    solve_t0 = cputime;
    [Uh, sysInfo]= dgNavierStokesSolve1(meshInfo,pde,option);
    sysInfo.SoverTime = cputime - solve_t0;
    disp(['the solve-func time: ',num2str(sysInfo.SoverTime)])
    
    %% compute the err
    dof_u1 = sysInfo.dof_u1;
    dof_p = sysInfo.dof_p;
    uh1 = Uh(1:dof_u1);
    uh2 = Uh(1+dof_u1:2*dof_u1);
    ph = Uh(1+2*dof_u1:end-1);
    
    %% other options
%     h(n) = sum(meshInfo.diameters)/length(elem);
%     h(n) = sum(dia)/length(elem);
    h(n) = 1./(sqrt(meshInfo.Nnodes)-1);
%     h_2 = sum(meshInfo.diameters)/length(elem);

    %% creat the structure of save mat file
    n_cycle = ['cycle_',n_str,'_'];
    dgNavierStokesStepChannel_output{n,1} = ...
        struct( ...
        [n_cycle,'basestype_u'], option.basesType_u, ...
        [n_cycle,'basestype_p'], option.basesType_p, ...
        [n_cycle,'dof_u1'], dof_u1, ...
        [n_cycle,'dof_p'], dof_p, ...
        [n_cycle,'uh1'], uh1, ...
        [n_cycle,'uh2'], uh2, ...
        [n_cycle,'ph'], ph, ...
        [n_cycle,'h'], h(n) ...
        );

    disp('------------------------------------------------------')
end % for n

% save result to mat file
saveFilename=[setpath_pwd,'/outputmat/',dgfunc_name,date,'.mat'];
save(saveFilename, 'dgNavierStokesStepChannel_output');


%% plot the results
degree_u = basesType2degreek(option.basesType_u);
degree_p = basesType2degreek(option.basesType_p);
Nelems = meshInfo.Nelems;
Nbases_u = nchoosek(degree_u+2,2);
Nbases_p = nchoosek(degree_p+2,2);
triCount = [1, 2, 3];
DGM = [];
DGT = [];
DGuh1 = [];
DGuh2 = [];
DGph = [];
for ii = 1:Nelems
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    singleElem = meshInfo.elem{ii,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges on singleElem.
    
    concavePointElem_ii = meshInfo.concavePointElem(ii); % [1 x 1], the concavePointElem of ii-the elem.
    coordv = meshInfo.node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    centroidElem = meshInfo.centroidElem(ii,:); % [1 x 2], the centroid(xing xin) (x-coord, y-coord) of element.
    
    coordTri0Elem = getcoordTri0Elem(singleNE, concavePointElem_ii, centroidElem, coordv); 
        %> if singleNE==3, coordTri0Elem, [3 x 2]. else coordTri0Elem, [3*singleNE x 2].
        
    elem_xT = centroidElem(1);
    elem_yT = centroidElem(2); 
    elem_hT = meshInfo.hElem(ii);
    
    %% Part II
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %-- on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.

        % get the bases values on quad points
        [trialPb_u, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, coordTri_nt(:,1), coordTri_nt(:,2), degree_u);
        [trialPb_p, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, coordTri_nt(:,1), coordTri_nt(:,2), degree_p);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases_trial].
        value_uh1_onTri = trialPb_u*uh1((ii-1)*Nbases_u+1:ii*Nbases_u);
        value_uh2_onTri = trialPb_u*uh2((ii-1)*Nbases_u+1:ii*Nbases_u);
        value_ph_onTri = trialPb_p*ph((ii-1)*Nbases_p+1:ii*Nbases_p);
        
        %--- temp 
        DGM_temp = [DGM; coordTri_nt];
        DGT_temp = [DGT; triCount];
        DGuh1_temp = [DGuh1; value_uh1_onTri];
        DGuh2_temp = [DGuh2; value_uh2_onTri];
        DGph_temp = [DGph; value_ph_onTri];
        
        %---- 
        DGM = DGM_temp;
        DGT = DGT_temp;
        DGuh1 = DGuh1_temp;
        DGuh2 = DGuh2_temp;
        DGph = DGph_temp;
        
        %---- update triCount
        triCount = triCount + 3;
    end % for nt
    %<<-- End Part II --------------------------------------------------------------------------------------
end % for ii

%---- plot speed
figure
trisurf(DGT,DGM(:,1),DGM(:,2),sqrt(DGuh1.^2+DGuh2.^2));
shading interp, xlabel('x'), ylabel('y'), colorbar , axis equal, axis off%, view(2)
title('Norm_{speed}');

%---- plot pressure
ph_min=min(DGph);
DGph=DGph+abs(ph_min)+2;
figure
trisurf(DGT,DGM(:,1),DGM(:,2),DGph);
shading interp, xlabel('x'), ylabel('y'), colorbar , axis equal, axis off%, view(2)
title('Pressure');


disp('------------------------------------------------------')

%% close the diary
%>>>>>>>>>  close the diary >>>>>>
diary off; % close the diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<
end % function dgPoissonEqn
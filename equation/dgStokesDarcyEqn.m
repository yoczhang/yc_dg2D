function [sysErr,sysTime] = dgStokesDarcyEqn(pde,option,varargin)
%
%   dgTensorStokesEqn solve the Tensor-form stokes equation by DG methods
%
%
%	YcZhang 25/8/2017
%
%   Last modified 26/8/2017
%

%% Default setting of mesh and pde data
if ~exist('option','var') 
    option = []; 
end
if ~exist('pde','var')
    pde = StokesDarcyData(0,7); % default data
end


%% Parameters
option = dgOption(option);
maxIt = option.maxIt;


%% Initialize err
sysErr = zeros(maxIt,1); % not used
sysTime = zeros(maxIt,1); % not used

phih_L2_error = zeros(maxIt,1); phih_H1_error = zeros(maxIt,1);
uh_L2_error = zeros(maxIt,1); uh_H1_error = zeros(maxIt,1);
ph_L2_error = zeros(maxIt,1);

phihL2rate = zeros(maxIt,1); phihH1rate = zeros(maxIt,1);
uhL2rate = zeros(maxIt,1); uhH1rate = zeros(maxIt,1);
phL2rate = zeros(maxIt,1);
h_S = zeros(maxIt,1); h_D = zeros(maxIt,1);
% Nnodes = zeros(maxIt,1);
% Ndof = zeros(maxIt,1);
dgStokesDarcy_output = cell(maxIt,1);

%% set path and file name
load('setpath_pwd.mat')
dgfunc_name = 'dgStokesDarcy_';

%>>>>>>>>>  creat log file >>>>>>>>
date = datestr(now,31); 
    %> capture the now time, 31 stands for the scheme: 2000-03-01 15:45:17
logFilename=[setpath_pwd,'/logs/',dgfunc_name,date,'_log.txt'];
diary(logFilename);
diary on; % begin diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<


disp('********************  dg Stokes Darcy  ********************')
disp('------------------------------------------------------')
disp('Stokes Darcy Equation:')
disp(['   phi = ',func2str(pde.phi)])
disp(['   u1 = ',func2str(pde.u1)])
disp(['   u2 = ',func2str(pde.u2)])
disp(['   p = ',func2str(pde.p)])
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
disp('Stokes eqn fluid kinematic viscosity:')
disp(['   mu = ', num2str(pde.mu)])
disp('Darcy eqn PermeabilityCoefficient:')
disp(['   PermeabiltiyCoeff (i.e. pde.K) = ', num2str(pde.K)])
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
disp('DG parameters: ')
disp(['   basesType_phi: ', option.basesType_Pp])
disp(['   basesType_u: ', option.basesType_Fu])
disp(['   basesType_p: ', option.basesType_Fp])
disp(['   penalty pars: p_epsilon: ', num2str(option.p_epsilon)])
disp(['   penalty pars: p_sigma: ', num2str(option.p_sigma)])
disp(['   penalty pars: Darcy_p_sigma: ', num2str(option.Darcy_p_sigma)])
disp(['   penalty pars: p_beta: ', num2str(option.p_beta)])
disp(['   maxIt: ', num2str(maxIt)])

%% DG methods
disp('------------------------------------------------------')
disp('Solving dgStokesDarcy: ')
disp('------------------------------------------------------')
% format long e
for n = 1:maxIt
    close all
    n_str = num2str(n);
    
    disp('')
    n_th_cycle = ['the ',n_str,'-th cycle(i.e. the ', n_str, '-th mesh)'];
    disp(n_th_cycle)
    
    
    if option.usingTensorStokes == 1 && strcmpi(pde.pdecase,'case1')
        %---- mesh 1
        %-------------------- Tri mesh ---------------------
        %--- get Stokes mesh information
        left=0; right=1; bottom=0; top=1;
        h_x = 1/2^(n+1); 
        %h_x = 1/4;
        h_y = h_x;
        h_partition=[h_x,h_y];
        [Snode, Selem] = generate_Tri_P_T(left,right,bottom,top,h_partition);
        
        %--- get Darcy mesh information
        left=0; right=1; bottom=-1; top=0;
        h_x = 1/2^(n+1); 
        %h_x = 1/4;
        h_y = h_x;
        h_partition=[h_x,h_y];
        [Dnode, Delem] = generate_Tri_P_T(left,right,bottom,top,h_partition);
        %-----------------------------------------------------
        
        %---- mesh 2
        %-------------------- NonRectangular mesh ---------------------
%         %--- get Stokes mesh information
%         %load('trapezoidal_2')
%         load('non-rectangular_2')
%         Snode = vertices'; Selem = elements';
%         %--- get Darcy mesh information
%         Dnode = vertices'; Delem = elements';
        %--------------------------------------------------------------------------
        
    elseif option.usingTensorStokes == 1 && strcmpi(pde.pdecase,'case2')
        %% get Stokes mesh information
        mesh_name = ['case2','_Stokes_Polygon_',n_str,'_1'];
        load(mesh_name);
        Snode = vertices'; Selem = elements';
        
        %% get Darcy mesh information
        %n_str = num2str(n+1);
        mesh_name = ['case2','_Darcy_Polygon_',n_str,'_1'];
        load(mesh_name);
        Dnode = vertices'; Delem = elements';
    elseif option.usingTensorStokes == 0 && ...
            (strcmpi(pde.pdecase,'case3') || strcmpi(pde.pdecase,'case6') ...
            || strcmpi(pde.pdecase,'case7') || strcmpi(pde.pdecase,'case8'))
        %---Sdomain: [0,1]x[1,2]; Ddomain: [0,1]x[0,1].
        %---- mesh 1
        %-------------------- Poly mesh ---------------------
%         %--- get Stokes mesh information
%         mesh_name = ['case2','_Stokes_Polygon_',n_str,'_1'];
%         load(mesh_name);
%         Snode = vertices'; Selem = elements';
%         
%         %--- get Darcy mesh information
%         %n_str = num2str(n+1);
%         mesh_name = ['case2','_Darcy_Polygon_',n_str,'_1'];
%         load(mesh_name);
%         Dnode = vertices'; Delem = elements';
        %-----------------------------------------------------
        
        %---- mesh 2
        %-------------------- Tri mesh ---------------------
        %--- get Stokes mesh information
        left=0; right=1; bottom=1; top=2;
        h_x = 1/2^(n+1); 
        %h_x = 1/4;
        h_y = h_x;
        h_partition=[h_x,h_y];
        [Snode, Selem] = generate_Tri_P_T(left,right,bottom,top,h_partition);
        
        %--- get Darcy mesh information
        left=0; right=1; bottom=0; top=1;
        h_x = 1/2^(n+1); 
        %h_x = 1/4;
        h_y = h_x;
        h_partition=[h_x,h_y];
        [Dnode, Delem] = generate_Tri_P_T(left,right,bottom,top,h_partition);
        %-----------------------------------------------------
    elseif option.usingTensorStokes == 0 && strcmpi(pde.pdecase,'case4')
        n_str_x = num2str(2^(n+1));
        n_str_y = num2str(2^n);
        mesh_name = ['S_quadmesh_',n_str_x,'times',n_str_y,'_[0_1]_[1-2_1]'];
        load(mesh_name);
        Snode = vertices'; Selem = elements';
        
        mesh_name = ['D_quadmesh_',n_str_x,'times',n_str_y,'_[0_1]_[0_1-2]'];
        load(mesh_name);
        Dnode = vertices'; Delem = elements';
    elseif option.usingTensorStokes == 0 && strcmpi(pde.pdecase,'case5')
        n_str_x = num2str(2^(n+1));
        n_str_y = num2str(2^n);
        mesh_name = ['S_quadmesh_',n_str_x,'times',n_str_y,'_[0_1]_[1-2_1]'];
        load(mesh_name);
        Snode = vertices'; Selem = elements';
        
        mesh_name = ['D_quadmesh_',n_str_x,'times',n_str_y,'_[0_1]_[0_1-2]'];
        load(mesh_name);
        Dnode = vertices'; Delem = elements';
    end % if
    
    %------------- Stokes meshInfo ----------------------
    StokesmeshInfo = polyMeshAuxStructure(Snode, Selem);
    %DiffusivityCoeffs = ones(size(StokesmeshInfo.elem,1),1);
    %StokesmeshInfo.DiffusivityCoeffs = DiffusivityCoeffs;
    %----------------------------------------------------------
    %-------------- Darcy meshInfo ----------------------
    DarcymeshInfo = polyMeshAuxStructure(Dnode, Delem);
    %--- get the enriched elems in the subfunction:
	%enrichedElem = getEnrichedElemForDarcyDomain(DarcymeshInfo.node, DarcymeshInfo.elem, 3);
    %DiffusivityCoeffs = ones(size(DarcymeshInfo.elem,1),1);
	%DiffusivityCoeffs(enrichedElem) = pde.DiffusivityCoeff;
	%DarcymeshInfo.DiffusivityCoeffs = DiffusivityCoeffs;
    %--- the Darcy permeability matrix coefficient.
    PermeabilityCoeffs = pde.K*ones(size(DarcymeshInfo.elem,1),1);
    DarcymeshInfo.PermeabilityCoeffs = PermeabilityCoeffs;
    %----------------------------------------------------------
    
    %interfacemeshInfo=test_coupled_interface_1(StokesmeshInfo, DarcymeshInfo, 'xcoord');
%     plotPolyMsh(StokesmeshInfo)
    %patchPlotMesh(StokesmeshInfo.node, StokesmeshInfo.elem)
%     plotPolyMsh(DarcymeshInfo)
    %patchPlotMesh(DarcymeshInfo.node, DarcymeshInfo.elem)

    %% solve equations
    solve_t0 = cputime;
    if option.usingTensorStokes == 1
        [Uh, sysInfo] = dgStokesDarcySolve_usingTensorStokes(StokesmeshInfo,DarcymeshInfo,pde,option);
    else
        [Uh, sysInfo] = dgStokesDarcySolve_usingScalarStokes(StokesmeshInfo,DarcymeshInfo,pde,option);
    end
    sysInfo.SoverTime = cputime - solve_t0;
    disp(['the solve-func time: ',num2str(sysInfo.SoverTime)])
    
    %% compute the err
    dof_phi = sysInfo.dof_phi;
    dof_u1 = sysInfo.dof_u1;
    dof_p = sysInfo.dof_p;
    phih = Uh(1:dof_phi);
    uh1 = Uh(dof_phi+1:dof_phi+dof_u1);
    uh2 = Uh(dof_phi+dof_u1+1:dof_phi+2*dof_u1);
    ph = Uh(dof_phi+2*dof_u1+1:end);
    err_t0 = cputime;
    [phih_L2_error(n), phih_H1_error(n)] = dgL2H1Error(pde.phi,pde.phix,pde.phiy,phih,DarcymeshInfo,sysInfo.Gaussformulas{1},option.basesType_Pp);
    [uh1_L2_error, uh1_H1_error] = dgL2H1Error(pde.u1,pde.u1x,pde.u1y,uh1,StokesmeshInfo,sysInfo.Gaussformulas{1},option.basesType_Fu);
    [uh2_L2_error, uh2_H1_error] = dgL2H1Error(pde.u2,pde.u2x,pde.u2y,uh2,StokesmeshInfo,sysInfo.Gaussformulas{1},option.basesType_Fu);
    [ph_L2_error(n), ~] = dgL2H1Error(pde.p,pde.px,pde.py,ph,StokesmeshInfo,sysInfo.Gaussformulas{1},option.basesType_Fp);
    uh_L2_error(n) = sqrt(uh1_L2_error^2+uh2_L2_error^2);
    uh_H1_error(n) = sqrt(uh1_H1_error^2+uh2_H1_error^2);
    sysInfo.ErrTime = cputime - err_t0;
    disp(['compute Err time: ',num2str(sysInfo.ErrTime)])
    
    disp(['uh_L2_error(',num2str(n),') = ',num2str(uh_L2_error(n))])
    disp(['uh_H1_error(',num2str(n),') = ',num2str(uh_H1_error(n))])
    disp(['ph_L2_error(',num2str(n),') = ',num2str(ph_L2_error(n))])
    disp(['phih_L2_error(',num2str(n),') = ',num2str(phih_L2_error(n))])
    disp(['phih_H1_error(',num2str(n),') = ',num2str(phih_H1_error(n))])
    
    %% other options
    h_S(n) = sum(StokesmeshInfo.hElem)/StokesmeshInfo.Nelems;
    h_D(n) = sum(DarcymeshInfo.hElem)/DarcymeshInfo.Nelems;


    %% creat the structure of save mat file
    n_cycle = ['cycle_',n_str,'_'];
    dgStokesDarcy_output{n,1} = ...
        struct( ...
        [n_cycle,'basestype_Pp'], option.basesType_Pp, ...
        [n_cycle,'basestype_Fu'], option.basesType_Fu, ...
        [n_cycle,'basestype_Fp'], option.basesType_Fp, ...
        [n_cycle,'dof_phi'], dof_phi, ...
        [n_cycle,'dof_u1'], dof_u1, ...
        [n_cycle,'dof_p'], dof_p, ...
        [n_cycle,'phih'], phih, ...
        [n_cycle,'uh1'], uh1, ...
        [n_cycle,'uh2'], uh2, ...
        [n_cycle,'ph'], ph, ...
        [n_cycle,'phih_L2_error'], phih_L2_error(n), ...
        [n_cycle,'phih_H1_error'], phih_H1_error(n), ...
        [n_cycle,'uh_L2_error'], uh_L2_error(n), ...
        [n_cycle,'uh_H1_error'], uh_H1_error(n), ...
        [n_cycle,'ph_L2_error'], ph_L2_error(n), ...
        [n_cycle,'h_S'], h_S(n), ...
        [n_cycle,'h_D'], h_D(n) ...
        );

    
    disp('------------------------------------------------------')
end % for n

% save result to mat file
saveFilename=[setpath_pwd,'/outputmat/',dgfunc_name,date,'.mat'];
save(saveFilename, 'dgStokesDarcy_output');


disp('------------------------------------------------------')
%%
phihL2rate(2:maxIt) = log(phih_L2_error(1:maxIt-1)./phih_L2_error(2:maxIt)) ...
    ./log(h_D(1:maxIt-1)./h_D(2:maxIt));
phihH1rate(2:maxIt) = log(phih_H1_error(1:maxIt-1)./phih_H1_error(2:maxIt)) ...
    ./log(h_D(1:maxIt-1)./h_D(2:maxIt));
uhL2rate(2:maxIt) = log(uh_L2_error(1:maxIt-1)./uh_L2_error(2:maxIt)) ...
    ./log(h_S(1:maxIt-1)./h_S(2:maxIt));
uhH1rate(2:maxIt) =  log(uh_H1_error(1:maxIt-1)./uh_H1_error(2:maxIt)) ...
    ./log(h_S(1:maxIt-1)./h_S(2:maxIt));
phL2rate(2:maxIt) =  log(ph_L2_error(1:maxIt-1)./ph_L2_error(2:maxIt)) ...
    ./log(h_S(1:maxIt-1)./h_S(2:maxIt));

disp('Table: Error')
colname = {'h_S  ', '   ||u-U_h||_0 ','   ||u-U_h||_1','   ||p-p_h||_0'};
disptable(colname, h_S,'%0.2e', uh_L2_error,'%0.5e', ...
    uh_H1_error,'%0.5e', ph_L2_error,'%0.5e');
colname = {'h_D  ', '   ||phi-phi_h||_0 ','   ||phi-phi_h||_1'};
disptable(colname, h_D,'%0.2e', phih_L2_error,'%0.5e', phih_H1_error,'%0.5e');

disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
disp('Table: Error rate')
colname = {'h_S', '   UhL2rate', '   UhH1rate', '   phL2rate'};
disptable(colname, h_S,'%0.2e', uhL2rate,'%0.4f', ...
    uhH1rate,'%0.4f', phL2rate,'%0.4f');
colname = {'h_D', '   phihL2rate', '   phihH1rate'};
disptable(colname, h_D,'%0.2e', phihL2rate,'%0.4f', phihH1rate,'%0.4f');
disp('------------------------------------------------------')

%>>>>>>>>>  close the diary >>>>>>
diary off; % close the diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<
end % function dgPoissonEqn
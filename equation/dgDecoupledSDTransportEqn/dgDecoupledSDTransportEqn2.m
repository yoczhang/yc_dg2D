function [sysErr,sysTime] = dgDecoupledSDTransportEqn2(pde,option,varargin)
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       dgDecoupledSDTransport2, 
%       for the transport eqn, have inflow and outflow boundaryEdges.
%   %---------------------------------------------
%
%   dgSDTransport solve the Stokes-Darcy-Transport equation by DG methods
%
%
%	YcZhang 2/11/2017
%
%   Last modified 2/11/2017
%

%% Default setting of mesh and pde data
if ~exist('option','var') 
    option = []; 
end
if ~exist('pde','var')
    pde = DecoupledSDTransportData2(0,1); % default data
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
ch_L2_error = zeros(maxIt,1); ch_H1_error = zeros(maxIt,1);

phihL2rate = zeros(maxIt,1); phihH1rate = zeros(maxIt,1);
uhL2rate = zeros(maxIt,1); uhH1rate = zeros(maxIt,1);
phL2rate = zeros(maxIt,1);
chL2rate = zeros(maxIt,1); chH1rate = zeros(maxIt,1);
h_S = zeros(maxIt,1); h_D = zeros(maxIt,1); h_c = zeros(maxIt,1);
% Nnodes = zeros(maxIt,1);
% Ndof = zeros(maxIt,1);
dgDecoupledSDTransport2_output = cell(maxIt,1);

%% set path and file name
load('setpath_pwd.mat')
dgfunc_name = 'dgDecoupledSDTransport2_';

%>>>>>>>>>  creat log file >>>>>>>>
date = datestr(now,31); 
    %> capture the now time, 31 stands for the scheme: 2000-03-01 15:45:17
logFilename=[setpath_pwd,'/logs/',dgfunc_name,date,'_log.txt'];
diary(logFilename);
diary on; % begin diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<


disp('********************  dg Stokes-Darcy Transport  ********************')
disp('------------------------------------------------------')
disp('Stokes-Darcy Transport Equation:')
disp(['   K = ',num2str(pde.K)])
disp(['   phi = ',func2str(pde.phi)])
disp(['   u1 = ',func2str(pde.u1)])
disp(['   u2 = ',func2str(pde.u2)])
disp(['   c = ',func2str(pde.c)])
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
disp('Stokes eqn fluid kinematic viscosity:')
disp(['   mu = ', num2str(pde.mu)])
disp('Darcy eqn PermeabilityCoefficient:')
disp(['   PermeabiltiyCoeff (i.e. pde.K) = ', num2str(pde.K)])
disp('Transport eqn DiffusivityCoefficient:')
disp(['   DiffusivityCoeff = ', num2str(pde.DiffusivityCoeff)])
disp(['startingTime: ',num2str(option.startingT)])
disp(['terminalTime: ',num2str(option.terminalT)])
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
disp('DG parameters: ')
disp(['   basesType_phi: ', option.basesType_Pp])
disp(['   basesType_u: ', option.basesType_Fu])
disp(['   basesType_p: ', option.basesType_Fp])
disp(['   basesType_c: ', option.basesType_c])
disp(['   penalty pars: p_epsilon: ', num2str(option.p_epsilon)])
disp(['   penalty pars: p_sigma: ', num2str(option.p_sigma)])
disp(['   penalty pars: p_beta: ', num2str(option.p_beta)])
disp(['   penalty pars: p_epsilon_c: ', num2str(option.p_epsilon_c)])
disp(['   penalty pars: p_sigma_c: ', num2str(option.p_sigma_c)])
disp(['   penalty pars: p_beta_c: ', num2str(option.p_beta_c)])
disp(['   maxIt: ', num2str(maxIt)])

%% DG methods
disp('------------------------------------------------------')
disp('Solving dgDecoupledSDTransport2: ')
disp('------------------------------------------------------')
% format long e
for n = 1:maxIt
    close all
    n_str = num2str(n);
    disp('')
    n_th_cycle = ['the ',n_str,'-th cycle(i.e. the ', n_str, '-th mesh)'];
    disp(n_th_cycle)
    
    %---- mesh 1--- S:[0,1]x[0.5,1]. D:[0,1]x[0,0.5].
	%-------------------- Quad mesh ---------------------
% 	% get Stokes mesh information
% 	Smesh_name = ['S_quadmesh_',n_str_x,'times',n_str_y,'_[0_1]_[1-2_1]'];
% 	load(Smesh_name);
%     Snode = vertices'; Selem = elements';
%         
% 	% get Darcy mesh information
% 	Dmesh_name = ['D_quadmesh_',n_str_x,'times',n_str_y,'_[0_1]_[0_1-2]'];
% 	load(Dmesh_name);
%     Dnode = vertices'; Delem = elements';
	%-----------------------------------------------------
    
    %---- mesh 2--- S:[0,1]x[1,2]. D:[0,1]x[0,1].
	%-------------------- Quad mesh ---------------------
	% get Stokes mesh information
    h_x = 1/2^(n+1); h_y = h_x;
	h_partition=[h_x,h_y];
	left=0;right=1;bottom=0;top=1;
	[M,T]=generate_quad_P_T(left,right,bottom,top,h_partition,1);
	Snode = M'; Selem = T';
    disp(['Stokes domain: [',num2str(left),',',num2str(right),']','x','[',num2str(bottom),',',num2str(top),'].'])
        
	% get Darcy mesh information
	left=0;right=1;bottom=-1;top=0;
	[M,T]=generate_quad_P_T(left,right,bottom,top,h_partition,1);
	Dnode = M'; Delem = T';
    disp(['Darcy domain: [',num2str(left),',',num2str(right),']','x','[',num2str(bottom),',',num2str(top),'].'])
	%-----------------------------------------------------
    
    
    %------------- Stokes meshInfo ----------------------
    StokesmeshInfo = polyMeshAuxStructure(Snode, Selem);
    DiffusivityCoeffs = ones(size(StokesmeshInfo.elem,1),1);
    StokesmeshInfo.DiffusivityCoeffs = pde.DiffusivityCoeff*DiffusivityCoeffs;
    %----------------------------------------------------------
    
    %-------------- Darcy meshInfo ----------------------
    DarcymeshInfo = polyMeshAuxStructure(Dnode, Delem);
    %--- get the enriched elems in the subfunction:
	enrichedElem = getEnrichedElemForDarcyDomain(DarcymeshInfo.node, DarcymeshInfo.elem, 3);
    DiffusivityCoeffs = ones(size(DarcymeshInfo.elem,1),1);
	DiffusivityCoeffs(enrichedElem) = pde.DiffusivityCoeff;
	DarcymeshInfo.DiffusivityCoeffs = DiffusivityCoeffs;
    % the Darcy permeability matrix coefficient.
    PermeabilityCoeffs = pde.K*ones(size(DarcymeshInfo.elem,1),1);
    DarcymeshInfo.PermeabilityCoeffs = PermeabilityCoeffs;
    %----------------------------------------------------------
    
    %interfacemeshInfo=test_coupled_interface_1(StokesmeshInfo, DarcymeshInfo, 'xcoord');
    %plotPolyMsh(StokesmeshInfo)
    %patchPlotMesh(StokesmeshInfo.node, StokesmeshInfo.elem)
    %plotPolyMsh(DarcymeshInfo)
    %patchPlotMesh(DarcymeshInfo.node, DarcymeshInfo.elem)

    %% solve equations
    solve_t0 = cputime;
    [SD_Uh, Ch, sysInfo]= dgDecoupledSDTransportSolve2(StokesmeshInfo,DarcymeshInfo,pde,option);
    sysInfo.SoverTime = cputime - solve_t0;
    disp(['the solve-func time: ',num2str(sysInfo.SoverTime)])
    
    %% compute the err
    dof_phi = sysInfo.dof_phi;
    dof_u1 = sysInfo.dof_u1;
    dof_p = sysInfo.dof_p;
    dof_SD = dof_phi + 2*dof_u1 + dof_p;
    dof_Sc = sysInfo.dof_Sc;
    dof_Dc = sysInfo.dof_Dc;
    phih = SD_Uh(1:dof_phi);
    uh1 = SD_Uh(dof_phi+1:dof_phi+dof_u1);
    uh2 = SD_Uh(dof_phi+dof_u1+1:dof_phi+2*dof_u1);
    ph = SD_Uh(dof_phi+2*dof_u1+1:dof_SD);
    ch_D = Ch(1:dof_Dc);
    ch_S = Ch(dof_Dc+1:dof_Dc+dof_Sc);
    
    err_t0 = cputime;
    [phih_L2_error(n), phih_H1_error(n)] = dgL2H1Error(pde.phi,pde.phix,pde.phiy,phih,DarcymeshInfo,sysInfo.Gaussformulas{1},option.basesType_Pp);
    [uh1_L2_error, uh1_H1_error] = dgL2H1Error(pde.u1,pde.u1x,pde.u1y,uh1,StokesmeshInfo,sysInfo.Gaussformulas{1},option.basesType_Fu);
    [uh2_L2_error, uh2_H1_error] = dgL2H1Error(pde.u2,pde.u2x,pde.u2y,uh2,StokesmeshInfo,sysInfo.Gaussformulas{1},option.basesType_Fu);
    [ph_L2_error(n), ~] = dgL2H1Error(pde.p,pde.px,pde.py,ph,StokesmeshInfo,sysInfo.Gaussformulas{1},option.basesType_Fp);
    [chD_L2_error, chD_H1_error] = dgTimeL2H1Error(sysInfo.terminalT,pde.c,pde.cx,pde.cy,ch_D,DarcymeshInfo,sysInfo.Gaussformulas{1},option.basesType_c);
    [chS_L2_error, chS_H1_error] = dgTimeL2H1Error(sysInfo.terminalT,pde.c,pde.cx,pde.cy,ch_S,StokesmeshInfo,sysInfo.Gaussformulas{1},option.basesType_c);
    uh_L2_error(n) = sqrt(uh1_L2_error^2+uh2_L2_error^2);
    uh_H1_error(n) = sqrt(uh1_H1_error^2+uh2_H1_error^2);
    ch_L2_error(n) = sqrt(chD_L2_error^2 + chS_L2_error^2);
    ch_H1_error(n) = sqrt(chD_H1_error^2 + chS_H1_error^2);
    sysInfo.ErrTime = cputime - err_t0;
    disp(['compute Err time: ',num2str(sysInfo.ErrTime)])
    
    disp(['uh_L2_error(',num2str(n),') = ',num2str(uh_L2_error(n))])
    disp(['uh_H1_error(',num2str(n),') = ',num2str(uh_H1_error(n))])
    disp(['ph_L2_error(',num2str(n),') = ',num2str(ph_L2_error(n))])
    disp(['phih_L2_error(',num2str(n),') = ',num2str(phih_L2_error(n))])
    disp(['phih_H1_error(',num2str(n),') = ',num2str(phih_H1_error(n))])
    disp(['ch_L2_error(',num2str(n),') = ',num2str(ch_L2_error(n))])
    disp(['ch_H1_error(',num2str(n),') = ',num2str(ch_H1_error(n))])
    
    %% other options
    h_S(n) = sum(StokesmeshInfo.hElem)/StokesmeshInfo.Nelems;
    h_D(n) = sum(DarcymeshInfo.hElem)/DarcymeshInfo.Nelems;
    h_c(n) = 1./(sqrt(StokesmeshInfo.Nnodes + DarcymeshInfo.Nnodes)-1);


    %% creat the structure of save mat file
    n_cycle = ['cycle_',n_str,'_'];
    dgDecoupledSDTransport2_output{n,1} = ...
        struct( ...
        [n_cycle,'StokesmeshInfo'], StokesmeshInfo, ...
        [n_cycle,'DarcymeshInfo'], DarcymeshInfo, ...
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
        [n_cycle,'ch_S'], ch_S, ...
        [n_cycle,'ch_D'], ch_D, ...
        [n_cycle,'phih_L2_error'], phih_L2_error(n), ...
        [n_cycle,'phih_H1_error'], phih_H1_error(n), ...
        [n_cycle,'uh_L2_error'], uh_L2_error(n), ...
        [n_cycle,'uh_H1_error'], uh_H1_error(n), ...
        [n_cycle,'ph_L2_error'], ph_L2_error(n), ...
        [n_cycle,'ch_L2_error'], ch_L2_error(n), ...
        [n_cycle,'ch_H1_error'], ch_H1_error(n), ...
        [n_cycle,'h_S'], h_S(n), ...
        [n_cycle,'h_D'], h_D(n), ...
        [n_cycle,'h_c'], h_c(n) ...
        );

    
    disp('------------------------------------------------------')
end % for n

% save result to mat file
saveFilename=[setpath_pwd,'/outputmat/',dgfunc_name,date,'.mat'];
save(saveFilename, 'dgDecoupledSDTransport2_output');


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
chL2rate(2:maxIt) =  log(ch_L2_error(1:maxIt-1)./ch_L2_error(2:maxIt)) ...
    ./log(h_c(1:maxIt-1)./h_c(2:maxIt));
chH1rate(2:maxIt) =  log(ch_H1_error(1:maxIt-1)./ch_H1_error(2:maxIt)) ...
    ./log(h_c(1:maxIt-1)./h_c(2:maxIt));

%---- error
disp('Table: Error')
colname = {'h_S  ', '   ||u-U_h||_0 ','   ||u-U_h||_1','   ||p-p_h||_0'};
disptable(colname, h_S,'%0.2e', uh_L2_error,'%0.5e', ...
    uh_H1_error,'%0.5e', ph_L2_error,'%0.5e');
colname = {'h_D  ', '   ||phi-phi_h||_0 ','   ||phi-phi_h||_1'};
disptable(colname, h_D,'%0.2e', phih_L2_error,'%0.5e', phih_H1_error,'%0.5e');
colname = {'h_c  ','   ||c-c_h||_0','   ||c-c_h||_1'};
disptable(colname, h_c,'%0.2e', ch_L2_error,'%0.5e', ch_H1_error,'%0.5e');

disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')

%---- rate
disp('Table: Error rate')
colname = {'h_S', '   UhL2rate', '   UhH1rate', '   phL2rate'};
disptable(colname, h_S,'%0.2e', uhL2rate,'%0.4f', ...
    uhH1rate,'%0.4f', phL2rate,'%0.4f');
colname = {'h_D', '   phihL2rate', '   phihH1rate'};
disptable(colname, h_D,'%0.2e', phihL2rate,'%0.4f', phihH1rate,'%0.4f');
colname = {'h_c', '   chL2rate', '   chH1rate'};
disptable(colname, h_c,'%0.2e', chL2rate,'%0.4f', chH1rate,'%0.4f');
disp('------------------------------------------------------')

%>>>>>>>>>  close the diary >>>>>>
diary off; % close the diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<
end % function dgPoissonEqn
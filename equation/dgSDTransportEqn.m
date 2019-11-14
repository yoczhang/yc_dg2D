function [sysErr,sysTime] = dgSDTransportEqn(pde,option,varargin)
%
%   dgSDTransport solve the Stokes-Darcy-Transport equation by DG methods
%
%
%	YcZhang 10/9/2017
%
%   Last modified 10/9/2017
%

%% Default setting of mesh and pde data
if ~exist('option','var') 
    option = []; 
end
if ~exist('pde','var')
    pde = SDTransportData(0,2); % default data
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
dgSDTransport_output = cell(maxIt,1);

%% set path and file name
load('setpath_pwd.mat')
dgfunc_name = 'dgSDTransport_';

%>>>>>>>>>  creat log file >>>>>>>>
date = datestr(now,31); 
    %> capture the now time, 31 stands for the scheme: 2000-03-01 15:45:17
logFilename=[setpath_pwd,'/logs/',dgfunc_name,date,'_log.txt'];
diary(logFilename);
diary on; % begin diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<


disp('********************  dg Stokes-Darcy Transport  ********************')
disp('------------------------------------------------------')
disp('poisson Equation:')
disp(['   k = ',pde.k])
disp(['   phi = ',func2str(pde.phi)])
disp(['   u1 = ',func2str(pde.u1)])
disp(['   u2 = ',func2str(pde.u2)])
disp('DG parameters: ')
disp(['   basesType_phi: ', option.basesType_Pp])
disp(['   basesType_u: ', option.basesType_Fu])
disp(['   basesType_p: ', option.basesType_Fp])
disp(['   basesType_c: ', option.basesType_c])
disp(['   penalty pars: p_epsilon: ', num2str(option.p_epsilon)])
disp(['   penalty pars: p_sigma: ', num2str(option.p_sigma)])
disp(['   penalty pars: p_beta: ', num2str(option.p_beta)])
disp(['   maxIt: ', num2str(maxIt)])

%% DG methods
disp('------------------------------------------------------')
disp('Solving dgSDTransport: ')
disp('------------------------------------------------------')
% format long e
for n = 1:maxIt
    close all
    n_str = num2str(n);
    n_str_x = num2str(2^(n+1));
    n_str_y = num2str(2^n);
    
    disp('')
    n_th_cycle = ['the ',n_str,'-th cycle(i.e. the ', n_str, '-th mesh)'];
    disp(n_th_cycle)
    
    %% get Stokes mesh information
    % S_quadmesh_4times2_[0_1]_[1-2_1]
    mesh_name = ['S_quadmesh_',n_str_x,'times',n_str_y,'_[0_1]_[1-2_1]'];
    load(mesh_name);
    StokesmeshInfo = polyMeshAuxStructure(vertices, elements);
    plotPolyMsh(StokesmeshInfo)
    %patchPlotMesh(StokesmeshInfo.node, StokesmeshInfo.elem)
    
    DiffusivityCoeffs = ones(size(StokesmeshInfo.elem,1),1);
    StokesmeshInfo.DiffusivityCoeffs = DiffusivityCoeffs;
    
    
    %% get Darcy mesh information
    mesh_name = ['D_quadmesh_',n_str_x,'times',n_str_y,'_[0_1]_[0_1-2]'];
    load(mesh_name);
    DarcymeshInfo = polyMeshAuxStructure(vertices, elements);
    plotPolyMsh(DarcymeshInfo)
    %patchPlotMesh(DarcymeshInfo.node, DarcymeshInfo.elem)
    
    %interfacemeshInfo=test_coupled_interface_1(StokesmeshInfo, DarcymeshInfo, 'xcoord');
    
    %--- get the enriched elems in the subfunction:
	enrichedElem = getEnrichedElemForDarcyDomain(DarcymeshInfo.node, DarcymeshInfo.elem, 1);
    
    DiffusivityCoeffs = ones(size(DarcymeshInfo.elem,1),1);
	DiffusivityCoeffs(enrichedElem) = 1;
	DarcymeshInfo.DiffusivityCoeffs = DiffusivityCoeffs;

    %% solve equations
    solve_t0 = cputime;
    [SD_Uh, Ch, sysInfo]= dgDecoupledSDTransportSolve(StokesmeshInfo,DarcymeshInfo,pde,option);
    %[SD_Uh, Ch, sysInfo]= dgCoupledSDTransportSolve(StokesmeshInfo,DarcymeshInfo,pde,option);
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
    
    
    %% other options
    h_S(n) = sum(StokesmeshInfo.hElem)/StokesmeshInfo.Nelems;
    h_D(n) = sum(DarcymeshInfo.hElem)/DarcymeshInfo.Nelems;
    h_c(n) = 1./(sqrt(StokesmeshInfo.Nnodes + DarcymeshInfo.Nnodes)-1);


    %% creat the structure of save mat file
    n_cycle = ['cycle_',n_str,'_'];
    dgSDTransport_output{n,1} = ...
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
        [n_cycle,'ch_L2_error'], ch_L2_error(n), ...
        [n_cycle,'ch_H1_error'], ch_H1_error(n), ...
        [n_cycle,'h_S'], h_S(n), ...
        [n_cycle,'h_D'], h_D(n) ...
        );

    
    disp('------------------------------------------------------')
end % for n

% save result to mat file
saveFilename=[setpath_pwd,'/outputmat/',dgfunc_name,date,'.mat'];
save(saveFilename, 'dgSDTransport_output');


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
    ./log(h_S(1:maxIt-1)./h_S(2:maxIt));
chH1rate(2:maxIt) =  log(ch_H1_error(1:maxIt-1)./ch_H1_error(2:maxIt)) ...
    ./log(h_S(1:maxIt-1)./h_S(2:maxIt));

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
function [sysErr,sysTime] = dgTransportEqn(pde,option,varargin)
%
%   dgParabolicEqn solve the parabolic equation by DG methods
% 
%
%	YcZhang 5/9/2017
%
%   Last modified 5/9/2017
%

%% Default setting of mesh and pde data
if ~exist('option','var') 
    option = []; 
end
if ~exist('pde','var')
    pde = transportData(0,2); % default data
end


%% Parameters
option = dgOption(option);
maxIt = option.maxIt;


%% Initialize err
sysErr = zeros(maxIt,1); % not used
sysTime = zeros(maxIt,1); % not used

uh_L2_error = zeros(maxIt,1); uh_H1_error = zeros(maxIt,1);
uhL2rate = zeros(maxIt,1); uhH1rate = zeros(maxIt,1);
h = zeros(maxIt,1); 
dgTransport_output = cell(maxIt,1);


%% set path and file name
load('setpath_pwd.mat')
dgfunc_name = 'dgTransport_';

%>>>>>>>>>  creat log file >>>>>>>>
date = datestr(now,31); 
    %> capture the now time, 31 stands for the scheme: 2000-03-01 15:45:17
logFilename=[setpath_pwd,'/logs/',dgfunc_name,date,'_log.txt'];
diary(logFilename);
diary on; % begin diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<

disp('********************  dg Transport  ********************')
disp('------------------------------------------------------')
disp('Transport Equation:')
disp(['   k_11 = ',func2str( pde.k_11)])
disp(['   k_22 = ',func2str( pde.k_22)])
disp(['   u = ',func2str(pde.u)])
disp('DG parameters: ')
disp(['   basesType_trial: ', option.basesType_trial])
disp(['   basesType_test: ', option.basesType_test])
disp(['   penalty pars: p_epsilon: ', num2str(option.p_epsilon)])
disp(['   penalty pars: p_sigma: ', num2str(option.p_sigma)])
disp(['   penalty pars: p_beta: ', num2str(option.p_beta)])
disp(['   maxIt: ', num2str(maxIt)])
disp('Time parameters: ')
disp(['   starting time: ', num2str(option.startingT)])
disp(['   terminal time: ', num2str(option.terminalT)])
disp(['   time sheme: theta:',num2str(option.theta)]);

%% DG methods
disp('------------------------------------------------------')
disp('Solving dgTransport: ')
disp('------------------------------------------------------')
format long e

if strcmp(option.verifyCovergence,'space')
    for n = 1:maxIt
%         n_str = num2str(n);
%         mesh_name = ['Polygon_',n_str];
%         load(mesh_name);
        n_str = num2str(4^n);
        mesh_name = ['quadmesh_',n_str,'elem'];
        load(mesh_name);

        disp('')
        n_th_cycle = [option.verifyCovergence,': the ',n_str,'-th cycle(i.e. the ', n_str, '-th mesh)'];
        disp(n_th_cycle)

        %% get mesh information
        node = vertices;
        elem = elements;
        meshInfo = polyMeshAuxStructure(node, elem);
        %plotPolyMsh(meshInfo)

        %% solve equations
        solve_t0 = cputime;
        [Uh, sysInfo]= dgTransportSolve(meshInfo,pde,option);
        sysInfo.SoverTime = cputime - solve_t0;
        disp(['the solve-func time: ',num2str(sysInfo.SoverTime)])

        %% compute the err
        err_t0 = cputime;
        [uh_L2_error(n), uh_H1_error(n)] = dgTimeL2H1Error(sysInfo.terminalT,pde.u,pde.ux,pde.uy,Uh,meshInfo,sysInfo.Gaussformulas{1},option.basesType_trial);
        sysInfo.ErrTime = cputime - err_t0;
        disp(['compute Err time: ',num2str(sysInfo.ErrTime)])

        %% other options
        h(n) = sum(meshInfo.hElem)/meshInfo.Nelems;

        %% creat the structure of save mat file
        n_cycle = ['cycle_',n_str,'_'];
        dgTransport_output{n,1} = ...
            struct( ...
            [n_cycle,'verifyCovergence'], option.verifyCovergence, ...
            [n_cycle,'basestype_trial'], option.basesType_trial, ...
            [n_cycle,'basestype_test'], option.basesType_test, ...
            [n_cycle,'dof_u'], sysInfo.dof_u, ...
            [n_cycle,'uh'], Uh, ...
            [n_cycle,'uh_L2_error'], uh_L2_error(n), ...
            [n_cycle,'uh_H1_error'], uh_H1_error(n), ...
            [n_cycle,'h'], h(n) ...
            );

        disp('------------------------------------------------------')
    end % for n
    disp('------------------------------------------------------')
    disp('for space error and convergence:')
end % if

if strcmp(option.verifyCovergence,'time')
    Dt = option.Dt;
    for n = 1:maxIt
        option.Dt = Dt/2^(n-1);
        n_str = num2str(n);
        mesh_name = 'quadmesh_16384elem';
        load(mesh_name);

        disp('')
        n_th_cycle = [option.verifyCovergence,': the ',n_str,'-th time cycle'];
        disp(n_th_cycle)

        %% get mesh information
        node = vertices;
        elem = elements;
        meshInfo = polyMeshAuxStructure(node, elem);
        %plotPolyMsh(meshInfo)

        %% solve equations
        solve_t0 = cputime;
        [Uh, sysInfo]= dgTransportSolve(meshInfo,pde,option);
        sysInfo.SoverTime = cputime - solve_t0;
        disp(['the solve-func time: ',num2str(sysInfo.SoverTime)])

        %% compute the err
        err_t0 = cputime;
        [uh_L2_error(n), uh_H1_error(n)] = dgTimeL2H1Error(sysInfo.terminalT,pde.u,pde.ux,pde.uy,Uh,meshInfo,sysInfo.Gaussformulas{1},option.basesType_trial);
        sysInfo.ErrTime = cputime - err_t0;
        disp(['compute Err time: ',num2str(sysInfo.ErrTime)])

        %% other options
        h(n) = option.Dt;

        %% creat the structure of save mat file
        n_cycle = ['cycle_',n_str,'_'];
        dgTransport_output{n,1} = ...
            struct( ...
            [n_cycle,'verifyCovergence'], option.verifyCovergence, ...
            [n_cycle,'basestype_trial'], option.basesType_trial, ...
            [n_cycle,'basestype_test'], option.basesType_test, ...
            [n_cycle,'dof_u'], sysInfo.dof_u, ...
            [n_cycle,'uh'], Uh, ...
            [n_cycle,'uh_L2_error'], uh_L2_error(n), ...
            [n_cycle,'uh_H1_error'], uh_H1_error(n), ...
            [n_cycle,'h'], h(n) ...
            );

        disp('------------------------------------------------------')
    end % for n
    disp('------------------------------------------------------')
    disp('for time error and convergence:')
end % 



%% save result to mat file
saveFilename=[setpath_pwd,'/outputmat/',dgfunc_name,date,'.mat'];
save(saveFilename, 'dgTransport_output');


%%
uhL2rate(2:maxIt) = log(uh_L2_error(1:maxIt-1)./uh_L2_error(2:maxIt)) ...
    ./log(h(1:maxIt-1)./h(2:maxIt));
uhH1rate(2:maxIt) =  log(uh_H1_error(1:maxIt-1)./uh_H1_error(2:maxIt)) ...
    ./log(h(1:maxIt-1)./h(2:maxIt));

disp('Table: Error')
colname = {'h  ', '   ||u-U_h||_0 ','   ||u-U_h||_1'};
disptable(colname, h,'%0.2e', uh_L2_error,'%0.5e', ...
    uh_H1_error,'%0.5e');

disp('Table: Error rate')
colname = {'h', '   UhL2rate', '   UhH1rate'};
disptable(colname, h,'%0.2e', uhL2rate,'%0.4f', ...
    uhH1rate,'%0.4f');
disp('------------------------------------------------------')

%>>>>>>>>>  close the diary >>>>>>
diary off; % close the diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<
end % function dgPoissonEqn
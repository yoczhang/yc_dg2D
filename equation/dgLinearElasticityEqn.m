function [sysErr,sysTime] = dgLinearElasticityEqn(pde,option,varargin)
%
%   dgLinearElasticityEqn solve the linear Elasticity equation by DG methods
%
%
%	YcZhang 22/9/2017
%
%   Last modified 22/9/2017
%

%% Default setting of mesh and pde data
if ~exist('option','var') 
    option = []; 
end
if ~exist('pde','var')
    pde = linearElasticityData(1,1); % default data
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
% Nnodes = zeros(maxIt,1);
% Ndof = zeros(maxIt,1);
dgLinearElasticity_output = cell(maxIt,1);

%% set path and file name
load('setpath_pwd.mat')
dgfunc_name = 'dgLinearElasticity_';

%>>>>>>>>>  creat log file >>>>>>>>
date = datestr(now,31); 
    %> capture the now time, 31 stands for the scheme: 2000-03-01 15:45:17
logFilename=[setpath_pwd,'/logs/',dgfunc_name,date,'_log.txt'];
diary(logFilename);
diary on; % begin diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<


disp('********************  dg Linear Elasticity  ********************')
disp('------------------------------------------------------')
disp('LinearElasticity Equation:')
disp(['   u1 = ',func2str(pde.u1)])
disp(['   u2 = ',func2str(pde.u2)])
disp('DG parameters: ')
disp(['   basesType_u: ', option.basesType_u])
disp(['   penalty pars: p_epsilon: ', num2str(option.p_epsilon)])
disp(['   penalty pars: p_sigma: ', num2str(option.p_sigma)])
disp(['   penalty pars: p_beta: ', num2str(option.p_beta)])
disp(['   maxIt: ', num2str(maxIt)])

%% DG methods
disp('------------------------------------------------------')
disp('Solving dgLinearElasticity: ')
disp('------------------------------------------------------')
% format long e
for n = 1:maxIt
    n_str = num2str(n);
    %mesh_name = ['Polygon_',n_str];
    %load(mesh_name);
    
    disp('')
    n_th_cycle = ['the ',n_str,'-th cycle(i.e. the ', n_str, '-th mesh)'];
    disp(n_th_cycle)
    
    %% get mesh information

    left=0; right=1; bottom=0; top=1;
	h_x = 1/2^(n+1); h_y = h_x;
	h_partition=[h_x,h_y];
	[node, elem] = generate_Tri_P_T(left,right,bottom,top,h_partition);
    
    meshInfo = polyMeshAuxStructure(node, elem);
    %plotPolyMsh(meshInfo)

    %% solve equations
    solve_t0 = cputime;
    [Uh, sysInfo]= dgLinearElasticitySolve(meshInfo,pde,option);
    sysInfo.SoverTime = cputime - solve_t0;
    disp(['the solve-func time: ',num2str(sysInfo.SoverTime)])
    
    %% compute the err
    dof_u1 = sysInfo.dof_u1;
    uh1 = Uh(1:dof_u1);
    uh2 = Uh(1+dof_u1:2*dof_u1);
    err_t0 = cputime;
    [uh1_L2_error, uh1_H1_error] = dgL2H1Error(pde.u1,pde.u1x,pde.u1y,uh1,meshInfo,sysInfo.Gaussformulas{1},option.basesType_u);
    [uh2_L2_error, uh2_H1_error] = dgL2H1Error(pde.u2,pde.u2x,pde.u2y,uh2,meshInfo,sysInfo.Gaussformulas{1},option.basesType_u);
    uh_L2_error(n) = sqrt(uh1_L2_error^2+uh2_L2_error^2);
    uh_H1_error(n) = sqrt(uh1_H1_error^2+uh2_H1_error^2);
    sysInfo.ErrTime = cputime - err_t0;
    disp(['compute Err time: ',num2str(sysInfo.ErrTime)])
    
    
    %% other options
    h(n) = 1./(sqrt(meshInfo.Nnodes)-1);


    %% creat the structure of save mat file
    n_cycle = ['cycle_',n_str,'_'];
    dgLinearElasticity_output{n,1} = ...
        struct( ...
        [n_cycle,'basestype_u'], option.basesType_u, ...
        [n_cycle,'dof_u1'], dof_u1, ...
        [n_cycle,'uh1'], uh1, ...
        [n_cycle,'uh2'], uh2, ...
        [n_cycle,'uh_L2_error'], uh_L2_error(n), ...
        [n_cycle,'uh_H1_error'], uh_H1_error(n), ...
        [n_cycle,'h'], h(n) ...
        );

    
    disp('------------------------------------------------------')
end % for n

% save result to mat file
saveFilename=[setpath_pwd,'/outputmat/',dgfunc_name,date,'.mat'];
save(saveFilename, 'dgLinearElasticity_output');


disp('------------------------------------------------------')
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
function [sysErr,sysTime] = hcDGPoissonEqn(pde,option,varargin)
%
%   dgPoissonEqn solve the poisson equation by DG methods
% 
%
%	YcZhang  Jul.08.2019
%
%   Last modified  Jul.08.2019
%

%% Default setting of mesh and pde data
if ~exist('option','var') 
    option = []; 
end
if ~exist('pde','var')
    pde = hcDGPoissonData(0,5); % default data
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


%% set path and file name

%>>>>>>>>>  creat log file >>>>>>>>
% load('setpath_pwd.mat',setpath_pwd)
% dgfunc_name = 'hcDGPoisson_';
% 
% date = datestr(now,31); 
%     %> capture the now time, 31 stands for the scheme: 2000-03-01 15:45:17
% logFilename=[setpath_pwd,'/logs/',dgfunc_name,date,'_log.txt'];
% diary(logFilename);
% diary on; % begin diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<

disp('********************  hcDG Poisson  ********************')
disp('------------------------------------------------------')
disp('poisson Equation:')
disp(['   k_11 = ',func2str( pde.k_11)])
disp(['   k_22 = ',func2str( pde.k_22)])
disp(['   u = ',func2str(pde.u)])
disp('DG parameters: ')
disp(['   maxIt: ', num2str(maxIt)])

%% DG methods
disp('------------------------------------------------------')
disp('Solving hcDGPoisson: ')
disp('------------------------------------------------------')
format long e
for n = 1:maxIt
    %-------------------- Poly mesh ---------------------
%     mesh_name = ['Polygon_',n_str];
%     load(mesh_name);
%     
%     disp('')
%     n_th_cycle = ['the ',n_str,'-th cycle(i.e. the ', n_str, '-th mesh)'];
%     disp(n_th_cycle)
%     node = vertices;
%     elem = elements;
    %------------------------------------------------------
    
    %-------------------- Tri mesh ---------------------
    left=0; right=1; bottom=0; top=1;
    h_x = 1/2^(n+1);
    h_y = h_x;
    h_partition=[h_x,h_y];
    [node, elem] = generate_Tri_P_T(left,right,bottom,top,h_partition);
    %-----------------------------------------------------
        
    
    %% get mesh information
    meshInfo = polyMeshAuxStructure(node, elem);
%     plotPolyMsh(meshInfo)

    %% solve equations
    solve_t0 = cputime;
    [Uh, sysInfo]= hcDGPoissonSolve(meshInfo,pde,option);
    sysInfo.SoverTime = cputime - solve_t0;
    disp(['the solve-func time: ',num2str(sysInfo.SoverTime)])
    
    %% compute the err
    err_t0 = cputime;
    [uh_L2_error(n), uh_H1_error(n)] = dgL2H1Error(pde.u,pde.ux,pde.uy,Uh,meshInfo,sysInfo.Gaussformulas{1},option.basesType);
    sysInfo.ErrTime = cputime - err_t0;
    disp(['compute Err time: ',num2str(sysInfo.ErrTime)])
    
    %% other options
    h(n) = 1./(sqrt(meshInfo.Nnodes)-1);

    disp('------------------------------------------------------')
end % for n

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
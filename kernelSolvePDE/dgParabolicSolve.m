function [Uh, sysInfo, system] = dgParabolicSolve(meshInfo,pde,option)
%
%   We let  Npoints denote the number of Gauss-Points on T,
%               NTbases denote the number of LOCAL bases on each element.
%
%   input:
%       meshInfo, mesh info structure.
%       pde, pde info structure.
%       option, some options about the euqation.
%
%   output:
%       Uh, [Nbases*Nelems x 1], matrix with uh.
%       system, {stiff-matrix,rhs-term}.
%       sysInfo.assembleLocalsolverTime,
%       sysInfo.solveGlobalUhatTime, 
%       sysInfo.solveLocalsolverTime,
%       sysInfo.hdgDiffusionSoverTime.
%
%
%
%   YcZhang 3/9/2017
%
%   Last modified 4/9/2017
%


if ~exist('option','var')
    option = dgOption(option);
end

trial_k = basesType2degreek(option.basesType_trial);
test_k = basesType2degreek(option.basesType_test);
dof_u = meshInfo.Nelems * (trial_k+2)*(trial_k+1)/2;


%%  generate the Gaussformulas
% 2d 
[Q1,Q2,W] = quadRule2D(2*trial_k+5); % Gauss-Points on Ref elem
formulaGauss2D = [Q1' Q2' W']; 
    %> In order to facilitate the calculation, we structure the 'formula' matrix,
    %> [Npoints x 3],
    %> the first column is the x-coordinates of all Gauss-Points,
    %> the second column is the y-coordinates of all Gauss-Points,
    %> the third is the weights of all Gauss-Points.

% 1d
% [q,w] =quadRule1D(2*trial_k+1); % Gauss-Points on [0,1]
[q,w] =quadRule1D(12); % Gauss-Points on [0,1]
formulaGauss1D = [q' w'];
    %> [Npoints x 2],
    %> the first column is the 1D coordinates of all Gauss-Points,
    %> the second is the weights of all Gauss-Points.

% get the uniform Gaussformulas, is a CELL structure data. 
Gaussformulas = {formulaGauss2D,... % 2d quad formula for ERRPRs and VAR COEFFs. 
    formulaGauss2D, ... % 2d quad formula for CONSTANT COEFFs.
    formulaGauss1D, ... % 1d quad formula.
    formulaGauss1D}; % 1d quad formula for ERRORs.


%% edge information
% we need to get the interior edges 
% and DEFINE the B.C. edges, here we set all the bd edges as the Dir edges.
interEdgeIndex = meshInfo.interEdgeIndex; % [Ninter x 1]
DirichletEdgeIndex = meshInfo.bdEdgeIndex; % [Ndir x 1]

if ~isfield(meshInfo, 'interEdgeIndex')
    meshInfo.interEdgeIndex = interEdgeIndex;
end

if ~isfield(meshInfo, 'DirichletEdgeIndex')
    meshInfo.DirichletEdgeIndex = DirichletEdgeIndex;
end

%% Dt setting
h = sum(meshInfo.hElem)/meshInfo.Nelems;
if option.TimeStepAccordingSpaceStep
    if option.theta == 1/2
        Dt = h^((trial_k+1)/2);
    else
        Dt = h^(trial_k+1);
    end
    if Dt >= (option.terminalT-option.startingT) % here we need the (NT = (...)/Dt) >=2.
        Dt = (option.terminalT-option.startingT)/2;
    end
else
    Dt = 1/10000;
end % if 1

NT = floor((option.terminalT-option.startingT)/Dt);
if NT <= 1
    NT = 2;
end
disp(['Dt = ',num2str(Dt)])
disp(['NT = ',num2str(NT)])

%% assemeble matrix and rhs
tic; %<<<<<<<<<<<<<< tic2
Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.k_11;
Coeffs_func{2} = pde.k_22;
Coeffs_func{3} = pde.funcZero;

%> here the notations according the 'pure parabolic problem' from Beatrice Riviere's DG book.
% the init matrix and rhs(t==0)
[M, A_uv, rhs0] = ...
    matElemCoeffsDtrialDtestParabolic(Coeffs_func,pde.u0,meshInfo,Gaussformulas{1}, trial_k, test_k);
interE = ...
    matInterEdgeAverJumpParabolic(Coeffs_func, meshInfo,option,Gaussformulas{3},trial_k,test_k);

DirichletE = ...
    matDirichletEdgeAverJumpParabolic(Coeffs_func, meshInfo, option, Gaussformulas{3}, trial_k, test_k);


sysInfo.AssembleTime = toc; %<<<<<<<<<<<<<< corresponding tic2
disp(['assemble time: ',num2str(sysInfo.AssembleTime)])

% system matrix
A = A_uv + interE + DirichletE;

% get the L2 projection of the initial analytical solution
r0 = M\rhs0;


%% solve the system
tic; %<<<<<<<<<<<<<< tic3

% LOOP
% Because the load vector depends on the time, we must write it into the time cycle
sysM = M + option.theta*Dt*A;
right_term = M + (option.theta-1)*Dt*A;
clear A M

n1=1; % control the disp()
for n = 0:NT-1
    if n==0 || n>=n1
        dispname1 = ['in the ', num2str(n),'-th time iteration, t=',num2str(n*Dt)];
        disp(dispname1)
        n1 = n1 + NT/10;
    end % if
    
    t_1 = n*Dt; % current time
    t_2 = (n+1)*Dt; % next time
    t = option.theta*t_2 + (1-option.theta)*t_1;
    
    rhs_ft = vecElemRhsFtParabolic(t, pde.f, meshInfo, Gaussformulas{1}, test_k);
    rhs_DiriBC = vecEdgeRhsDirichletParabolic(t, Coeffs_func, pde.gD, ...
        meshInfo, option, Gaussformulas{3}, test_k);
    
    Rhs = right_term*r0 + Dt*(rhs_ft + rhs_DiriBC);
    Uh = sysM\Rhs;
    clear rhs_ft rhs_DiriBC Rhs
    
    % reassign the r0
    r0 = Uh;
end % for n




sysInfo.SolveSystemTime = toc; %<<<<<<<<<<<<<< corresponding tic3
disp(['solve the system matrix time: ',num2str(sysInfo.SolveSystemTime)])

%%
sysInfo.dof_u = dof_u;
sysInfo.Gaussformulas = Gaussformulas; 
sysInfo.terminalT = option.startingT + NT*Dt;
disp(['terminal Time in program:',num2str(sysInfo.terminalT)])

end % function
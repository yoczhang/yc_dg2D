function [Uh, sysInfo, system] = dgPoissonSolve(meshInfo,pde,option)
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
%   YcZhang 14/5/2017
%
%   Last modified 14/5/2017
%


if ~exist('option','var')
    option = dgOption(option);
    export = 0;
else
    export = option.export;
        %> in the dgOption function, option.exoprt default is 0(i.e FALSE).
        %> if the export is true, then export the variables: system, solvers. 
end

trial_k = basesType2degreek(option.basesType_trial);
test_k = basesType2degreek(option.basesType_test);
dof_u = meshInfo.Nelems * (trial_k+2)*(trial_k+1)/2;


%%  generate the Gaussformulas
% 2d 
[Q1,Q2,W] = quadRule2D(2*trial_k+1); % Gauss-Points on Ref elem
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

%% assemeble matrix and rhs
tic; %<<<<<<<<<<<<<< tic2
Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.k_11;
Coeffs_func{2} = pde.k_22;
Coeffs_func{3} = pde.func0;

[A, rhs_fh] = ...
    matElemCoeffsDtrialDtestPoisson(Coeffs_func,pde.f,meshInfo,Gaussformulas{1}, trial_k, test_k);
interE = ...
    matInterEdgeAverJumpPoisson(Coeffs_func, meshInfo,option,Gaussformulas{3},trial_k,test_k);
[DirichletE, rhs_uDterm] = ...
    matDirichletEdgeAverJumpPoisson(Coeffs_func,pde.gD,meshInfo,option,Gaussformulas{3},trial_k,test_k);

sysM = A + interE + DirichletE;
rhs = rhs_fh + rhs_uDterm;

% %----------------
% A1_yc = A1;
% A2_yc = A2;
% A3_yc = -interE1+interE2+interE3-dirE1+dirE2+dirE3;
% A3_bd_yc = -dirE1+dirE2+dirE3;
% b1_yc = rhs_fh{1};
% b2_yc = rhs_uDterm{1};
% 
% save M_yc A1_yc A2_yc A3_yc b1_yc b2_yc A3_bd_yc
% %----------------

sysInfo.AssembleTime = toc; %<<<<<<<<<<<<<< corresponding tic2
disp(['assemble time: ',num2str(sysInfo.AssembleTime)])

%% solve the system
tic; %<<<<<<<<<<<<<< tic3
if export
    system = {sysM, rhs};
    Uh = [];
    return
else
    system = [];
end 

Uh = sysM\rhs;

sysInfo.SolveSystemTime = toc; %<<<<<<<<<<<<<< corresponding tic3
disp(['solve the system matrix time: ',num2str(sysInfo.SolveSystemTime)])

%%
sysInfo.dof_u = dof_u;
sysInfo.Gaussformulas = Gaussformulas; 

end % function
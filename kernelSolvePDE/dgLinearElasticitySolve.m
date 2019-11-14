function [Uh, sysInfo, system] = dgLinearElasticitySolve(meshInfo,pde,option)
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
%   YcZhang 22/9/2017
%
%   Last modified 22/9/2017
%


if ~exist('option','var')
    option = dgOption(option);
    export = 0;
else
    export = option.export;
        %> in the dgOption function, option.exoprt default is 0(i.e FALSE).
        %> if the export is true, then export the variables: system, solvers. 
end

lambda = pde.lambda;
mu = pde.mu;

degree_u = basesType2degreek(option.basesType_u);
dof_u1 = meshInfo.Nelems * (degree_u+2)*(degree_u+1)/2;

%%  generate the Gaussformulas
% 2d 
[Q1,Q2,W] = quadRule2D(2*degree_u+5); % Gauss-Points on Ref elem
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
tic; %<<<<<<<<<<<<<< tic 1

Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.funcOne;
Coeffs_func{2} = pde.funcOne;
Coeffs_func{3} = pde.funcZero;

[elem_uxvx,elem_uxvy,elem_uyvx,elem_uyvy,rhsf1,rhsf2] = ...
    matElemCoeffsDuDvLinearElasticity(Coeffs_func,pde.f1,pde.f2, meshInfo,Gaussformulas{1}, degree_u);

[Inter_u1v1, Inter_u1v2, Inter_u2v1, Inter_u2v2] = ...
    matInterEdgeAverJumpLinearElasticity(Coeffs_func, pde, meshInfo, option, formulaGauss1D, degree_u);

[Diri_u1v1,Diri_u1v2,Diri_u2v1,Diri_u2v2,rhs_gD1,rhs_gD2] = ...
    matDirichletEdgeAverJumpLinearElasticity(Coeffs_func, pde, pde.gD1, pde.gD2, meshInfo, option, formulaGauss1D, degree_u);

elem_u1v1 = (lambda+2*mu)*elem_uxvx+mu*elem_uyvy;
elem_u1v2 = lambda*elem_uxvy+mu*elem_uyvx;
elem_u2v1 = mu*elem_uxvy+lambda*elem_uyvx;
elem_u2v2 = mu*elem_uxvx + (lambda+2*mu)*elem_uyvy;

sysM = [(elem_u1v1+Inter_u1v1+Diri_u1v1),      (elem_u2v1+Inter_u2v1+Diri_u2v1);
    (elem_u1v2+Inter_u1v2+Diri_u1v2),     (elem_u2v2+Inter_u2v2+Diri_u2v2)];

% sysM = [nu*(elem_u1v1+Inter_u1v1+Diri_u1v1),      nu*(elem_u2v1+Inter_u2v1+Diri_u2v1),     -elem_pv1+Inter_pv1+Diri_pv1;
%     nu*(elem_u1v2+Inter_u1v2+Diri_u1v2),     nu*(elem_u2v2+Inter_u2v2+Diri_u2v2),      -elem_pv2+Inter_pv2+Diri_pv2;
%     -(-elem_pv1+Inter_pv1+Diri_pv1)',     -(-elem_pv2+Inter_pv2+Diri_pv2)',     zeroPP];

Rhs = [rhsf1+rhs_gD1;
    rhsf2+rhs_gD2];

% % save elemMatInfo lambda mu elem_uxvx elem_uxvy elem_uyvx elem_uyvy rhsf1 rhsf2
% % save interEdgeMatInfo Inter_u1v1 Inter_u1v2 Inter_u2v1 Inter_u2v2
% % save DirEdgeMatInfo Diri_u1v1 Diri_u1v2 Diri_u2v1 Diri_u2v2 rhs_gD1 rhs_gD2
% % save systemMatInfo elem_u1v1 elem_u1v2 elem_u2v1 elem_u2v2 sysM Rhs

% disp the condition number of the system matrix
disp(['condition number of sysM: ', num2str(condest(sysM))])

sysInfo.AssembleTime = toc; %<<<<<<<<<<<<<< corresponding tic1
disp(['assemble time: ',num2str(sysInfo.AssembleTime)])

%% solve the system
tic; %<<<<<<<<<<<<<< tic2
if export
    system = {sysM, Rhs};
    Uh = [];
    return
else
    system = [];
end 

Uh = sysM\Rhs;

sysInfo.SolveSystemTime = toc; %<<<<<<<<<<<<<< corresponding tic2
disp(['solve the system matrix time: ',num2str(sysInfo.SolveSystemTime)])

%%
sysInfo.dof_u1 = dof_u1;
sysInfo.Gaussformulas = Gaussformulas; 

end % function
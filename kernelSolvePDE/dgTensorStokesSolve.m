function [Uh, sysInfo, system] = dgTensorStokesSolve(meshInfo,pde,option)
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
%   YcZhang 25/5/2017
%
%   Last modified 26/5/2017
%


if ~exist('option','var')
    option = dgOption(option);
    export = 0;
else
    export = option.export;
        %> in the dgOption function, option.exoprt default is 0(i.e FALSE).
        %> if the export is true, then export the variables: system, solvers. 
end

mu = pde.mu;

degree_u = basesType2degreek(option.basesType_u);
degree_p = basesType2degreek(option.basesType_p);
dof_u1 = meshInfo.Nelems * (degree_u+2)*(degree_u+1)/2;
dof_p = meshInfo.Nelems * (degree_p+2)*(degree_p+1)/2;

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
% interEdgeIndex = meshInfo.interEdgeIndex; % [Ninter x 1]
% DirichletEdgeIndex = meshInfo.bdEdgeIndex; % [Ndir x 1]

meshInfo = getTensorStokesBoundaryInfo(meshInfo);


if ~isfield(meshInfo, 'interEdgeIndex')
    meshInfo.interEdgeIndex = interEdgeIndex;
end

if ~isfield(meshInfo, 'DirichletEdgeIndex')
    meshInfo.DirichletEdgeIndex = DirichletEdgeIndex;
end

%% assemeble matrix and rhs
tic; %<<<<<<<<<<<<<< tic 1
% zeroUU = sparse(dof_u1,dof_u1);
% zeroUP = sparse(dof_u1,dof_p);
zeroPP = sparse(dof_p,dof_p);


Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.k_11;
Coeffs_func{2} = pde.k_22;
Coeffs_func{3} = pde.funcZero;

[elem_uxvx,elem_uxvy,elem_uyvx,elem_uyvy,elem_pv1,elem_pv2,rhsf1,rhsf2,lambda_p] = ...
    matElemCoeffsDuDvTensorStokes(Coeffs_func,pde.f1,pde.f2, meshInfo,Gaussformulas{1}, degree_u, degree_p);

[Inter_u1v1, Inter_u1v2, Inter_u2v1, Inter_u2v2, Inter_pv1, Inter_pv2] = ...
    matInterEdgeAverJumpTensorStokes(Coeffs_func, meshInfo, option, formulaGauss1D, degree_u, degree_p);

[Diri_u1v1,Diri_u1v2,Diri_u2v1,Diri_u2v2,Diri_pv1,Diri_pv2,rhs_gD1,rhs_gD2,rhs_div] = ...
    matDirichletEdgeAverJumpTensorStokes(Coeffs_func, pde.gD1, pde.gD2, meshInfo, option, formulaGauss1D, degree_u, degree_p);

 [Krhs_gN1,Krhs_gN2] = ...
    vecNeumannEdgeAverJumpTensorStokes(pde, meshInfo, formulaGauss1D, degree_u);

elem_u1v1 = 2*(elem_uxvx + 0.5*elem_uyvy);
elem_u1v2 = 2*(0.5*elem_uyvx);
elem_u2v1 = 2*(0.5*elem_uxvy);
elem_u2v2 = 2*(0.5*elem_uxvx + elem_uyvy);

sysM = [mu*(elem_u1v1+Inter_u1v1+Diri_u1v1),      mu*(elem_u2v1+Inter_u2v1+Diri_u2v1),     -elem_pv1+Inter_pv1+Diri_pv1,     zeros(dof_u1,1);
    mu*(elem_u1v2+Inter_u1v2+Diri_u1v2),     mu*(elem_u2v2+Inter_u2v2+Diri_u2v2),      -elem_pv2+Inter_pv2+Diri_pv2,     zeros(dof_u1,1);
    -(-elem_pv1+Inter_pv1+Diri_pv1)',     -(-elem_pv2+Inter_pv2+Diri_pv2)',     zeroPP,     lambda_p;
    zeros(1, dof_u1),       zeros(1, dof_u1),       lambda_p',      0];

Rhs = [rhsf1+mu*rhs_gD1+Krhs_gN1;
    rhsf2+mu*rhs_gD2+Krhs_gN2;
    -rhs_div;
    0];

% sysM = [nu*(elem_u1v1+Inter_u1v1+Diri_u1v1),      nu*(elem_u2v1+Inter_u2v1+Diri_u2v1),     -elem_pv1+Inter_pv1+Diri_pv1;
%     nu*(elem_u1v2+Inter_u1v2+Diri_u1v2),     nu*(elem_u2v2+Inter_u2v2+Diri_u2v2),      -elem_pv2+Inter_pv2+Diri_pv2;
%     -(-elem_pv1+Inter_pv1+Diri_pv1)',     -(-elem_pv2+Inter_pv2+Diri_pv2)',     zeroPP];



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
sysInfo.dof_p = dof_p;
sysInfo.Gaussformulas = Gaussformulas; 

end % function




%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function meshInfo = getTensorStokesBoundaryInfo(meshInfo)
%
%
%
%
%
%   YcZhang 24/10/2017
%
%   Last modified 24/10/2017
%


% the domain is [0,1]x[0,1]
bdEdgeIndex = meshInfo.bdEdgeIndex; % here, all the bdEdge 
DirichletEdgeIndex_1 = bdEdgeIndex( abs(meshInfo.baryEdge(bdEdgeIndex,1)-0) <= 5e-8 ); % left bd
DirichletEdgeIndex_2 = bdEdgeIndex( abs(meshInfo.baryEdge(bdEdgeIndex,2)-0) <= 5e-8 ); % bottom bd
meshInfo.DirichletEdgeIndex = union(DirichletEdgeIndex_1,DirichletEdgeIndex_2);
meshInfo.DirichletEdgeIndex = bdEdgeIndex;

% % % there is no Neumann boudary in the Stokes/Navier-Stokes problem
NeumannEdgeIndex_1 = bdEdgeIndex( abs(meshInfo.baryEdge(bdEdgeIndex,2)-1) <= 5e-8 ); % top bd
NeumannEdgeIndex_2 = bdEdgeIndex( abs(meshInfo.baryEdge(bdEdgeIndex,1)-1) <= 5e-8 ); % right bd
meshInfo.NeumannEdgeIndex = union(NeumannEdgeIndex_1,NeumannEdgeIndex_2);
meshInfo.NeumannEdgeIndex = [];
end % function 
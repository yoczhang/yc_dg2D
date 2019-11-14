function [Uh, sysInfo] = dgVaryingDiffusivityDarcySolve(meshInfo,pde,option)
%
%   %---------------------------------------
%       dgVaryingDiffusivityDarcySolve
%       This test case from: Discontinuous Galerkin methods for
%       Varying Diffusivity Darcy equations for varying and vanishing
%       diffusivity, the 1-th numericak case.
%       domain: [0,2]x[0,1], and left boundary is the inflow, 
%       right boundary is the outflow, bottom and top have no flow.
%   %---------------------------------------
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
%   YcZhang 2/10/2017
%
%   Last modified 2/10/2017
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

% bd_case = 'allDirichlet';
bd_case = 'partDirichlet';
meshInfo = getVaryingDiffusivityDarcyBoundaryInfo(meshInfo,bd_case);


%% assemeble matrix and rhs
tic; %<<<<<<<<<<<<<< tic2
Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.k_11;
Coeffs_func{2} = pde.k_22;


%   %-----------------------------------------------------
%       here the notations according the 'parabolic problems with convection' from Beatrice Riviere's DG book.
%       the init matrix and rhs(t==0)
%       here we use the prefix 'D_' denote the matrix of diffusion-term, and 'C_'
%       denote the matrix of covction-term.
%   %-----------------------------------------------------
[A, rhs_fh] = ...
    matElemCoeffsDtrialDtestVaryingDiffusivityDarcy(Coeffs_func,pde.funcZero,meshInfo,Gaussformulas{1}, trial_k, test_k);

if option.Adap == 1
%     interE = ...
%         matInterEdgeAverJumpAdapTransport(Coeffs_func, meshInfo,option,Gaussformulas{3},trial_k,test_k);
    interE = ...
        matInterEdgeAverJumpVaryingDiffusivityDarcy(Coeffs_func, meshInfo,option,Gaussformulas{3},trial_k,test_k);
else
    interE = ...
        matInterEdgeAverJumpVaryingDiffusivityDarcy_noAdap(Coeffs_func, meshInfo,option,Gaussformulas{3},trial_k,test_k);
end

[DirichletE, rhs_uDterm] = ...
    matDirichletEdgeAverJumpVaryingDiffusivityDarcy(Coeffs_func,pde.funcOne,meshInfo,option,Gaussformulas{3},trial_k,test_k);

sysM = A + interE + DirichletE;
rhs = rhs_fh + rhs_uDterm;

sysInfo.AssembleTime = toc; %<<<<<<<<<<<<<< corresponding tic2
disp(['assemble time: ',num2str(sysInfo.AssembleTime)])


%% solve the system
tic; %<<<<<<<<<<<<<< tic3
Uh = sysM\rhs;

sysInfo.SolveSystemTime = toc; %<<<<<<<<<<<<<< corresponding tic3
disp(['solve the system matrix time: ',num2str(sysInfo.SolveSystemTime)])

%%
sysInfo.dof_u = dof_u;
sysInfo.Gaussformulas = Gaussformulas; 

end % function





%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function meshInfo = getVaryingDiffusivityDarcyBoundaryInfo(meshInfo,bd_case)
%
%
%
%
%
%   YcZhang 2/10/2017
%
%   Last modified 2/10/2017
%

if strcmpi(bd_case,'allDirichlet')
    interEdgeIndex = meshInfo.interEdgeIndex; % [Ninter x 1]
    DirichletEdgeIndex = meshInfo.bdEdgeIndex; % [Ndir x 1]

    meshInfo.interEdgeIndex = interEdgeIndex;
    meshInfo.DirichletEdgeIndex = DirichletEdgeIndex;
end % if

if strcmpi(bd_case,'partDirichlet')
    % beacuse we have the inflow vector_u = [1; 0.5], and the domain is
    % [0,2] x [0,1], so we let:
    % y=0 and x=0 are the inflow boundary,
    % y=1 and x=1 are the outflow boundary.
    
    interEdgeIndex = meshInfo.interEdgeIndex; % [Ninter x 1]
    
    bdEdgeIndex = meshInfo.bdEdgeIndex;
    
    DirichletEdgeIndex_1 = bdEdgeIndex( abs(meshInfo.baryEdge(bdEdgeIndex,1)-0) <= 5e-8 );
    DirichletEdgeIndex_2 = bdEdgeIndex( abs(meshInfo.baryEdge(bdEdgeIndex,1)-2) <= 5e-8 );
    DirichletEdgeIndex = [DirichletEdgeIndex_1; DirichletEdgeIndex_2];
    
    meshInfo.interEdgeIndex = interEdgeIndex;
    meshInfo.DirichletEdgeIndex = DirichletEdgeIndex;
    
end % if


end % function 
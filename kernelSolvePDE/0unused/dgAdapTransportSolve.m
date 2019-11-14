function [Uh, sysInfo] = dgAdapTransportSolve(meshInfo,pde,option)
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
%   YcZhang 23/9/2017
%
%   Last modified 23/9/2017
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

bd_case = 'allDirichlet';
% bd_case = 'haveInflow';
meshInfo = getTransportBoundaryInfo(meshInfo,bd_case);

%% Dt setting
Dt = option.Dt;
NT = floor((option.terminalT-option.startingT)/Dt);
if NT <= 1
    NT = 2;
end

% %--------------------------------
% h = sum(meshInfo.hElem)/meshInfo.Nelems;
% if option.theta == 1/2
%     Dt = h^((trial_k+1)/2);
% else
%     Dt = h^(trial_k+1);
% end
% NT = floor((option.terminalT-option.startingT)/Dt);
% %--------------------------------

disp(['Dt = ',num2str(Dt)])
disp(['NT = ',num2str(NT)])

%% assemeble matrix and rhs
tic; %<<<<<<<<<<<<<< tic2
Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.k_11;
Coeffs_func{2} = pde.k_22;
Coeffs_func{3} = pde.vector_u_1;
Coeffs_func{4} = pde.vector_u_2;

%   %-----------------------------------------------------
%       here the notations according the 'parabolic problems with convection' from Beatrice Riviere's DG book.
%       the init matrix and rhs(t==0)
%       here we use the prefix 'D_' denote the matrix of diffusion-term, and 'C_'
%       denote the matrix of covction-term.
%   %-----------------------------------------------------
[M, D_elem_dudv, C_elem_udv, rhs0] = ...
    matElemCoeffsDtrialDtestAdapTransport(Coeffs_func,pde.u0,meshInfo,Gaussformulas{1}, trial_k, test_k);
[D_interE, D_upwindE, C_interE_upwind] = ...
    matInterEdgeAverJumpAdapTransport(Coeffs_func, meshInfo,option,Gaussformulas{3},trial_k,test_k);

if strcmpi(bd_case,'allDirichlet')
    [D_matDirichletEdge, C_matDirichletEdge] = ...
        matDirichletEdgeAverJumpAdapTransport(Coeffs_func, meshInfo, option, Gaussformulas{3}, trial_k, test_k);
    
    % system matrix
    A = D_elem_dudv + D_interE + D_upwindE + D_matDirichletEdge - C_elem_udv + C_interE_upwind + C_matDirichletEdge;
    clear D_elem_dudv D_interE D_matDirichletEdge C_elem_udv C_interE_upwind C_matDirichletEdge
else
    C_matOutflow = ...
        matOutflowEdgeAdapTransport(Coeffs_func, meshInfo, Gaussformulas{3}, trial_k, test_k);
    % system matrix
    A = D_elem_dudv + D_interE + D_upwindE - C_elem_udv + C_interE_upwind + C_matOutflow;
    %clear D_elem_dudv D_interE C_elem_udv C_interE_upwind C_matOutflow
end

sysInfo.AssembleTime = toc; %<<<<<<<<<<<<<< corresponding tic2
disp(['assemble time: ',num2str(sysInfo.AssembleTime)])



% get the L2 projection of the initial analytical solution
r0 = M\rhs0;
clear rhs0


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
    
    rhs_ft = vecElemRhsFtAdapTransport(t, pde.f_rhs, meshInfo, Gaussformulas{1}, test_k);
    
    if strcmpi(bd_case,'allDirichlet')
        rhs_DiriBC = vecEdgeRhsDirichletAdapTransport(t, Coeffs_func, pde.gD, ...
            meshInfo, option, Gaussformulas{3}, test_k);
        Rhs = right_term*r0 + Dt*(rhs_ft + rhs_DiriBC);
    else
        rhs_Inflow = vecInflowEdgeAdapTransport(t, Coeffs_func, meshInfo, pde, Gaussformulas{3}, test_k);
        Rhs = right_term*r0 + Dt*(rhs_ft - rhs_Inflow);
    end % if
    

    Uh = sysM\Rhs;
    clear rhs_ft rhs_DiriBC rhs_Inflow Rhs
    
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





%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function meshInfo = getTransportBoundaryInfo(meshInfo,bd_case)
%
%
%
%
%
%   YcZhang 5/9/2017
%
%   Last modified 5/9/2017
%

if strcmpi(bd_case,'allDirichlet')
    interEdgeIndex = meshInfo.interEdgeIndex; % [Ninter x 1]
    DirichletEdgeIndex = meshInfo.bdEdgeIndex; % [Ndir x 1]

    meshInfo.interEdgeIndex = interEdgeIndex;
    meshInfo.DirichletEdgeIndex = DirichletEdgeIndex;
end % if

if strcmpi(bd_case,'haveInflow')
    % beacuse we have the inflow vector_u = [1; 0.5], and the domain is
    % [0,1] x [0,1], so we let:
    % y=0 and x=0 are the inflow boundary,
    % y=1 and x=1 are the outflow boundary.
    
    bdEdgeIndex = meshInfo.bdEdgeIndex;
    
    inflowEdgeIndex_1 = bdEdgeIndex( meshInfo.baryEdge(bdEdgeIndex,1)==1 );
    inflowEdgeIndex_2 = bdEdgeIndex( meshInfo.baryEdge(bdEdgeIndex,2)==1 );
    inflowEdgeIndex = [inflowEdgeIndex_1; inflowEdgeIndex_2];
    
    outflowEdgeIndex_1 = bdEdgeIndex( meshInfo.baryEdge(bdEdgeIndex,1)==0 );
    outflowEdgeIndex_2 = bdEdgeIndex( meshInfo.baryEdge(bdEdgeIndex,2)==0 );
    outflowEdgeIndex = [outflowEdgeIndex_1; outflowEdgeIndex_2];
    
    meshInfo.inflowEdgeIndex = inflowEdgeIndex;
    meshInfo.outflowEdgeIndex = outflowEdgeIndex;
    
end % if


end % function 
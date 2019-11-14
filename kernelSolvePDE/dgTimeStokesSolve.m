function [Uh, sysInfo] = dgTimeStokesSolve(StokesmeshInfo,pde,option)
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
%	YcZhang 08/27/2018
%
%   Last modified 08/27/2018
%


if ~exist('option','var')
    option = dgOption(option);
end

mu = pde.mu;

degree_u = basesType2degreek(option.basesType_u);
degree_p = basesType2degreek(option.basesType_p);
dof_u1 = StokesmeshInfo.Nelems * (degree_u+2)*(degree_u+1)/2;
dof_p = StokesmeshInfo.Nelems * (degree_p+2)*(degree_p+1)/2;

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
bd_case = 'allDirichlet';
StokesmeshInfo = ...
    getStokesBoundaryInfo(StokesmeshInfo, bd_case);

%% Dt setting
% Dt = option.Dt;
% NT = floor((option.terminalT-option.startingT)/Dt);
% if NT <= 1
%     NT = 2;
% end

%--------------------------------
h = sum(StokesmeshInfo.hElem)/StokesmeshInfo.Nelems;
if option.theta == 1/2
    Dt = h^((degree_u+1)/2);
else
    Dt = h^(degree_u+1);
end
if Dt >= (option.terminalT-option.startingT) % here we need the (NT = (...)/Dt) >=2.
	Dt = (option.terminalT-option.startingT)/2;
end

NT = floor((option.terminalT-option.startingT)/Dt);
if NT <= 1
    NT = 2;
end
%--------------------------------

disp(['Dt = ',num2str(Dt)])
disp(['NT = ',num2str(NT)])

%% assemeble matrix and rhs
zeroUU = sparse(dof_u1,dof_u1);
% zeroUP = sparse(dof_u1,dof_p);
zeroPP = sparse(dof_p,dof_p);

Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.funcOne;
Coeffs_func{2} = pde.funcOne;
Coeffs_func{3} = pde.funcZero;

tic %<<<<<<<<<<<<<< tic1, assemble StokesDarcy matrix

%% tensor Stokes and interface for Stokes
% tensor Stokes
[elem_u0v0,elem_p0q0,elem_u1v1,elem_pv1,elem_pv2,rhs_u10,rhs_u20,rhs_p0,lambda_p] = ...
    matElemDuDvTimeStokes(Coeffs_func,pde.u10,pde.u20, pde.p0, StokesmeshInfo,Gaussformulas{1}, degree_u, degree_p);

[Inter_u1v1, Inter_pv1, Inter_pv2] = ...
    matInterEdgeAverJumpTimeStokes(Coeffs_func, StokesmeshInfo, option, formulaGauss1D, degree_u, degree_p);

[Diri_u1v1,Diri_pv1,Diri_pv2,~,~,~] = ...
    matDirEdgeAverJumpTimeStokes(Coeffs_func, pde.gD1, pde.gD2, StokesmeshInfo, option, formulaGauss1D, degree_u, degree_p);

% the Stokes system matrix and rhs
Stokes_stiff = [mu*(elem_u1v1+Inter_u1v1+Diri_u1v1),      zeroUU,     -elem_pv1+Inter_pv1+Diri_pv1,     zeros(dof_u1,1);
    zeroUU,     mu*(elem_u1v1+Inter_u1v1+Diri_u1v1),      -elem_pv2+Inter_pv2+Diri_pv2,     zeros(dof_u1,1);
    -(-elem_pv1+Inter_pv1+Diri_pv1)',     -(-elem_pv2+Inter_pv2+Diri_pv2)',     zeroPP,     lambda_p;
    zeros(1, dof_u1),       zeros(1, dof_u1),       lambda_p',      0];

% % StokesRhs = [rhs_u10+mu*rhs_gD1;
% %     rhs_u20+mu*rhs_gD2;
% %     -rhs_div];

%% get the StokeDarcy system matrix and rhs
Stokes_mass = sparse(2*dof_u1+dof_p+1, 2*dof_u1+dof_p+1);

Stokes_mass(1:dof_u1, 1:dof_u1) = elem_u0v0;
Stokes_mass(dof_u1+1:2*dof_u1, dof_u1+1:2*dof_u1) = elem_u0v0;


%%
sysInfo.AssembleTime = toc; %<<<<<<<<<<<<<< corresponding tic1
disp(['assemble time: ',num2str(sysInfo.AssembleTime)])

%% solve the system
tic; %<<<<<<<<<<<<<< tic2

% LOOP
% Because the load vector depends on the time, we must write it into the time cycle
u1h_0 = elem_u0v0\rhs_u10; 
u2h_0 = elem_u0v0\rhs_u20; 
ph_0 = elem_p0q0\rhs_p0; 
Uh0 = [u1h_0;u2h_0;ph_0;0];

sysM = Stokes_mass + option.theta*Dt*Stokes_stiff;
right_term = Stokes_mass + (option.theta-1)*Dt*Stokes_stiff;

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
    
    [S_rhs_f1, S_rhs_f2] = vecElemRhsFtTimeStokes(t, pde.f1, pde.f2, StokesmeshInfo, Gaussformulas{1}, degree_u);
    vecRhs_f = [S_rhs_f1;S_rhs_f2;zeros(dof_p,1);0];
    
    if strcmpi(bd_case,'allDirichlet')
        [S_rhs_gD1, S_rhs_gD2, S_rhs_div] = vecEdgeRhsDirTimeStokes(t, Coeffs_func, pde.gD1, pde.gD2, ...
            StokesmeshInfo, option, Gaussformulas{3}, degree_u, degree_p);
        rhs_DiriBC = [mu*S_rhs_gD1;mu*S_rhs_gD2;-S_rhs_div;0];
        
        S_sysRhs = right_term*Uh0 + Dt*(vecRhs_f + rhs_DiriBC);
    else
        error('error in dgTimeStokesDarcySolve_usingTensorStokes.m')
    end % if
    

    Uh = sysM\S_sysRhs; 
    
    % reassign the r0
    Uh0 = Uh;
end % for n

Uh = Uh(1:end-1);

sysInfo.SolveSystemTime = toc; %<<<<<<<<<<<<<< corresponding tic2
disp(['solve the system matrix time: ',num2str(sysInfo.SolveSystemTime)])

%%
sysInfo.dof_u1 = dof_u1;
sysInfo.dof_p = dof_p;
sysInfo.Gaussformulas = Gaussformulas;
sysInfo.terminalT = option.startingT + NT*Dt;
disp(['terminal Time in program:',num2str(sysInfo.terminalT)])

end % function


%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function StokesmeshInfo = ...
    getStokesBoundaryInfo(StokesmeshInfo, bd_case)
%
%
%
%
%
%	YcZhang 08/23/2018
%
%   Last modified 08/23/2018
%

if strcmpi(bd_case, 'allDirichlet')
    %> Stokes domain: [0,1]x[0,1];
    
    %% bdEdge setting
    % Stokes domain bdEdge setting
    S_DirichletEdgeIndex = StokesmeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 

    StokesmeshInfo.DirichletEdgeIndex = S_DirichletEdgeIndex;
    
else
    error('there is no case in getStokesBoundaryInfo ');
    
end % if 


end % function 


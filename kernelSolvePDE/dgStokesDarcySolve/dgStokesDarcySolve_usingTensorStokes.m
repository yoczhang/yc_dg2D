function [Uh, sysInfo, system] = dgStokesDarcySolve_usingTensorStokes(StokesmeshInfo,DarcymeshInfo,pde,option)
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
%   YcZhang 28/5/2017
%
%   Last modified 28/5/2017
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

degree_phi = basesType2degreek(option.basesType_Pp);
degree_u = basesType2degreek(option.basesType_Fu);
degree_p = basesType2degreek(option.basesType_Fp);
dof_phi = DarcymeshInfo.Nelems * (degree_phi+2)*(degree_phi+1)/2;
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
[StokesmeshInfo, DarcymeshInfo, interfacemeshInfo] = ...
    getStokesDarcyBoundaryInfo(StokesmeshInfo, DarcymeshInfo, pde);

% plot_interface(StokesmeshInfo, DarcymeshInfo);

%% assemeble matrix and rhs
tic; %<<<<<<<<<<<<<< tic 1

% zeroUU = sparse(dof_u1,dof_u1);
% zeroUP = sparse(dof_u1,dof_p);
zeroPP = sparse(dof_p,dof_p);

Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.funcOne;
Coeffs_func{2} = pde.funcOne;
Coeffs_func{3} = pde.funcZero;

Kappa = pde.Kappa; % this is the coefficient of (P(uh),P(vh))_interface


%% tensor Stokes and interface for Stokes
% tensor Stokes
[elem_uxvx,elem_uxvy,elem_uyvx,elem_uyvy,elem_pv1,elem_pv2,rhsf1,rhsf2,~] = ...
    matElemCoeffsDuDvTensorStokes_SD(Coeffs_func,pde.f1,pde.f2, StokesmeshInfo,Gaussformulas{1}, degree_u, degree_p);

[Inter_u1v1, Inter_u1v2, Inter_u2v1, Inter_u2v2, Inter_pv1, Inter_pv2] = ...
    matInterEdgeAverJumpTensorStokes_SD(Coeffs_func, StokesmeshInfo, option, formulaGauss1D, degree_u, degree_p);

[Diri_u1v1,Diri_u1v2,Diri_u2v1,Diri_u2v2,Diri_pv1,Diri_pv2,rhs_gD1,rhs_gD2,rhs_div] = ...
    matDirichletEdgeAverJumpTensorStokes_SD(Coeffs_func, pde.gD1, pde.gD2, StokesmeshInfo, option, formulaGauss1D, degree_u, degree_p);

elem_u1v1 = 2*(elem_uxvx + 0.5*elem_uyvy);
elem_u1v2 = 2*(0.5*elem_uyvx);
elem_u2v1 = 2*(0.5*elem_uxvy);
elem_u2v2 = 2*(0.5*elem_uxvx + elem_uyvy);

% the interface integration for Stokes 
[interface_u1v1, interface_u2v1, interface_u1v2, interface_u2v2] = ...
    matInterfaceOfStokes(Coeffs_func, StokesmeshInfo, interfacemeshInfo, formulaGauss1D, degree_u);

% the Stokes system matrix and rhs
StokesM = [mu*(elem_u1v1+Inter_u1v1+Diri_u1v1)+Kappa*interface_u1v1,      mu*(elem_u2v1+Inter_u2v1+Diri_u2v1)+Kappa*interface_u2v1,     -elem_pv1+Inter_pv1+Diri_pv1;
    mu*(elem_u1v2+Inter_u1v2+Diri_u1v2)+Kappa*interface_u1v2,     mu*(elem_u2v2+Inter_u2v2+Diri_u2v2)+Kappa*interface_u2v2,      -elem_pv2+Inter_pv2+Diri_pv2;
    -(-elem_pv1+Inter_pv1+Diri_pv1)',     -(-elem_pv2+Inter_pv2+Diri_pv2)',     zeroPP];

StokesRhs = [rhsf1+mu*rhs_gD1;
    rhsf2+mu*rhs_gD2;
    -rhs_div];


%% the Darcy 
[Darcy_A, Darcy_rhs_fh] = ...
    matElemCoeffsDtrialDtestPoisson_SD(Coeffs_func,pde.fp,DarcymeshInfo,Gaussformulas{1}, degree_phi, degree_phi);
Darcy_interE = ...
    matInterEdgeAverJumpPoisson_SD(Coeffs_func, DarcymeshInfo,option,Gaussformulas{3},degree_phi,degree_phi);
[Darcy_DirichletE, Darcy_rhs_uDterm] = ...
    matDirichletEdgeAverJumpPoisson_SD(Coeffs_func,pde.gDp,DarcymeshInfo,option,Gaussformulas{3},degree_phi,degree_phi);
% % rhsNeumann = vecNeumannEdgePoisson_SD(Coeffs_func,pde,DarcymeshInfo,Gaussformulas{3},degree_phi);

DarcyM = Darcy_A + Darcy_interE + Darcy_DirichletE;
% DarcyRhs = Darcy_rhs_fh + Darcy_rhs_uDterm + rhsNeumann;
DarcyRhs = Darcy_rhs_fh + Darcy_rhs_uDterm;

%% the Stoke-Darcy coupled interface integration
[interface_uphi_n1, interface_uphi_n2, interface_phiu_n1, interface_phiu_n2] = ...
    matInterfaceStokesDarcy(Coeffs_func, StokesmeshInfo, DarcymeshInfo, interfacemeshInfo, formulaGauss1D, degree_u, degree_phi);


%% get the StokeDarcy system matrix and rhs
systemM = sparse(dof_phi+2*dof_u1+dof_p, dof_phi+2*dof_u1+dof_p);

systemM(1:dof_phi, 1:dof_phi) = DarcyM;
systemM(1:dof_phi, dof_phi+1:dof_phi+dof_u1) = interface_uphi_n1;
systemM(1:dof_phi, dof_phi+dof_u1+1:dof_phi+2*dof_u1) = interface_uphi_n2;
systemM(dof_phi+1:dof_phi+dof_u1, 1:dof_phi) = interface_phiu_n1;
systemM(dof_phi+dof_u1+1:dof_phi+2*dof_u1, 1:dof_phi) = interface_phiu_n2;
systemM(dof_phi+1:end, dof_phi+1:end) = StokesM;

systemRhs = [DarcyRhs; StokesRhs];

%%
% disp the condition number of the system matrix
disp(['condition number of systemM: ', num2str(condest(systemM))])

sysInfo.AssembleTime = toc; %<<<<<<<<<<<<<< corresponding tic1
disp(['assemble time: ',num2str(sysInfo.AssembleTime)])

%% solve the system
tic; %<<<<<<<<<<<<<< tic2
if export
    system = {systemM, systemRhs};
    Uh = [];
    return
else
    system = [];
end 

Uh = systemM\systemRhs;

sysInfo.SolveSystemTime = toc; %<<<<<<<<<<<<<< corresponding tic2
disp(['solve the system matrix time: ',num2str(sysInfo.SolveSystemTime)])

%%
sysInfo.dof_phi = dof_phi;
sysInfo.dof_u1 = dof_u1;
sysInfo.dof_p = dof_p;
sysInfo.Gaussformulas = Gaussformulas;

end % function


%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function [StokesmeshInfo, DarcymeshInfo, interfacemeshInfo] = ...
    getStokesDarcyBoundaryInfo(StokesmeshInfo, DarcymeshInfo, pde)
%
%
%
%
%
%   YcZhang 28/5/2017
%
%   Last modified 28/5/2017
%

if strcmpi(pde.pdecase,'case1') 
    %> in pde case1, the interface is setted by y==0.
    %% bdEdge setting
    % Stokes domain bdEdge setting
    S_interEdgeIndex = StokesmeshInfo.interEdgeIndex; % [Ninter x 1]

    S_DirichletEdgeIndex = StokesmeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 
    S_interfaceEdgeIndex = S_DirichletEdgeIndex( abs(StokesmeshInfo.baryEdge(S_DirichletEdgeIndex,2)-0)<5e-7 );
    S_DirichletEdgeIndex = setdiff(S_DirichletEdgeIndex, S_interfaceEdgeIndex);

    StokesmeshInfo.interEdgeIndex = S_interEdgeIndex;
    StokesmeshInfo.DirichletEdgeIndex = S_DirichletEdgeIndex;
    StokesmeshInfo.interfaceEdgeIndex = S_interfaceEdgeIndex;

    % Darcy domain bdEdge setting
    D_interEdgeIndex = DarcymeshInfo.interEdgeIndex; % [Ninter x 1]

    D_bdEdgeIndex = DarcymeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 
    D_interfaceEdgeIndex = D_bdEdgeIndex( abs(DarcymeshInfo.baryEdge(D_bdEdgeIndex,2)-0)<5e-8 );
    D_DirichletEdgeIndex = D_bdEdgeIndex( abs(DarcymeshInfo.baryEdge(D_bdEdgeIndex,2)-0)>=5e-8 );
    D_NeumannEdgeIndex = setdiff(D_bdEdgeIndex, union(D_interfaceEdgeIndex,D_DirichletEdgeIndex));

    DarcymeshInfo.interEdgeIndex = D_interEdgeIndex;
    DarcymeshInfo.interfaceEdgeIndex =D_interfaceEdgeIndex;
    DarcymeshInfo.DirichletEdgeIndex = D_DirichletEdgeIndex;
    DarcymeshInfo.NeumannEdgeIndex = D_NeumannEdgeIndex;
    
    %% get the interface infromation
    interfacemeshInfo=getCoupledEdgeInfo(StokesmeshInfo, S_interfaceEdgeIndex, ...
        DarcymeshInfo, D_interfaceEdgeIndex, 'xcoord');
    
elseif strcmpi(pde.pdecase,'case2')
    %> in pde case2, the interface is setted by y==1.
    
    %% bdEdge setting
    % Stokes domain bdEdge setting
    S_interEdgeIndex = StokesmeshInfo.interEdgeIndex; % [Ninter x 1]

    S_DirichletEdgeIndex = StokesmeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 
    S_interfaceEdgeIndex = S_DirichletEdgeIndex( StokesmeshInfo.baryEdge(S_DirichletEdgeIndex,2)==1 );
    S_DirichletEdgeIndex = setdiff(S_DirichletEdgeIndex, S_interfaceEdgeIndex);

    StokesmeshInfo.interEdgeIndex = S_interEdgeIndex;
    StokesmeshInfo.DirichletEdgeIndex = S_DirichletEdgeIndex;
    StokesmeshInfo.interfaceEdgeIndex = S_interfaceEdgeIndex;

    % Darcy domain bdEdge setting
    D_interEdgeIndex = DarcymeshInfo.interEdgeIndex; % [Ninter x 1]

    D_bdEdgeIndex = DarcymeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 
    D_interfaceEdgeIndex = D_bdEdgeIndex( DarcymeshInfo.baryEdge(D_bdEdgeIndex,2)==1 );
    D_DirichletEdgeIndex = D_bdEdgeIndex( DarcymeshInfo.baryEdge(D_bdEdgeIndex,2)~=1 );
    D_NeumannEdgeIndex = setdiff(D_bdEdgeIndex, union(D_interfaceEdgeIndex,D_DirichletEdgeIndex));

    DarcymeshInfo.interEdgeIndex = D_interEdgeIndex;
    DarcymeshInfo.interfaceEdgeIndex =D_interfaceEdgeIndex;
    DarcymeshInfo.DirichletEdgeIndex = D_DirichletEdgeIndex;
    DarcymeshInfo.NeumannEdgeIndex = D_NeumannEdgeIndex;
    
    %% get the interface infromation
    interfacemeshInfo=getCoupledEdgeInfo(StokesmeshInfo, S_interfaceEdgeIndex, ...
        DarcymeshInfo, D_interfaceEdgeIndex, 'xcoord');
    
elseif strcmpi(pde.pdecase,'case4')
        %> in pde case4, the interface is setted by y==1/2.
    
    %% bdEdge setting
    % Stokes domain bdEdge setting
    S_interEdgeIndex = StokesmeshInfo.interEdgeIndex; % [Ninter x 1]

    S_DirichletEdgeIndex = StokesmeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 
    S_interfaceEdgeIndex = S_DirichletEdgeIndex( abs(StokesmeshInfo.baryEdge(S_DirichletEdgeIndex,2)-1/2)<5e-7 );
    S_DirichletEdgeIndex = setdiff(S_DirichletEdgeIndex, S_interfaceEdgeIndex);

    StokesmeshInfo.interEdgeIndex = S_interEdgeIndex;
    StokesmeshInfo.DirichletEdgeIndex = S_DirichletEdgeIndex;
    StokesmeshInfo.interfaceEdgeIndex = S_interfaceEdgeIndex;

    % Darcy domain bdEdge setting
    D_interEdgeIndex = DarcymeshInfo.interEdgeIndex; % [Ninter x 1]

    D_bdEdgeIndex = DarcymeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 
    D_interfaceEdgeIndex = D_bdEdgeIndex( abs(DarcymeshInfo.baryEdge(D_bdEdgeIndex,2)-1/2)<5e-8 );
    D_DirichletEdgeIndex = D_bdEdgeIndex( abs(DarcymeshInfo.baryEdge(D_bdEdgeIndex,2)-1/2)>=5e-8 );
    D_NeumannEdgeIndex = setdiff(D_bdEdgeIndex, union(D_interfaceEdgeIndex,D_DirichletEdgeIndex));

    DarcymeshInfo.interEdgeIndex = D_interEdgeIndex;
    DarcymeshInfo.interfaceEdgeIndex =D_interfaceEdgeIndex;
    DarcymeshInfo.DirichletEdgeIndex = D_DirichletEdgeIndex;
    DarcymeshInfo.NeumannEdgeIndex = D_NeumannEdgeIndex;
    
    %% get the interface infromation
    interfacemeshInfo=getCoupledEdgeInfo(StokesmeshInfo, S_interfaceEdgeIndex, ...
        DarcymeshInfo, D_interfaceEdgeIndex, 'xcoord');
    
end % if 


end % function 




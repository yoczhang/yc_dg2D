b2function [SD_Uh0, Ch, sysInfo, system] = dgDecoupledSDTransportSolve(StokesmeshInfo,DarcymeshInfo,pde,option)
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
%   YcZhang 27/9/2017
%
%   Last modified 27/10/2017
%


if ~exist('option','var')
    option = dgOption(option);
    export = 0;
else
    export = option.export;
        %> in the dgOption function, option.exoprt default is 0(i.e FALSE).
        %> if the export is true, then export the variables: system, solvers. 
end

nu = 1.0;

degree_phi = basesType2degreek(option.basesType_Pp);
degree_u = basesType2degreek(option.basesType_Fu);
degree_p = basesType2degreek(option.basesType_Fp);
degree_c = basesType2degreek(option.basesType_c);
dof_phi = DarcymeshInfo.Nelems * (degree_phi+2)*(degree_phi+1)/2;
dof_u1 = StokesmeshInfo.Nelems * (degree_u+2)*(degree_u+1)/2;
dof_p = StokesmeshInfo.Nelems * (degree_p+2)*(degree_p+1)/2;
dof_Sc = StokesmeshInfo.Nelems * (degree_c+2)*(degree_c+1)/2;
    %> the dof_c on Stokes domain.
dof_Dc = DarcymeshInfo.Nelems * (degree_c+2)*(degree_c+1)/2;
    %> the dof_c on Darcy domain.
dof_c = dof_Sc + dof_Dc;
dof_SD = dof_phi + 2*dof_u1 + dof_p;
% system_dof = dof_SD + dof_c;

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
% zeroUC = sparse(dof_SD,dof_c);

Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.funcOne;
Coeffs_func{2} = pde.funcOne;
Coeffs_func{3} = pde.funcZero;

Ctau = pde.c_tau; % this is the coefficient of (P(uh),P(vh))_interface


%% ---- tensor Stokes and interface for Stokes
% tensor Stokes
[S_elem_uxvx,S_elem_uxvy,S_elem_uyvx,S_elem_uyvy,S_elem_pv1,S_elem_pv2,S_rhsf1,S_rhsf2,~] = ...
    matElemCoeffsDuDvSDTransport_TensorStokes(Coeffs_func,pde.f1,pde.f2, StokesmeshInfo,Gaussformulas{1}, degree_u, degree_p);

[S_Inter_u1v1, S_Inter_u1v2, S_Inter_u2v1, S_Inter_u2v2, S_Inter_pv1, S_Inter_pv2] = ...
    matInterEdgeAverJumpSDTransport_TensorStokes(Coeffs_func, StokesmeshInfo, option, formulaGauss1D, degree_u, degree_p);

[S_Diri_u1v1,S_Diri_u1v2,S_Diri_u2v1,S_Diri_u2v2,S_Diri_pv1,S_Diri_pv2,S_rhs_gD1,S_rhs_gD2,S_rhs_div] = ...
    matDirichletEdgeAverJumpSDTransport_TensorStokes(Coeffs_func, pde.gD1, pde.gD2, StokesmeshInfo, option, formulaGauss1D, degree_u, degree_p);

S_elem_u1v1 = 2*(S_elem_uxvx + 0.5*S_elem_uyvy);
S_elem_u1v2 = 2*(0.5*S_elem_uyvx);
S_elem_u2v1 = 2*(0.5*S_elem_uxvy);
S_elem_u2v2 = 2*(0.5*S_elem_uxvx + S_elem_uyvy);

% the interface integration for Stokes 
[interface_u1v1, interface_u2v1, interface_u1v2, interface_u2v2] = ...
    matInterfaceOfStokes(Coeffs_func, StokesmeshInfo, interfacemeshInfo, formulaGauss1D, degree_u);

% the Stokes system matrix and rhs
StokesM = [nu*(S_elem_u1v1+S_Inter_u1v1+S_Diri_u1v1)+Ctau*interface_u1v1,      nu*(S_elem_u2v1+S_Inter_u2v1+S_Diri_u2v1)+Ctau*interface_u2v1,     -S_elem_pv1+S_Inter_pv1+S_Diri_pv1;
    nu*(S_elem_u1v2+S_Inter_u1v2+S_Diri_u1v2)+Ctau*interface_u1v2,     nu*(S_elem_u2v2+S_Inter_u2v2+S_Diri_u2v2)+Ctau*interface_u2v2,      -S_elem_pv2+S_Inter_pv2+S_Diri_pv2;
    -(-S_elem_pv1+S_Inter_pv1+S_Diri_pv1)',     -(-S_elem_pv2+S_Inter_pv2+S_Diri_pv2)',     zeroPP];

StokesRhs = [S_rhsf1+nu*S_rhs_gD1;
    S_rhsf2+nu*S_rhs_gD2;
    -S_rhs_div];


%% ---- the Darcy 
[Darcy_A, Darcy_rhs_fh] = ...
    matElemCoeffsDtrialDtestSDTransport_Poisson(Coeffs_func,pde.fp,DarcymeshInfo,Gaussformulas{1}, degree_phi, degree_phi);
Darcy_interE = ...
    matInterEdgeAverJumpSDTransport_Poisson(Coeffs_func, DarcymeshInfo,option,Gaussformulas{3},degree_phi,degree_phi);
[Darcy_DirichletE, Darcy_rhs_uDterm] = ...
    matDirichletEdgeAverJumpSDTransport_Poisson(Coeffs_func,pde.gDp,DarcymeshInfo,option,Gaussformulas{3},degree_phi,degree_phi);
% rhsNeumann = vecNeumannEdgePoisson(Coeffs_func,pde,DarcymeshInfo,Gaussformulas{3},degree_phi);

DarcyM = Darcy_A + Darcy_interE + Darcy_DirichletE;
% DarcyRhs = Darcy_rhs_fh + Darcy_rhs_uDterm + rhsNeumann;
DarcyRhs = Darcy_rhs_fh + Darcy_rhs_uDterm;


%% ---- the Stoke-Darcy coupled interface integration
[interface_uphi_n1, interface_uphi_n2, interface_phiu_n1, interface_phiu_n2] = ...
    matSDTransport_InterfaceStokesDarcy(Coeffs_func, StokesmeshInfo, DarcymeshInfo, interfacemeshInfo, formulaGauss1D, degree_u, degree_phi);


%% ---- Transport matrix and Transport coupled Stokes, Darcy matrix
% and we will using 
%   'S' standsfor Stokes equation or domain, 
%   'D' standsfor Darcy equation or domain, 
%   'T' standsfor Transport euqation.

% the Stokes domain



%% ---- get the StokeDarcy system matrix and rhs
systemM_SD = sparse(dof_phi+2*dof_u1+dof_p, dof_phi+2*dof_u1+dof_p);

systemM_SD(1:dof_phi, 1:dof_phi) = DarcyM;
systemM_SD(1:dof_phi, dof_phi+1:dof_phi+dof_u1) = interface_uphi_n1;
systemM_SD(1:dof_phi, dof_phi+dof_u1+1:dof_phi+2*dof_u1) = interface_uphi_n2;
systemM_SD(dof_phi+1:dof_phi+dof_u1, 1:dof_phi) = interface_phiu_n1;
systemM_SD(dof_phi+dof_u1+1:dof_phi+2*dof_u1, 1:dof_phi) = interface_phiu_n2;
systemM_SD(dof_phi+1:end, dof_phi+1:end) = StokesM;

systemRhs_SD = [DarcyRhs; StokesRhs];

%--- get the StokesDarcy Uh0
SD_Uh0 = systemM_SD\systemRhs_SD;

% systemM = sparse(system_dof,system_dof);
lastStep_DarcyUh = SD_Uh0(1:dof_phi);
lastStep_StokesUh = SD_Uh0(dof_phi+1:dof_phi+2*dof_u1);
lastStep_Ch = zeros(dof_c,1);

%--- using the SD_Uh0 to get the: inflow and outflow edges of Transport
[StokesmeshInfo, DarcymeshInfo] = ...
    getTransportBoundaryInfo(StokesmeshInfo, DarcymeshInfo, ...
    lastStep_DarcyUh, lastStep_StokesUh, degree_phi, degree_u);

%% % transport time loop
% Dt setting
Dt = option.Dt;
NT = floor((option.terminalT-option.startingT)/Dt);
if NT <= 1
    NT = 2;
end
disp(['Dt = ',num2str(Dt)])
disp(['NT = ',num2str(NT)])

%
[C_G, C_GxxPlusGyy, C_u1G0xPlusu2G0y, ~, ~, Transp_rhs0] = ...
    matElemCoeffsDtrialDtestSDTransport_Transport(...
    Coeffs_func, pde.c0, DarcymeshInfo, StokesmeshInfo, ...
	lastStep_DarcyUh, lastStep_StokesUh, lastStep_Ch, ...
	Gaussformulas{1}, degree_phi, degree_u, degree_c, degree_c);

[C_diff_K_noAdapE, C_diff_K_AdapE, C_conv_K, ~, ~] = ...
	matInterEdgeAverJumpSDTransport_Transport(...%Coeffs_func, ...
	DarcymeshInfo, StokesmeshInfo, interfacemeshInfo, ...
	lastStep_DarcyUh, lastStep_StokesUh, lastStep_Ch, ...
	option, formulaGauss1D, degree_phi, degree_u, degree_c, degree_c);

[C_K_outflow, ~, ~] = ...
	matOutflowEdgeSDTransport_Transport(... %
	DarcymeshInfo, StokesmeshInfo, ...
	lastStep_DarcyUh, lastStep_StokesUh, lastStep_Ch, ...
	formulaGauss1D, degree_phi, degree_u, degree_c, degree_c ); % need to check the code again, 24/10/2017

Ch0 = C_G\Transp_rhs0;

A = C_GxxPlusGyy + C_diff_K_noAdapE + C_diff_K_AdapE ...
    - C_u1G0xPlusu2G0y + C_conv_K + C_K_outflow;

sysM = C_G + option.theta*Dt*A;
rhs_term = C_G + (option.theta-1)*Dt*A;


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
    
    rhs_ft = vecElemRhsFtSDTransport_Transport(t, ...
        pde, lastStep_DarcyUh, lastStep_StokesUh,...
        DarcymeshInfo, StokesmeshInfo, ...
        formulaGauss2D, degree_phi, degree_u, degree_c);
    
    % --- in flow edges
	%pde.f_inflow = @(t,x,y) 1 + 0.*x;
	rhs_Inflow = vecInflowEdgeSDTransport_Transport(t, pde, ...
        DarcymeshInfo, StokesmeshInfo, ...
        lastStep_DarcyUh, lastStep_StokesUh, ...
        formulaGauss1D, degree_phi, degree_u, degree_c); 
	%max(rhs_Inflow)
	%min(rhs_Inflow)
	sysRhs = rhs_term*Ch0 + Dt*(rhs_ft - rhs_Inflow);
    

    Ch = sysM\sysRhs;
    clear rhs_ft rhs_DiriBC rhs_Inflow Rhs
    
    % reassign the r0
    Ch0 = Ch;
end % for n


%% ---- 
% disp the condition number of the system matrix
disp(['condition number of systemM: ', num2str(condest(systemM_SD))])

sysInfo.AssembleTime = toc; %<<<<<<<<<<<<<< corresponding tic1
disp(['assemble time: ',num2str(sysInfo.AssembleTime)])

%% solve the system
tic; %<<<<<<<<<<<<<< tic2
if export
    system = {systemM_SD, systemRhs_SD};
    Ch = [];
    return
else
    system = [];
end 

Ch = systemM_SD\systemRhs_SD;

sysInfo.SolveSystemTime = toc; %<<<<<<<<<<<<<< corresponding tic2
disp(['solve the system matrix time: ',num2str(sysInfo.SolveSystemTime)])

%%
sysInfo.dof_phi = dof_phi;
sysInfo.dof_u1 = dof_u1;
sysInfo.dof_p = dof_p;
sysInfo.dof_Sc = dof_Sc;
sysInfo.dof_Dc = dof_Dc;
sysInfo.Gaussformulas = Gaussformulas;
sysInfo.terminalT = option.startingT + NT*Dt;

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
    
elseif strcmpi(pde.pdecase,'case2')
    %> in pde case2, the interface is setted by y==1/2.
    
    %% bdEdge setting
    % Stokes domain bdEdge setting
    S_interEdgeIndex = StokesmeshInfo.interEdgeIndex; % [Ninter x 1]

    S_DirichletEdgeIndex = StokesmeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 
    S_interfaceEdgeIndex = S_DirichletEdgeIndex( abs(StokesmeshInfo.baryEdge(S_DirichletEdgeIndex,2)-1/2)<=5e-7 );
    S_DirichletEdgeIndex = setdiff(S_DirichletEdgeIndex, S_interfaceEdgeIndex);
    
    %--- set inflow and outflow edges for Transport eqn.
    S_inflowEdgeIndex = S_DirichletEdgeIndex( abs(StokesmeshInfo.baryEdge(S_DirichletEdgeIndex,1)-0)<=5e-7 );
    S_outflowEdgeIndex = S_DirichletEdgeIndex( abs(StokesmeshInfo.baryEdge(S_DirichletEdgeIndex,1)-1)<=5e-7 );
    
    StokesmeshInfo.interEdgeIndex = S_interEdgeIndex;
    StokesmeshInfo.DirichletEdgeIndex = S_DirichletEdgeIndex;
    StokesmeshInfo.interfaceEdgeIndex = S_interfaceEdgeIndex;
    
    %--- set inflow and outflow edges for Transport eqn.
    StokesmeshInfo.inflowEdgeIndex = S_inflowEdgeIndex;
    StokesmeshInfo.outflowEdgeIndex = S_outflowEdgeIndex;

    % Darcy domain bdEdge setting
    D_interEdgeIndex = DarcymeshInfo.interEdgeIndex; % [Ninter x 1]

    D_bdEdgeIndex = DarcymeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 
    D_interfaceEdgeIndex = D_bdEdgeIndex( abs(DarcymeshInfo.baryEdge(D_bdEdgeIndex,2)-1/2)<= 5e-7 );
    D_DirichletEdgeIndex = D_bdEdgeIndex( abs(DarcymeshInfo.baryEdge(D_bdEdgeIndex,2)-1/2) > 5e-7 );
    D_NeumannEdgeIndex = setdiff(D_bdEdgeIndex, union(D_interfaceEdgeIndex,D_DirichletEdgeIndex));
    
    %--- set inflow and outflow edges for Transport eqn.
    D_inflowEdgeIndex = D_bdEdgeIndex( abs(DarcymeshInfo.baryEdge(D_bdEdgeIndex,1)-0)<= 5e-7 );
    D_outflowEdgeIndex = D_bdEdgeIndex( abs(DarcymeshInfo.baryEdge(D_bdEdgeIndex,1)-1)<= 5e-7 );

    DarcymeshInfo.interEdgeIndex = D_interEdgeIndex;
    DarcymeshInfo.interfaceEdgeIndex =D_interfaceEdgeIndex;
    DarcymeshInfo.DirichletEdgeIndex = D_DirichletEdgeIndex;
    DarcymeshInfo.NeumannEdgeIndex = D_NeumannEdgeIndex;
    
    %--- set inflow and outflow edges for Transport eqn.
    DarcymeshInfo.inflowEdgeIndex = D_inflowEdgeIndex;
    DarcymeshInfo.outflowEdgeIndex = D_outflowEdgeIndex;
    
    %% get the interface infromation
    interfacemeshInfo=getCoupledEdgeInfo(StokesmeshInfo, S_interfaceEdgeIndex, ...
        DarcymeshInfo, D_interfaceEdgeIndex, 'xcoord');
    
    
end % if 


end % function 

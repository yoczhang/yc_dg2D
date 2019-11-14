function [SD_Uh, Ch, sysInfo] = dgTimeDecoupledSDTSolve1(StokesmeshInfo,DarcymeshInfo,pde,option)
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       dgTimeDecoupledSDT1, 
%       for the transport eqn, have NO inflow and outflow boundaryEdges,
%       i.e. for the transport eqn, all the boundaryEdges are Dirichlet Edges.
%   %---------------------------------------------
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
%   YcZhang 18/11/2017
%
%   Last modified 18/11/2017
%


if ~exist('option','var')
    option = dgOption(option);
end

mu = pde.mu;
Kappa = pde.Kappa; % this is the coefficient of (P(uh),P(vh))_interface

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
% dof_c = dof_Sc + dof_Dc;
% dof_SD = dof_phi + 2*dof_u1 + dof_p;
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
bd_case = 'allDirichlet';
[StokesmeshInfo, DarcymeshInfo, interfacemeshInfo] = ...
    getStokesDarcyBoundaryInfo(StokesmeshInfo, DarcymeshInfo, bd_case, pde);

% plot_interface(StokesmeshInfo, DarcymeshInfo);

%% Dt setting
% %-------------- plan 1 ------------------
Dt = option.Dt;
NT = floor((option.terminalT-option.startingT)/Dt);
if NT <= 1
    NT = 2;
end
% %-----------------------------------------

% %-------------- plan 2 ------------------
% h = sum(StokesmeshInfo.hElem)/StokesmeshInfo.Nelems;
% if option.theta == 1/2
%     Dt = h^((degree_u+1)/2);
% else
%     Dt = h^(degree_u+1);
% end
% NT = floor((option.terminalT-option.startingT)/Dt);
% %----------------------------------------

disp(['Dt = ',num2str(Dt)])
disp(['NT = ',num2str(NT)])

%% assemeble matrix and rhs
% zeroUU = sparse(dof_u1,dof_u1);
% zeroUP = sparse(dof_u1,dof_p);
zeroPP = sparse(dof_p,dof_p);

Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.funcOne;
Coeffs_func{2} = pde.funcOne;
Coeffs_func{3} = pde.funcZero;

tic %<<<<<<<<<<<<<< tic1, assemble StokesDarcy matrix

%% tensor Stokes and interface for Stokes
% tensor Stokes
[elem_u0v0,elem_p0q0,elem_uxvx,elem_uxvy,elem_uyvx,elem_uyvy,elem_pv1,elem_pv2,rhs_u10,rhs_u20,rhs_p0,~] = ...
    matElemCoeffsDuDvTensorStokes_timeSDT(Coeffs_func,pde.u10,pde.u20, pde.p0, StokesmeshInfo,Gaussformulas{1}, degree_u, degree_p);

[Inter_u1v1, Inter_u1v2, Inter_u2v1, Inter_u2v2, Inter_pv1, Inter_pv2] = ...
    matInterEdgeAverJumpTensorStokes_timeSDT(Coeffs_func, StokesmeshInfo, option, formulaGauss1D, degree_u, degree_p);

[Diri_u1v1,Diri_u1v2,Diri_u2v1,Diri_u2v2,Diri_pv1,Diri_pv2,~,~,~] = ...
    matDirichletEdgeAverJumpTensorStokes_timeSDT(Coeffs_func, pde.gD1, pde.gD2, StokesmeshInfo, option, formulaGauss1D, degree_u, degree_p);

elem_u1v1 = 2*(elem_uxvx + 0.5*elem_uyvy);
elem_u1v2 = 2*(0.5*elem_uyvx);
elem_u2v1 = 2*(0.5*elem_uxvy);
elem_u2v2 = 2*(0.5*elem_uxvx + elem_uyvy);

% the interface integration for Stokes 
[interface_u1v1, interface_u2v1, interface_u1v2, interface_u2v2] = ...
    matInterfaceOfStokes_timeSDT(Coeffs_func, StokesmeshInfo, interfacemeshInfo, formulaGauss1D, degree_u);

% the Stokes system matrix and rhs
StokesM = [mu*(elem_u1v1+Inter_u1v1+Diri_u1v1)+Kappa*interface_u1v1,      mu*(elem_u2v1+Inter_u2v1+Diri_u2v1)+Kappa*interface_u2v1,     -elem_pv1+Inter_pv1+Diri_pv1;
    mu*(elem_u1v2+Inter_u1v2+Diri_u1v2)+Kappa*interface_u1v2,     mu*(elem_u2v2+Inter_u2v2+Diri_u2v2)+Kappa*interface_u2v2,      -elem_pv2+Inter_pv2+Diri_pv2;
    -(-elem_pv1+Inter_pv1+Diri_pv1)',     -(-elem_pv2+Inter_pv2+Diri_pv2)',     zeroPP];

%% the Darcy 
[Darcy_u0v0, Darcy_DuDv, Darcy_rhs_phi0] = ...
    matElemCoeffsDtrialDtestPoisson_timeSDT(Coeffs_func,pde.phi0,DarcymeshInfo,Gaussformulas{1}, degree_phi, degree_phi);
Darcy_interE = ...
    matInterEdgeAverJumpPoisson_timeSDT(Coeffs_func, DarcymeshInfo,option,Gaussformulas{3},degree_phi,degree_phi);
[Darcy_DirichletE, ~] = ...
    matDirichletEdgeAverJumpPoisson_timeSDT(Coeffs_func,pde.gDp,DarcymeshInfo,option,Gaussformulas{3},degree_phi,degree_phi);
% % rhsNeumann = vecNeumannEdgePoisson_timeSD(Coeffs_func,pde,DarcymeshInfo,Gaussformulas{3},degree_phi);

DarcyM = Darcy_DuDv + Darcy_interE + Darcy_DirichletE;

%% the Stoke-Darcy coupled interface integration
[interface_uphi_n1, interface_uphi_n2, interface_phiu_n1, interface_phiu_n2] = ...
    matInterfaceStokesDarcy(Coeffs_func, StokesmeshInfo, DarcymeshInfo, interfacemeshInfo, formulaGauss1D, degree_u, degree_phi);


%% get the StokeDarcy system matrix and rhs
SD_M = sparse(dof_phi+2*dof_u1+dof_p, dof_phi+2*dof_u1+dof_p);
SD_A = sparse(dof_phi+2*dof_u1+dof_p, dof_phi+2*dof_u1+dof_p);

SD_M(1:dof_phi, 1:dof_phi) = Darcy_u0v0;
SD_M(dof_phi+1:dof_phi+dof_u1, dof_phi+1:dof_phi+dof_u1) = elem_u0v0;
SD_M(dof_phi+dof_u1+1:dof_phi+2*dof_u1, dof_phi+dof_u1+1:dof_phi+2*dof_u1) = elem_u0v0;
SD_A(1:dof_phi, 1:dof_phi) = DarcyM;
SD_A(1:dof_phi, dof_phi+1:dof_phi+dof_u1) = interface_uphi_n1;
SD_A(1:dof_phi, dof_phi+dof_u1+1:dof_phi+2*dof_u1) = interface_uphi_n2;
SD_A(dof_phi+1:dof_phi+dof_u1, 1:dof_phi) = interface_phiu_n1;
SD_A(dof_phi+dof_u1+1:dof_phi+2*dof_u1, 1:dof_phi) = interface_phiu_n2;
SD_A(dof_phi+1:end, dof_phi+1:end) = StokesM;

%--- Time-Stokes-Darcy, SD_Uh0 and systemMartix setting
u1h_0 = elem_u0v0\rhs_u10; 
u2h_0 = elem_u0v0\rhs_u20; 
ph_0 = elem_p0q0\rhs_p0; 
phi_0 = Darcy_u0v0\Darcy_rhs_phi0;
SD_Uh0 = [phi_0;u1h_0;u2h_0;ph_0];

SD_sysM = SD_M + option.theta*Dt*SD_A;
SD_rhs_term = SD_M + (option.theta-1)*Dt*SD_A;

% disp the condition number of the system matrix
disp(['condition number of Stokes-Darcy Matrix: ', num2str(condest(SD_sysM))])
sysInfo.AssembleTime = toc; %<<<<<<<<<<<<<< corresponding tic1
disp(['assemble Stokes-Darcy Matrix time: ',num2str(sysInfo.AssembleTime)])

%% Transport eqn setting
%--- Transport, Ch0 and systemMatrix setting
lastStep_DarcyUh = SD_Uh0(1:dof_phi);
lastStep_StokesUh = SD_Uh0(dof_phi+1:dof_phi+2*dof_u1);

% % [C_G, C_GxxPlusGyy, C_u1G0xPlusu2G0y, Transp_rhs0] = ...
% %     matElemCoeffsDtrialDtestTransport_timeSDT(...
% %     pde, Coeffs_func, DarcymeshInfo, StokesmeshInfo, ...
% % 	lastStep_DarcyUh, lastStep_StokesUh, ...
% % 	Gaussformulas{1}, degree_phi, degree_u, degree_c, degree_c);
% % 
% % [C_diff_K_noAdapInterE, C_diff_K_AdapInterE, C_conv_K_InterE] = ...
% % 	matInterEdgeAverJumpTransport_timeSDT(...
% % 	pde, DarcymeshInfo, StokesmeshInfo, interfacemeshInfo, ...
% % 	lastStep_DarcyUh, lastStep_StokesUh, ...
% % 	option, formulaGauss1D, degree_phi, degree_u, degree_c, degree_c);
% % 
% % [diff_DiriK, conv_DiriK] = ...
% %     matDirichletEdgeAverJumpTransport_timeSDT(...
% %     pde, DarcymeshInfo, StokesmeshInfo, ...
% %     lastStep_DarcyUh, lastStep_StokesUh, option, formulaGauss1D, ...
% %     degree_phi, degree_u, degree_c, degree_c);
% % 
% % Ch0 = C_G\Transp_rhs0;
% % 
% % C_A = C_GxxPlusGyy + C_diff_K_noAdapInterE + C_diff_K_AdapInterE + diff_DiriK ...
% %     - C_u1G0xPlusu2G0y + C_conv_K_InterE + conv_DiriK;
% % 
% % C_sysM = C_G + option.theta*Dt*C_A;
% % C_rhs_term = C_G + (option.theta-1)*Dt*C_A;


%% solve the system
tic; %<<<<<<<<<<<<<< tic2
%% time loop 
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
    
    %% Time-Stokes-Darcy setting
    [S_rhs_f1, S_rhs_f2] = vecElemRhsFtStokes_timeSD(t, pde.f1, pde.f2, StokesmeshInfo, Gaussformulas{1}, degree_u);
    D_rhs_fp = vecElemRhsFtDarcy_timeSD(t, pde.fp, DarcymeshInfo, Gaussformulas{1}, degree_phi);
    SD_vecRhs_f = [D_rhs_fp;S_rhs_f1;S_rhs_f2;zeros(dof_p,1)];
    
    if strcmpi(bd_case,'allDirichlet')
        [S_rhs_gD1, S_rhs_gD2, S_rhs_div] = vecEdgeRhsDirichletStokes_timeSD(t, Coeffs_func, pde.gD1, pde.gD2, ...
            StokesmeshInfo, option, Gaussformulas{3}, degree_u, degree_p);
        D_rhs_gD1 = vecEdgeRhsDirichletDarcy_timeSD(t, Coeffs_func, pde.gDp, ...
            DarcymeshInfo, option, Gaussformulas{3}, degree_phi);
        SD_rhs_DiriBC = [D_rhs_gD1;mu*S_rhs_gD1;mu*S_rhs_gD2;-S_rhs_div];
        SD_sysRhs = SD_rhs_term*SD_Uh0 + Dt*(SD_vecRhs_f + SD_rhs_DiriBC);
    else
        error('error in dgTimeStokesDarcySolve_usingTensorStokes.m')
    end % if
    
    SD_Uh = SD_sysM\SD_sysRhs; 
    
    % reassign the r0
    SD_Uh0 = SD_Uh;
    
    %% Transport setting
    [C_G, C_GxxPlusGyy, C_u1G0xPlusu2G0y, Transp_rhs0] = ...
        matElemCoeffsDtrialDtestTransport_timeSDT(...
        pde, Coeffs_func, DarcymeshInfo, StokesmeshInfo, ...
        lastStep_DarcyUh, lastStep_StokesUh, ...
        Gaussformulas{1}, degree_phi, degree_u, degree_c, degree_c);
    [C_diff_K_noAdapInterE, C_diff_K_AdapInterE, C_conv_K_InterE] = ...
        matInterEdgeAverJumpTransport_timeSDT(...
        pde, DarcymeshInfo, StokesmeshInfo, interfacemeshInfo, ...
        lastStep_DarcyUh, lastStep_StokesUh, ...
        option, formulaGauss1D, degree_phi, degree_u, degree_c, degree_c);
    [diff_DiriK, conv_DiriK] = ...
        matDirichletEdgeAverJumpTransport_timeSDT(...
        pde, DarcymeshInfo, StokesmeshInfo, ...
        lastStep_DarcyUh, lastStep_StokesUh, option, formulaGauss1D, ...
        degree_phi, degree_u, degree_c, degree_c);
    C_A = C_GxxPlusGyy + C_diff_K_noAdapInterE + C_diff_K_AdapInterE + diff_DiriK ...
        - C_u1G0xPlusu2G0y + C_conv_K_InterE + conv_DiriK;
    C_sysM = C_G + option.theta*Dt*C_A;
    C_rhs_term = C_G + (option.theta-1)*Dt*C_A;

    C_rhs_fc = vecElemRhsFtSDTransport1_Transport(t, ...
        pde, lastStep_DarcyUh, lastStep_StokesUh,...
        DarcymeshInfo, StokesmeshInfo, ...
        formulaGauss2D, degree_phi, degree_u, degree_c);
    C_rhs_gD = ...
        vecRhsDirichletEdgeSDTransport1_Transport( t, ...
        Coeffs_func, pde.c, DarcymeshInfo, StokesmeshInfo, ...
        option, formulaGauss1D, degree_c);
    if n==0, Ch0 = C_G\Transp_rhs0; end
    C_sysRhs = C_rhs_term*Ch0 + Dt*(C_rhs_fc + C_rhs_gD);
    Ch = C_sysM\C_sysRhs;
    
    % reassign the r0 and LastSetpUh
    Ch0 = Ch;
    
    lastStep_DarcyUh = SD_Uh0(1:dof_phi);
    lastStep_StokesUh = SD_Uh0(dof_phi+1:dof_phi+2*dof_u1);
end % for n

sysInfo.SolveSystemTime = toc; %<<<<<<<<<<<<<< corresponding tic2
disp(['solve the Time-Stokes-Darcy-Transport system: ',num2str(sysInfo.SolveSystemTime)])

%%
sysInfo.dof_phi = dof_phi;
sysInfo.dof_u1 = dof_u1;
sysInfo.dof_p = dof_p;
sysInfo.dof_Sc = dof_Sc;
sysInfo.dof_Dc = dof_Dc;
sysInfo.Gaussformulas = Gaussformulas;
sysInfo.terminalT = option.startingT + NT*Dt;
disp(['terminal Time in program:',num2str(sysInfo.terminalT)])

end % function


%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function [StokesmeshInfo, DarcymeshInfo, interfacemeshInfo] = ...
    getStokesDarcyBoundaryInfo(StokesmeshInfo, DarcymeshInfo, bd_case, pde)
%
%
%
%
%
%   YcZhang 28/5/2017
%
%   Last modified 28/5/2017
%

if strcmpi(pde.pdecase,'case1') && strcmpi(bd_case, 'allDirichlet')
    %> in pde case1, the interface is setted by y==0. 
    %> Stokes domain: [0,1]x[0,1]; Darcy domain: [0,1]x[-1,0]
    
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
    
elseif strcmpi(pde.pdecase,'case2') && strcmpi(bd_case, 'allDirichlet')
    %> in pde case2, the interface is setted by y==1.
    %> Stokes domain: [0,1]x[1,2]; Darcy domain: [0,1]x[0,1]
    
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
	%> Stokes domain: [0,1]x[1/2,1]; Darcy domain: [0,1]x[0,1/2]
    
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
else
    error('there is no case in getStokesDarcyBoundaryInfo ');
    
end % if 


end % function 


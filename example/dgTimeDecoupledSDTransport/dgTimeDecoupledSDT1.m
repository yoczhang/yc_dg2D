%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       dgTimeDecoupledSDT1, 
%       for the transport eqn, have NO inflow and outflow boundaryEdges,
%       i.e. for the transport eqn, all the boundaryEdges are Dirichlet Edges.
%   %---------------------------------------------
%
%
%   RATE OF CONVERGENCE OF DG METHOD FOR THE StokesDarcyTransport EQNs IN 2D
%
%
%
%	YcZhang 18/11/2017
%
%   Last modified 18/11/2017
%


clc
clearvars; close all;
%% Setting
% program setting
option.L0 = 1; % to uniformrefine the mesh and generate an initial mesh 
option.maxIt = 4; % to define the maximum circulation, to generate the Rate.
option.printlevel = 1;
option.plotflag = 1;

% bases element setting
option.basesType_Fu = 'P1';
option.basesType_Fp = 'P0';
option.basesType_Pp = 'P1';
option.basesType_c = 'P1';
    %> here, F standsfor fluid region, P standsfor porous region, and T standsfor Transport eq.
    
%% Time setting
% choose the time scheme
option.theta = 1; 
    %> theta = 0, forward Euler (explicit Euler),
    %> theta = 1, backward Euler (implicit Euler),
    %> theta = 1/2, C-N.
    
% set the terminal time
option.verifyCovergence = 'space'; 
    % 'space' stands for we will vertify SPACE convergence order, so we will take Dt=1/1000,
    % 'time' stands for we will veritfy TIME convergence order, so we will take very small h, and large Dt
if strcmp(option.verifyCovergence, 'space')
    option.Dt = 1/10000;
    option.startingT = 0;
    option.terminalT = 0.01;
end


%% penalty pars setting
option.p_epsilon = -1; % -1, SIPG; 1, NIPG; 0, IIPG; This is for the StokesDarcy eqn.
option.p_beta = 1;

option.p_epsilon_c = -1; % this is for the Transport eqn.
option.p_beta_c = 1;

%> for the NIPG
if option.p_epsilon == 1
    option.p_sigma = 1;
end 

if option.p_epsilon_c == 1
    option.p_sigma_c = 19;
end 

%> for the SIPG or IIPG
%--- for the Stokes-Darcy
if option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 1
    option.p_sigma = 19;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 2
    option.p_sigma = 18;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 3
    option.p_sigma = 36;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 4
    option.p_sigma = 78;
end 

%--- for the Transport eqn
if option.p_epsilon_c ~= 1 && basesType2degreek(option.basesType_c) == 1
    option.p_sigma_c = 19;
elseif option.p_epsilon_c ~= 1 && basesType2degreek(option.basesType_c) == 2
    option.p_sigma_c = 96;
elseif option.p_epsilon_c ~= 1 && basesType2degreek(option.basesType_c) == 3
    option.p_sigma_c = 36;
elseif option.p_epsilon_c ~= 1 && basesType2degreek(option.basesType_c) == 4
    option.p_sigma_c = 78;
end 


%--- for the Transport Adap Diffusivity Edges
option.p_epsilon1 = option.p_epsilon_c;
option.p_sigma1 = option.p_sigma_c;
 
%% pde setting
pde = TimeDecoupledSDTData1(0,1); % input: (coeff_case,u_case)

option.Darcy_p_sigma = pde.K; % this is to control the penalty coefficient of Daryc eqn.

%% equation
dgTimeDecoupledSDTEqn1(pde,option);

%-- end function


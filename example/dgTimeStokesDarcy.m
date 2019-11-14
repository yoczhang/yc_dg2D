%
%   RATE OF CONVERGENCE OF DG METHOD FOR THE Non-stationary StokesDarcy EQNs IN 2D
%
%
%
%	YcZhang 7/11/2017
%
%   Last modified 7/11/2017
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
option.basesType_Fu = 'P2';
option.basesType_Fp = 'P1';
option.basesType_Pp = 'P2';
    %> here, F standsfor fluid region, and P standsfor porous region.

%% penalty pars setting
option.p_epsilon = -1; % -1, SIPG; 1, NIPG; 0, IIPG;
option.p_beta = 1;

%> for the NIPG
if option.p_epsilon == 1
    option.p_sigma = 9;
end 

%> for the SIPG or IIPG
if option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 1
    option.p_sigma = 139;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 2
    option.p_sigma = 136;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 3
    option.p_sigma = 36;
    elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 4
    option.p_sigma = 78;
end 

%% setting the Tensor or Saclar Stokes
option.usingTensorStokes = 1;
 
%% Time setting
%--- choose the time scheme
option.theta = 1; 
    %> theta = 0, forward Euler (explicit Euler),
    %> theta = 1, backward Euler (implicit Euler),
    %> theta = 1/2, C-N.
    
%--- set the terminal time
option.verifyCovergence = 'space'; 
    % 'space' stands for we will vertify SPACE convergence order, so we will take Dt=1/1000,
    % 'time' stands for we will veritfy TIME convergence order, so we will take very small h, and large Dt
if strcmp(option.verifyCovergence, 'space')
    option.Dt = 1/1000;
    option.startingT = 0;
    option.terminalT = 0.02;
end

%% pde setting
pde = timeStokesDarcyData(1);

option.Darcy_p_sigma = pde.K; % this is to control the penalty coefficient of Daryc eqn.
%option.Darcy_p_sigma = 1; % this is to control the penalty coefficient of Daryc eqn.

%% equation
dgTimeStokesDarcyEqn(pde,option);

%-- end function


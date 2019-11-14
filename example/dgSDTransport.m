%
%   RATE OF CONVERGENCE OF DG METHOD FOR THE StokesDarcyTransport EQNs IN 2D
%
%
%
%	YcZhang 10/9/2017
%
%   Last modified 10/9/2017
%


clc
clearvars; close all;
%% Setting
% program setting
option.L0 = 1; % to uniformrefine the mesh and generate an initial mesh 
option.maxIt = 1; % to define the maximum circulation, to generate the Rate.
option.printlevel = 1;
option.plotflag = 1;

% bases element setting
option.basesType_Fu = 'P2';
option.basesType_Fp = 'P1';
option.basesType_Pp = 'P2';
option.basesType_c = 'P2';
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
    option.Dt = 1/1000;
    option.startingT = 0;
    option.terminalT = 0.5;
end


%% penalty pars setting
option.p_epsilon = -1; % -1, SIPG; 1, NIPG; 0, IIPG;

option.p_beta = 1;

%> for the NIPG
if option.p_epsilon == 1
    option.p_sigma = 1;
end 

%> for the SIPG or IIPG
if option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 1
    option.p_sigma = 10;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 2
    option.p_sigma = 18;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 3
    option.p_sigma = 36;
    elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 4
    option.p_sigma = 78;
end 

%> for the Adap Diffusivity Edges
option.p_epsilon1 = option.p_epsilon;
option.p_sigma1 = option.p_sigma;
 
%% pde setting
pde = SDTransportData(0,2); % input: (coeff_case,u_case)

%% equation
dgSDTransportEqn(pde,option);

%-- end function


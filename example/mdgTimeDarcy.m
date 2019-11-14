%
%   This mdgTimeDarcy.m is acorrding the paper: 
%       Mixed Discontinuous Galerkin for the Time-dependent Darcy's Problem .
%
%
%
%	YcZhang 08/27/2018 -- MM/DD/YYYY
%
%   Last modified 08/27/2018 -- MM/DD/YYYY
%


clc
clearvars; close all;
%% Setting
% program setting
option.L0 = 1; % to uniformrefine the mesh and generate an initial mesh 
option.maxIt = 3; % to define the maximum circulation, to generate the Rate.
option.printlevel = 1;
option.plotflag = 1;

% bases element setting
option.basesType_u = 'P1';
option.basesType_p = 'P0';

%% penalty pars setting
option.p_epsilon = -1; % -1, SIPG; 1, NIPG; 0, IIPG;
option.p_beta = 1;

%> for the NIPG
if option.p_epsilon == 1
    option.p_sigma = 1;
end 

%> for the SIPG or IIPG
if option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 1
    option.p_sigma = 6;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 2
    option.p_sigma = 9;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 3
    option.p_sigma = 36;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 4
    option.p_sigma = 78;
end 

%% Time setting
%--- choose the time scheme
option.theta = 1/1; 
    %> theta = 0, forward Euler (explicit Euler),
    %> theta = 1, backward Euler (implicit Euler),
    %> theta = 1/2, C-N.
    
%--- set the terminal time
option.startingT = 0;
option.terminalT = 1; 
option.TimeStepAccordingSpaceStep = true; 
    % 'true' stands for time step according the space step, else 'false',
    % we will choose Dt = 1/10000.

%% pde setting
pde = mdgTimeDarcyData(5); % input: (coeff_case, u_case)

%% equation
mdgTimeDarcyEqn(pde,option);

%-- end function


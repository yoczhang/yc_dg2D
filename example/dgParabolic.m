%
%   RATE OF CONVERGENCE OF DG METHOD FOR THE PARABOLIC EQNs IN 2D
%
%   This example is to show the rate of convergence of IPDG
%	approximation of the Parabolic equation on the unit square:
%
%       u_t - div(K \nabla u) = f in (0,T) x \Omega,
%                                  u = g_D on (0,T) x \partial\Omega,
%                                  u = u_0 on {0} x \Omega.
%
%
%
%
%	YcZhang 3/9/2017
%
%   Last modified 3/9/2017
%


clc
clearvars; close all;

%% program setting
option.L0 = 1; % to uniformrefine the mesh and generate an initial mesh 
option.maxIt = 4; % to define the maximum circulation, to generate the Rate.
option.printlevel = 1;
option.plotflag = 1;

%% DG setting
% bases element setting
option.basesType_trial = 'P1';
option.basesType_test = 'P1';

%% penalty pars setting
option.p_epsilon = -1; % -1, SIPG; 1, NIPG; 0, IIPG;
option.p_beta = 1;

%> for the NIPG
if option.p_epsilon == 1
    option.p_sigma = 1;
end 

%> for the SIPG or IIPG
if option.p_epsilon ~= 1 && basesType2degreek(option.basesType_trial) == 1
    option.p_sigma = 6;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_trial) == 2
    option.p_sigma = 18;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_trial) == 3
    option.p_sigma = 36;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_trial) == 4
    option.p_sigma = 78;
end 

%% Time setting
% choose the time scheme
option.theta = 1/2; 
    %> theta = 0, forward Euler (explicit Euler),
    %> theta = 1, backward Euler (implicit Euler),
    %> theta = 1/2, C-N.
    
% set the terminal time
option.startingT = 0;
option.terminalT = 0.5; 
option.TimeStepAccordingSpaceStep = false; 
    % 'true' stands for time step according the space step, else 'false',
    % we will choose Dt = 1/10000.


%% pde setting
pde = parabolicData(0,4);
    % input: (coeff_case, u_case)
    
dgParabolicEqn(pde,option);
    
    
%%  
% end Parabolic problem function
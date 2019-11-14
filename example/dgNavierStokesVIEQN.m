%
%   RATE OF CONVERGENCE OF DG METHOD FOR THE 
%   Navier-Stokes Variational Inequality IN 2D
%
%
%
%
%
%	YcZhang 24/10/2017
%
%   Last modified 24/10/2017
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
option.basesType_u = 'P2';
option.basesType_p = 'P1';

%% penalty pars setting
option.p_epsilon = 1; % -1, SIPG; 1, NIPG; 0, IIPG;
option.p_beta = 1;

%> for the NIPG
if option.p_epsilon == 1
    option.p_sigma = 496;
end

%> for the SIPG or IIPG
if option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 1
    option.p_sigma = 20;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 2
    option.p_sigma = 469;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 3
    option.p_sigma = 56;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 4
    option.p_sigma = 98;
end 

 
%% pde setting
pde = navierstokesVIEQNData(0,1); 
    %> input: (coeff_case, u_case)


%% equation
% 1. Navier-Stokes, to get the convergence order
dgNavierStokesVIEQNEqn(pde,option);


%-- end function


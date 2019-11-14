%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       Using the Saclar-Stokes, 
%       the dgNavierStokes2 using the Tensor-Stokes.
%   %---------------------------------------------
%
%   RATE OF CONVERGENCE OF DG METHOD FOR THE Navier-Stokes EQNs IN 2D
%
%   This example is to show the rate of convergence of DG
%	approximation of the Navier-Stokes equation.
%
%
%	YcZhang 17/10/2017
%
%   Last modified 17/10/2017
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
    option.p_sigma = 18;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 3
    option.p_sigma = 36;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 4
    option.p_sigma = 78;
end 

 
%% pde setting
pde = navierstokesData1(0,1); 
    %> input: (coeff_case, u_case)
    %> Be careful, the u_case 4 is the NavierStokesStepChannel example,
    %> and u_case 2 is not accuracy.


%% equation
% 1. Navier-Stokes, to get the convergence order
dgNavierStokesEqn1(pde,option);

% 2. NS step channel problem
% dgNavierStokesStepChannelEqn1(pde,option);

%-- end function


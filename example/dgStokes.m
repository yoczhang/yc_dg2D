%
%   RATE OF CONVERGENCE OF DG METHOD FOR THE Stokes EQNs IN 2D
%
%   This example is to show the rate of convergence of HDG
%	approximation of the Diffusion equation on the unit square:
%
%       -div(*grad u) = f in \Omega,
%                        u = gD  on \Gamma_D.
%     grad u \cdot n = (gN1,gN2)\cdot n on \Gamma_N.
%
%
%
%
%	YcZhang 12/8/2017
%
%   Last modified 12/8/2017
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
    option.p_sigma = 96;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 2
    option.p_sigma = 18;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 3
    option.p_sigma = 36;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 4
    option.p_sigma = 78;
end 

 
%% pde setting
pde = stokesData(0,2); % input: (coeff_case, u_case)

%% equation
dgStokesEqn(pde,option);

%-- end function


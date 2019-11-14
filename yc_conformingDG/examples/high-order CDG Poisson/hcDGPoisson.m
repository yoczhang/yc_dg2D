%
%   RATE OF CONVERGENCE OF high-order conforming DG METHOD FOR THE DIFFUSION EQNs IN 2D
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
%	YcZhang Jul.08.2019
%
%   Last modified Jul.08.2019
%


clc
clearvars; close all;
%% Setting
% program setting
option.L0 = 1; % to uniformrefine the mesh and generate an initial mesh 
option.maxIt = 6; % to define the maximum circulation, to generate the Rate.

% bases element setting
option.basesType = 'P1';

 
%% pde setting
pde = hcDGPoissonData(0,5);
    % input: (coeff_case, u_case)


%% equation
hcDGPoissonEqn(pde,option);

%%-- end function
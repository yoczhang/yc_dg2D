%
%   RATE OF CONVERGENCE OF DG METHOD FOR THE Varying Diffusivity DARCY EQNs IN 2D
%
%   This example is to show the rate of convergence of IPDG
%	approximation of the Transport equation on the unit square:
%
%   This example (according to LiRui's adivce), is just be designed to check the Poisson eqns. And the
%   results show that there is no overshoot and undershoot in this example.
%
%   In this example, we using the dg method for Darcy(i.e Poisson) eqns, and 
%   for varying and vanishing diffusivity, we use the adaptive-dg.
%   The methods are according to:
%   1. B.Riviere, Discontinuous Galerkin Methods for Solving Elliptic and
%       Parabolic Equations. --- 4.3.3 An improved DG method.
%   2. J.Proft and B.Riviere, Discontinuous Galerkin methods for
%       convection-diffusion equations for varying and vanishing diffusivity.
%
%
%
%	YcZhang 2/10/2017
%
%   Last modified 2/10/2017
%


clc
clearvars; close all;

%% program setting
option.L0 = 1; % to uniformrefine the mesh and generate an initial mesh 
option.maxIt = 1; % to define the maximum circulation, to generate the Rate.
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

%> for the Adap Diffusivity Edges
option.p_epsilon1 = option.p_epsilon;
option.p_sigma1 = option.p_sigma;

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
    option.terminalT = 3;
elseif strcmp(option.verifyCovergence, 'time')
    option.Dt = 1/10;
    option.startingT = 0;
    option.terminalT = 0.5;
end

%% pde setting
pde = varyingDiffusivityDarcyData(0,3);
    % input: (coeff_case, u_case)
    
dgVaryingDiffusivityDarcyEqn(pde,option);
    
    
%%  
% end Parabolic problem function
%
%   RATE OF CONVERGENCE OF DG METHOD FOR THE StokesDarcy EQNs IN 2D
%
%
%
%	YcZhang 27/8/2017
%
%   Last modified 27/8/2017
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
option.p_epsilon = 1; % -1, SIPG; 1, NIPG; 0, IIPG;
option.p_beta = 1;

%> for the NIPG
if option.p_epsilon == 1
    option.p_sigma = 0;
end 

%> for the SIPG or IIPG
if option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 1
    option.p_sigma = 19;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 2
    option.p_sigma = 196;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 3
    option.p_sigma = 36;
    elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_Fu) == 4
    option.p_sigma = 78;
end 

%% setting the Tensor or Saclar Stokes
option.usingTensorStokes = 0;
 
%% pde setting
if option.usingTensorStokes == 1
    pde = StokesDarcyData(1); % input: (u_case)
    disp('StokesDarcy using tensor Stokes eqn')
else
    pde = StokesDarcyDataUsingScalarStokes(8); % input: (u_case)
    disp('StokesDarcy using scalar Stokes eqn')
end

option.Darcy_p_sigma = pde.K; % this is to control the penalty coefficient of Daryc eqn.
% option.Darcy_p_sigma = 1; % this is to control the penalty coefficient of Daryc eqn.

%% equation
dgStokesDarcyEqn(pde,option);

%-- end function


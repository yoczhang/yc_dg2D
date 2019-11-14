%
%   RATE OF CONVERGENCE OF DG METHOD FOR THE TensorStokes EQNs IN 2D
%
%
%
%	YcZhang 25/8/2017
%
%   Last modified 26/8/2017
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
option.basesType_u = 'P4';
option.basesType_p = 'P3';

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
    option.p_sigma = 19;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 3
    option.p_sigma = 36;
elseif option.p_epsilon ~= 1 && basesType2degreek(option.basesType_u) == 4
    option.p_sigma = 78;
end 

 
%% pde setting
pde = tensorStokesData(0,7); % input: (coeff_case, u_case)

%% equation
dgTensorStokesEqn(pde,option);

%-- end function


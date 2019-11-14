% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function thetaC = dirMaximumCircumferentialStress(K1,K2,omega)
% This function will calculate the direction of crack growth according to
% the maximum circumferential stress criterion.

thetaC = 2*atan(1/4*((K1/K2-sign(K2)*sqrt((K1/K2)^2+8))))+omega;            % Angle of future crack growth (local coordinates)
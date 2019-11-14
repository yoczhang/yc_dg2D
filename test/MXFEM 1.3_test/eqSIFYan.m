% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function dK = eqSIFYan(K1,K2,thetaC)
% This function will calculate the equivalent stress intensity factor as
% defined by Yan.

dK = 1/2*cos(thetaC/2)*(K1*(1+cos(thetaC))-3*K2*sin(thetaC));               % Effective deltaK by Yan    (1992)
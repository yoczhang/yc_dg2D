% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function thetaC = dirCriticalPlane(K1,K2,s,omega)
% This function will calculate the direction of crack growth according to
% the critical plane method.

beta   = 1/2*atan(2*K2/K1);

if     s <= 1
    gamma = 1/2*acos( (-2+sqrt(4-4*(1/s^2-3)*(5-1/s^2-4*s^2)))/(2*(s-1/s^2-4*s^2)) );
elseif s >  1
    gamma = 0;
end

thetaC = beta+gamma+omega;                                                  % Angle of future crack growth (local coordinates)
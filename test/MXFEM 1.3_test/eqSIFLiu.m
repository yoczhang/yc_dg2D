% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function dK = eqSIFLiu(K1,K2,s,thetaC)
% This function will calculate the equivalent stress intensity factor as
% defined by Liu.

k1 = K1/2*(1+cos(2*thetaC))+K2*sin(2*thetaC);
k2 = K1/2*sin(2*thetaC)+K2*cos(2*thetaC);
kH = K1/3;

if     s <= 1
    gamma = 1/2*arccos( (-2+sqrt(4-4*(1/s^2-3)*(5-1/s^2-4*s^2)))/(2*(s-1/s^2-4*s^2)) );
    A     = 0;
    B     = sqrt((s*cos(2*gamma))^2+(sin(2*gamma))^2);
elseif s >  1
    A     = 9*(s^2-1);
    B     = s;
end

dK = 1/B*sqrt(k1^2+(k2/s)^2+A*kH^2);                                        % Effective deltaK by Liu (1999)
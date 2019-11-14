% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function dK = eqSIFRhee(K1,K2)
% This function will calculate the equivalent stress intensity factor as
% defined by Rhee.

dK = sqrt(K1^2+K2^2);                                                       % Effective deltaK by Rhee   (1987)
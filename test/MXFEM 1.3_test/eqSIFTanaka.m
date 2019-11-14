% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function dK = eqSIFTanaka(K1,K2)
% This function will calculate the equivalent stress intensity factor as
% defined by Tanaka.

dK = (K1^4+8*K2^4)^(1/4);                                                   % Effective deltaK by Tanaka (1974)
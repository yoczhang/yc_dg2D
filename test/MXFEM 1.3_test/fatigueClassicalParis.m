% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function da = fatigueClassicalParis(dN,C,m,dK)
% This function will calculate the crack growth increment according to the
% classical Paris model. Assumed units are Pascals.

da = C*dN*(dK/1E6)^m;                                                       % Crack growth increment from Paris model
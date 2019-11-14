% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [da,Kminp,Kol,aOL] = fatigueModifiedParis(dN,C,m,Kmax,Kmin,Kminp,Kol,aOL,R,b,b1,n)
% This function will calculate the crack growth increment according to the
% classical Paris model. Assumed units are Pascals.

global MAT 

if     (-5.0 <= R) && (R < 0.0)
    Mr = (1-R)^(-b1);
elseif  (0.0 <= R) && (R < 0.5)
    Mr = (1-R)^(-b);
elseif  (0.5 <= R) && (R <= 1.0)
    Mr = (1.05-1.4*R+0.6*R^2)^(-b);
end

crackCoord2Length;
v     = MAT(2);
Sy    = MAT(10);
if MAT(5) == 1
    alpha = pi/8;
elseif MAT(5) == 2
    alpha = (1-2*v)^2*pi/8;
end
ry    = alpha*(Kmax/Sy)^2;
rOL   = alpha*(Kol/Sy)^2;
Ku    = Kminp-Kmin;
rd    = alpha*(Ku/Sy)^2;
a     = a(end);

if     (a+ry) <  (aOL+rOL-rd)
    Mp = (ry/(aOL+rOL-a-rd))^n;
elseif (a+ry) >= (aOL+rOL-rd)
    Mp = 1;
end

K1t   = MAT(8);
da    = C*dN*( (Mr*Mp*(Kmax-Kmin)/1E6)^m - (K1t/1E6)^m );                    % Crack growth increment from Paris model
Kminp = Kmin;
if ry > rOL
    Kol   = Kmax;
    aOL   = a;
end
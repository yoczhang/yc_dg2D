% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [uu, vv, Eps11, Eps12, Eps22, Sigma11, Sigma12, Sigma22] = ...
    ExactSol_Hole(xx, yy, aa, kappa, lambda, mu)

% Get the exact displacements, strains and stresses according to the 
% structure-with-hole test case.

% Get polar coordinates (rr, th) out of (xx, yy).
rr = sqrt(xx*xx+yy*yy);
if xx == 0 
    if yy == 0
        disp('Theta is undetermined for (x,y) = (0,0)!')
        th = 0;
    elseif yy > 0
        th = 0.5*pi;
    else
        th = 1.5*pi;
    end
elseif yy == 0
    if xx > 0
        th = 0;
    else
        th = pi;
    end
elseif xx > 0
    th = atan(yy/xx);
elseif xx < 0
    th = pi+atan(yy/xx);
else
    error('Internal error!')
end

% Get exact solution.
uu = aa/(8*mu) * (rr/aa*(kappa+1)*cos(th) + 2*aa/rr*((1+kappa)*cos(th)+cos(3*th)) - ...
    2*(aa*aa*aa)/(rr*rr*rr)*cos(3*th));
vv = aa/(8*mu) * (rr/aa*(kappa-3)*sin(th) + 2*aa/rr*((1-kappa)*sin(th)+sin(3*th)) - ...
    2*(aa*aa*aa)/(rr*rr*rr)*sin(3*th));

Sigma11 = 1 - (aa*aa)/(rr*rr)*(3/2*cos(2*th)+cos(4*th)) + 3/2*(aa*aa*aa*aa)/(rr*rr*rr*rr)*cos(4*th);
Sigma22 =   - (aa*aa)/(rr*rr)*(1/2*cos(2*th)-cos(4*th)) - 3/2*(aa*aa*aa*aa)/(rr*rr*rr*rr)*cos(4*th);
Sigma12 =   - (aa*aa)/(rr*rr)*(1/2*sin(2*th)+sin(4*th)) + 3/2*(aa*aa*aa*aa)/(rr*rr*rr*rr)*sin(4*th);

Sigma = [Sigma11 Sigma12; Sigma12 Sigma22];
Eps = -0.25*lambda/(mu*(lambda+mu))*trace(Sigma)*[1 0; 0 1] + 1/(2*mu)*Sigma;
Eps11 = Eps(1,1);
Eps12 = Eps(1,2);
Eps22 = Eps(2,2);

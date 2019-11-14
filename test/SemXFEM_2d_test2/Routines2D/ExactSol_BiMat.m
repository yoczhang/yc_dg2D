% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [uX, uY, Eps11, Eps12, Eps22, Sigma11, Sigma12, Sigma22] = ...
    ExactSol_BiMat(xx, yy, aa, bb, delta, mu1, mu2, lambda1, lambda2)

% Get the exact displacements, strains and stresses according to the 
% bi-material-structure test case.

rr = sqrt(xx*xx + yy*yy);
if rr <= aa
    uR = ( (1-bb*bb/(aa*aa))*delta+bb*bb/(aa*aa) ) * rr; 
    eps_RR = (1-bb*bb/(aa*aa))*delta+bb*bb/(aa*aa);
    eps_TT = (1-bb*bb/(aa*aa))*delta+bb*bb/(aa*aa);
    sig_RR = 2*mu1*eps_RR+lambda1*(eps_RR+eps_TT);
    sig_TT = 2*mu1*eps_TT+lambda1*(eps_RR+eps_TT);
else
    uR = (rr-bb*bb/rr)*delta+bb*bb/rr;
    eps_RR = (1+bb*bb/(rr*rr))*delta-bb*bb/(rr*rr);
    eps_TT = (1-bb*bb/(rr*rr))*delta+bb*bb/(rr*rr);
    sig_RR = 2*mu2*eps_RR+lambda2*(eps_RR+eps_TT);
    sig_TT = 2*mu2*eps_TT+lambda2*(eps_RR+eps_TT);
end

uT = 0;
EPS_RT = [eps_RR 0; 0 eps_TT];
SIG_RT = [sig_RR 0; 0 sig_TT];

T  = [xx/rr  yy/rr; -yy/rr xx/rr]; % = [cos(th)  sin(th); -sin(th) cos(th)]
Tt = [xx/rr -yy/rr;  yy/rr xx/rr]; % = [cos(th) -sin(th);  sin(th) cos(th)]
uX = Tt(1, :) * [uR; uT];
uY = Tt(2, :) * [uR; uT];

% Make transformation of stress tensor from polar to (x,y)-coordinates.
EPS_XY = Tt * EPS_RT * T;
SIG_XY = Tt * SIG_RT * T;

Eps11   = EPS_XY(1, 1);
Eps22   = EPS_XY(2, 2);
Eps12   = EPS_XY(1, 2);
Sigma11 = SIG_XY(1, 1);
Sigma22 = SIG_XY(2, 2);
Sigma12 = SIG_XY(1, 2);

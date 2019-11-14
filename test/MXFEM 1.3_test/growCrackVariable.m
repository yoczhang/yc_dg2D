% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [exit,Kminp,Kol,aOL] = growCrackVariable(K1,K2,omega,i,Kminp,Kol,aOL)
% This function will grow the crack in the direction determined by the
% critical plane method according to the modified Paris model.

global CRACK MAT GROW

Kc   = MAT(7);                                                              % Critical stress intensity factor
K1t  = MAT(8);                                                              % Threshold Mode I stress intensity factor
K2t  = MAT(9);                                                              % Threshold Mode II stress intensity factor
s    = K1t/K2t;

C    = GROW(3);
m    = GROW(4);
b    = GROW(5);
b1   = GROW(6);
n    = GROW(7);

stressHistory = inputLoadHistory(ceil(i/2));
dN = stressHistory(7);

exit = 0;
nCT  = length(omega);                                                       % Number of crack tips in domain

for iCT = 1:nCT    
    
    thetaC = dirCriticalPlane(K1(i),K2(i),s,omega(iCT));                          % Crack growth direction from critical plane
%     thetaC = dirMaximumCircumferentialStress(K1(i),K2(i),omega(iCT));
    Ka = eqSIFLiu(K1(i),K2(i),s,thetaC);                                       % Effective deltaK by Liu (1999)
%     Ka = eqSIFTanaka(K1(i),K2(i));
    thetaC = dirCriticalPlane(K1(i-1),K2(i-1),s,omega(iCT));                          % Crack growth direction from critical plane    
%     thetaC = dirMaximumCircumferentialStress(K1(i),K2(i),omega(iCT));
    Kb = eqSIFLiu(K1(i-1),K2(i-1),s,thetaC);                                 % Effective deltaK by Liu (1999)
% Kb = eqSIFTanaka(K1(i-1),K2(i-1));
if abs(K1(i)) > abs(K1(i-1))
    thetaC = dirMaximumCircumferentialStress(abs(K1(i)-K1(i-1)),K2(i)-K2(i-1),omega(iCT));
%     thetaC = dirCriticalPlane(abs(K1(i)-K1(i-1)),K2(i)-K2(i-1),s,omega(iCT))                          % Crack growth direction from critical plane    
else
    thetaC = dirMaximumCircumferentialStress(abs(K1(i-1)-K1(i)),K2(i-1)-K2(i),omega(iCT));
%     thetaC = dirCriticalPlane(abs(K1(i-1)-K1(i)),K2(i-1)-K2(i),s,omega(iCT))                          % Crack growth direction from critical plane    
end
    
    dK = abs(Ka-Kb);
    
    
%     if dK < Kc
        nSeg = size(CRACK,1);
        Kmax = max(Ka,Kb);
        Kmin = min(Ka,Kb);
        R    = Kmin/Kmax;
        [da,Kminp,Kol,aOL] = fatigueModifiedParis(dN,C,m,Kmax,Kmin,Kminp,Kol,aOL,R,b,b1,n);
        
        if da > 0
%           thetaC;
            da
	else
	    da = 0;
	end

            dx = da*cos(thetaC);                                            % Increment of growth in the x-direction
            dy = da*sin(thetaC);                                            % Increment of growth in the y-direction

            % Add new crack portion to exisiting portion
            if iCT == 1
                CRACK(nSeg+1,:) = [CRACK(nSeg,1)+dx,CRACK(nSeg,2)+dy];
            elseif iCT == 2
                newCR = zeros(nSeg+1,2);
                newCR(2:nSeg+1,:) = CRACK;
                newCR(1,:) = [CRACK(1,1)+dx,CRACK(1,2)+dy];
                CRACK = newCR;
            end
        
        
%     else
%         exit = 1;
%     end
end
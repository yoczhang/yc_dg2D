% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function exit = growCrack(KI,KII,omega)
% This function will grow the crack in the direction determined by the
% maximum circumferential stress criterion by a magnitude determined by 
% the input file.

global CRACK MAT GROW


Kc   = MAT(7);                                                              % Critical stress intensity factor
K1t  = MAT(8);                                                              % Threshold Mode I stress intensity factor
K2t  = MAT(9);                                                              % Threshold Mode II stress intensity factor
s    = K1t/K2t;
exit = 'NO';
nCT  = length(omega);                                                       % Number of crack tips in domain

for iCT = 1:nCT
    K1     = KI(iCT);                                                       % Mode I stress intensity factor for ith tip
    K2     = KII(iCT);                                                      % Mode II stress intensity factor for ith tip
    
%     thetaC = dirMaximumCircumferentialStress(K1,K2,omega(iCT));             % Crack growth direction from maximum circumferential stress 
    thetaC = dirCriticalPlane(K1,K2,s,omega(iCT));                          % Crack growth direction from critical plane method
    
    dK1 = ;
    dK2 = ;
    
%     dK = eqSIFRhee(dK1,dK2);                                                  % Effective deltaK by Rhee   (1987)
%     dK = eqSIFTanaka(dK1,dK2);                                                % Effective deltaK by Tanaka (1974)
%     dK = eqSIFYan(dK1,dK2,thetaC);                                            % Effective deltaK by Yan    (1992)
    dK = eqSIFLiu(dK1,dK2,s,thetaC);                                          % Effective deltaK by Liu (1999)

    if dK < Kc
        nSeg = size(CRACK,1);
        
        if length(GROW) == 2                                                % Constant deltaA used to grow crack
            da = GROW(2);                                                   % User-defined crack growth increment
        elseif length(GROW) == 4                                            % Classical Paris model used to grow crack
            da = fatigueClassicalParis(dN,C,m,dK);                          % Crack growth increment from Paris model
        elseif length(GROW) == 7
            da = fatigueModifiedParis(dN,C,m,Kmax,Kmin,Kminp,Kol,aOL,R,b,b1,n);
        end
        
        if da > 0
            dx = da*cos(thetaC);                                            % Increment of growth in the x-direction
            dy = da*sin(thetaC);                                            % Increment of growth in the y-direction
        else
            dx = 0;
            dy = 0;
        end
        
        % Add new crack portion to exisiting portion
        if iCT == 1
            CRACK(nSeg+1,:) = [CRACK(nSeg,1)+dx,CRACK(nSeg,2)+dy];
        elseif iCT == 2
            newCR = zeros(nSeg+1,2);
            newCR(2:nSeg+1,:) = CRACK;
            newCR(1,:) = [CRACK(1,1)+dx,CRACK(1,2)+dy];
            CRACK = newCR;            
        end
    else
        exit = 'YES';
    end
end
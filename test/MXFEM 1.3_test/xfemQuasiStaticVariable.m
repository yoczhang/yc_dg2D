% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

clear all; close all; clc; format compact; tic; global CRACK NODES K1 K2

%%%% Pre-Processing %%%%%
inputQuasiStatic;                                                           % Define the geometry, materials, discontinuities

iter = numIterations; tFS = zeros(iter,1); tAS = tFS;
for i = 1:iter
    
    %%%%% Processing %%%%%
    if i == 1, connectivity; pHDOF = []; else NODES(:,4:29) = 0; end        % Define connectivity
    omega           = levelSet(i);                                          % Create phi and psi, define enriched elements
    [DOF,DISP]      = calcDOF;                                              % Total degrees of freedom
    [updElem,IElem] = enrElem(i,pHDOF);                                     % Find enriched elements, inclusion elements

    if i == 1,  globalK = stiffnessMatrix(omega,DOF,iter,updElem,IElem);    % Construct global stiffness matrix        
    else globalK = updateStiffness(globalK,omega,DOF,updElem,pHDOF); end    % Update the global stiffness matrix
    
    globalF         = forceVector(DOF,i);                                   % Construct global force vector
    freeDOF         = boundaryCond(DOF);                                    % Solve for the degrees of freedom

    DISP(freeDOF,:) = globalK(freeDOF,freeDOF)\globalF(freeDOF,:);          % Find the nodal displacement values
    
    %%% Post-Processing %%%%%
    if i == iter, plotMain(omega,DISP); end                                 % Make plots
    if isempty(CRACK) == 0
        pHDOF    = 2*max(NODES(:,2));                                       % Maximum constant DOF at current iteration 
        [KI,KII] = JIntegral(omega,DISP);                                   % Calculate the stress intensity factors
        K1(i) = KI; K2(i) = KII;
        if mod(i,2) == 0
            if i == 2, Kminp = 0; Kol = 0; aOL = 0; end
            [exit,Kminp,Kol,aOL] = growCrackVariable(K1,K2,omega,i,Kminp,Kol,aOL);   % Advance crack for quasi-static growth
            
            if exit == 1
                disp('WARNING: No crack growth, iterations exited early.')
                plotMain(omega,DISP); break
            end
        end
    end
    
    disp(['Iteration ',num2str(i),' completed. Elapsed time is ',num2str(toc,'%0.4f'),'.'])    
end
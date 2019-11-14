function L2err = dgTimeL2JumpError(t,P1,P2,Ph1,Ph2,meshInfo,formulaGauss1D,basesType_trial)
%
%   We want to compute the L2 and H1 error on ELEMENTs.
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases denote the number of LOCAL bases on each K of Th.
%
%   input:
%       P, vectorized function of two variables (x,y), i.e. the true solutions.
%       Ph, discontinuous Pk function, [NTbases*Nelems x 1]
%       meshInfo, the mesh information.
%       formulaGauss2D, the 2d Gauss quadrature formula, size: a matrix, Npoints x 3,
%               the first column is the x-coordinates of all Gauss-Points,
%               the second column is the y-coordinates of all Gauss-Points,
%               the third is the weights of all Gauss-Points.
%       basesType_trial, polynomial degree.
%
%   output:
%       L2err, \| P-Ph \|_{L^2(Th)}, a scalar.
%       H2err, \| P-Ph \|_{H^1(Th)}, a scalar.
%
%
%   YcZhang 12/8/2017
%
%   Last modified 12/8/2017
%

% mesh information: interior edges 
% Nelems = meshInfo.Nelems;
interEdges = meshInfo.interEdgeIndex; 
trial_k = basesType2degreek(basesType_trial);
Nbases_trial = (trial_k+1)*(trial_k+2)/2;

L2err = 0;
mapRefE2PhyE = @(s, P1, P2) s.*P1+(1-s).*P2; % meshInfo.mapRefE2PhyE 

for CurrEdge = 1:length(interEdges)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
        % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    interEindx = interEdges(CurrEdge); 
    ePoint1 = meshInfo.node(meshInfo.edge(interEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(interEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = meshInfo.areaEdge(interEindx);
    
    phyGpointsX = mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = areaEdge * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
    elem1 = meshInfo.edge2elem(interEindx,1);
    elem2 = meshInfo.edge2elem(interEindx,2);
    local_e1 = meshInfo.edge2elem(interEindx,3);
    local_e2 = meshInfo.edge2elem(interEindx,4);
        %> [n1, n2]=edge2elem(k,1:2) stands for the n1-th elem and n2-th elem share the k-th edge.
        %> local_e1=edge2elem(k,3) stands for the local edge index in n1-th elem of the k-th edge .
        %> also, local_e2=edge2elem(k,4) stands for the local edge index in n2-th elem of the k-th edge .
    elemP = min(elem1, elem2); % let the smaller number of element as the element +.
    elemM = max(elem1, elem2); % let the larger number of element as the element -.
    
    if elemP == elem1
        eIndx_elemP = local_e1;
        eIndx_elemM = local_e2;
    else
        eIndx_elemP = local_e2;
        eIndx_elemM = local_e1;
    end 
    
    elemP_nu = meshInfo.nuEdge0Elem{elemP}(:,eIndx_elemP); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleP.
    elemM_nu = meshInfo.nuEdge0Elem{elemM}(:,eIndx_elemM); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleM.

    elemP_xT = meshInfo.centroidElem(elemP,1);  
    elemP_yT = meshInfo.centroidElem(elemP,2);
    elemP_hT = meshInfo.hElem(elemP);
        %> the centroid(xing xin) and diameter of eleP.
    elemM_xT = meshInfo.centroidElem(elemM,1);  
    elemM_yT = meshInfo.centroidElem(elemM,2);
    elemM_hT = meshInfo.hElem(elemM);
        %> the centroid(xing xin) and diameter of eleM.
        
    Ph1OnElem_P = Ph1((elemP-1)*Nbases_trial+1 : elemP*Nbases_trial,1); % [NTg1bases_trial x 1]
    Ph1OnElem_M = Ph1((elemM-1)*Nbases_trial+1 : elemM*Nbases_trial,1); % [NTg1bases_trial x 1]
    Ph2OnElem_P = Ph2((elemP-1)*Nbases_trial+1 : elemP*Nbases_trial,1); % [NTg1bases_trial x 1]
    Ph2OnElem_M = Ph2((elemM-1)*Nbases_trial+1 : elemM*Nbases_trial,1); % [NTg1bases_trial x 1]
    
    %-- subPart -- get the GROUP1 bases on quad
    [elemP_Pb,  ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, trial_k);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
    [elemM_Pb,  ~, ~] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, trial_k);
        %> eleMtrialPb, eleMtrialPbx, eleMtrialPby, [Npoints x NTg1bases_trial]
        
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
   %--- the funcValue may be chosen by case.
    valueP1 = P1(t,phyGpointsX, phyGpointsY); % [Npoints x 1] 
    valueP2 = P2(t,phyGpointsX, phyGpointsY); % [Npoints x 1] 
    
    Err_temp = (valueP1 - elemP_Pb*Ph1OnElem_P)*elemP_nu(1) + (valueP2 - elemP_Pb*Ph2OnElem_P)*elemP_nu(2) ... % e_u1^{+}*n1^{+}+e_u2^{+}*n2^{+}
        + (valueP1 - elemM_Pb*Ph1OnElem_M)*elemM_nu(1) + (valueP2 - elemM_Pb*Ph2OnElem_M)*elemM_nu(2);% e_u1^{-}*n1^{-}+e_u2^{-}*n2^{-}
    L2err = L2err +  phyGweights'*Err_temp.^2/areaEdge;

end % for ii
L2err = sqrt(L2err);

end % function dgL2H1Error


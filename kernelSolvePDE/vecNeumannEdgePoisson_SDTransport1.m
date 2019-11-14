function rhsNeumann = vecNeumannEdgePoisson_SDTransport1(Coeffs_func,pde,meshInfo,formulaGauss1D,degreek_test)
%
%   This function is to compute the Neumann edges of the rhs
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases_trial denote the number of LOCAL trial bases on each K of Th.
%               NTbases_test denote the number of LOCAL test bases on each K of Th.
%
%   input:
%       Coeffs, the cell-type, here 
%                   Coeffs{1}, and Coeffs{2} are the K=[k11, k12; k21, k22], 
%                   Coeffs{3} is the (\int_T Coeffs{1}uv dx),
%                   let k12=0, k21=0, the Coeffs{1}=k11, Coeffs{2}= k22. 
%                   Coeffs{1} is the (\int_T Coeffs{1} u_x v dx). 
%                   Coeffs{2} is the (\int_T Coeffs{2} u_y v dx).
%       pde, we will need the analysis solution ux, uy to get the (\nabla u \cdot nu).
%       meshInfo, the mesh information.
%       formulaGauss1D, the 1d Gauss quadrature formula, size: a matrix, [Npoints x 2],
%               the first column is the all Gauss-Points on [0,1],
%               the second is the weights of all Gauss-Points.
%       degreek_trial, the trial polynomial degree k.
%       degreek_test, the test polynomial degree k.
%
%   output:
%       rhsNeumann, the value \int_Ne (\nabla u \cdot nu) v dx.
%
%
%   YcZhang 30/8/2017
%
%   Last modified 30/8/2017
%

% mesh information
Nelems = meshInfo.Nelems;
NeumannEdgeIndex = meshInfo.NeumannEdgeIndex;

% bases setting
NTbases = (degreek_test+1)*(degreek_test+2)/2;

% coefficient setting
Coeff_k11 = Coeffs_func{1};
Coeff_k22 = Coeffs_func{2};

% output vec setting
rhsNeumann = zeros(Nelems*NTbases, 1);

for ii = 1:length(NeumannEdgeIndex)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    NEindx = NeumannEdgeIndex(ii); 
    ePoint1 = meshInfo.node(meshInfo.edge(NEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(NEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    
    phyGpointsX = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = meshInfo.areaEdge(NEindx) * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
        
    elem1 = meshInfo.edge2elem(NEindx,1);
    elem2 = meshInfo.edge2elem(NEindx,2);
    local_e1 = meshInfo.edge2elem(NEindx,3);
    local_e2 = meshInfo.edge2elem(NEindx,4);
        %> [n1, n2]=edge2elem(k,1:2) stands for the n1-th elem and n2-th elem share the k-th edge.
        %> local_e1=edge2elem(k,3) stands for the local edge index in n1-th elem of the k-th edge .
        %> also, local_e2=edge2elem(k,4) stands for the local edge index in n2-th elem of the k-th edge .
    elemP = max(elem1, elem2); % for the boundary edge, there will be one zero of (ele1,ele2), for let the bigger number of element as the element +.
    
    if elemP == elem1
        eIndx_elemP = local_e1;
    else
        eIndx_elemP = local_e2;
    end 
    
    elemP_nu = meshInfo.nuEdge0Elem{elemP}(:,eIndx_elemP); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleP.
    
    elemP_xT = meshInfo.baryElem(elemP,1);  
    elemP_yT = meshInfo.baryElem(elemP,2);
    elemP_hT = meshInfo.hElem(elemP);
        %> the centroid(xing xin) and diameter of eleP.
        
	%-- get the Diffusivity coefficients of elemP and elemM
    elemP_Permeability = meshInfo.PermeabilityCoeffs(elemP);
        
    %-- get the value of Coeffs func at Gauss Points
    valueCoeff_k11 = Coeff_k11(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valueCoeff_k22 = Coeff_k22(phyGpointsX, phyGpointsY); % [Npoints x 1]
    value_phix = pde.phix(phyGpointsX, phyGpointsY); % [Npoints x 1]
    value_phiy = pde.phiy(phyGpointsX, phyGpointsY); % [Npoints x 1]
    elemP_value_Permeability = elemP_Permeability * ones(length(phyGweights),1); % [Npoints x 1]
    
    %-- subPart -- get the bases on quad
    [elemP_testPb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_test);
        %> elemP_testPb, [Npoints x NTg1bases_test].

    %<<-- End Part I ---------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
    
    %% Part II, compute rhs Neumann
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    % get the location 
    rhs_Row = ( (elemP-1)*NTbases+1 : elemP*NTbases )';
    
    % the value rhs Neumann
    rhsNeu_temp = elemP_testPb' * (elemP_nu(1)*phyGweights.*(elemP_value_Permeability.*valueCoeff_k11.*value_phix) ...
        + elemP_nu(2)*phyGweights.*(elemP_value_Permeability.*valueCoeff_k22.*value_phiy));
    
    % add to the rhsNeumann
    rhsNeumann(rhs_Row,1) = rhsNeumann(rhs_Row,1) + rhsNeu_temp;
    
end % for 

end % function
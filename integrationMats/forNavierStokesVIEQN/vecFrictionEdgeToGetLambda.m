function Newlambda = vecFrictionEdgeToGetLambda(fric_g1, fric_g2, lastStepUh, ...
    meshInfo, formulaGauss1D, degreek_u)
%
%
%
%
%
%
%
%
%   YcZhang 25/10/2017
%
%   Last modified 25/10/2017
%
%

% mesh information: interior edges 
Nelems = meshInfo.Nelems;
FrictionEdges_1 = meshInfo.FrictionEdgeIndex_1; % here we need the FrictionEdgeIndex_1 information.
FrictionEdges_2 = meshInfo.FrictionEdgeIndex_2; % here we need the FrictionEdgeIndex_1 information.

NTbases_u = (degreek_u+1)*(degreek_u+2)/2;
Newlambda = zeros(size(formulaGauss1D,1),1);


for CurrEdge = 1:length(FrictionEdges_1)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    FrictionEindx = FrictionEdges_1(CurrEdge); 
    ePoint1 = meshInfo.node(meshInfo.edge(FrictionEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(FrictionEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    
    phyGpointsX = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    
    elem1 = meshInfo.edge2elem(FrictionEindx,1);
    elem2 = meshInfo.edge2elem(FrictionEindx,2);
    local_e1 = meshInfo.edge2elem(FrictionEindx,3);
    local_e2 = meshInfo.edge2elem(FrictionEindx,4);
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
    elemP_tau = [-elemP_nu(2),elemP_nu(1)];
    
    elemP_xT = meshInfo.centroidElem(elemP,1);  
    elemP_yT = meshInfo.centroidElem(elemP,2);
    elemP_hT = meshInfo.hElem(elemP);
        %> the centroid(xing xin) and diameter of eleP.
        
    %-- get the value of Coeffs func at Gauss Points
    valuefunc_fric_g1 = fric_g1(phyGpointsX, phyGpointsY); % [Npoints x 1]
    
    Elem_uh1 = lastStepUh((elemP-1)*NTbases_u+1:elemP*NTbases_u);
    Elem_uh2 = lastStepUh(Nelems*NTbases_u+(elemP-1)*NTbases_u+1:Nelems*NTbases_u+elemP*NTbases_u);
    
    %-- subPart -- get the bases on quad
    [u_elemP_Pb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_u);
        %> u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, [Npoints x NTbases_u]
    

    %<<-- End Part I ---------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
    
    %% Part II, compute the 
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %    
    %-- subPart 3 -- the rhs Dirichlet edges
    
    Newlambda = Newlambda + elemP_tau(1) * valuefunc_fric_g1.*(u_elemP_Pb * Elem_uh1) ...
        + elemP_tau(2) * valuefunc_fric_g1.*(u_elemP_Pb * Elem_uh2);
    
    %<<-- End Part II --------------------------------------------------------------------------------------
end % for FrictionEdgeIndex_1


for CurrEdge = 1:length(FrictionEdges_2)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    FrictionEindx = FrictionEdges_2(CurrEdge); 
    ePoint1 = meshInfo.node(meshInfo.edge(FrictionEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(FrictionEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    
    phyGpointsX = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    
    elem1 = meshInfo.edge2elem(FrictionEindx,1);
    elem2 = meshInfo.edge2elem(FrictionEindx,2);
    local_e1 = meshInfo.edge2elem(FrictionEindx,3);
    local_e2 = meshInfo.edge2elem(FrictionEindx,4);
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
    elemP_tau = [-elemP_nu(2),elemP_nu(1)];
    
    elemP_xT = meshInfo.centroidElem(elemP,1);  
    elemP_yT = meshInfo.centroidElem(elemP,2);
    elemP_hT = meshInfo.hElem(elemP);
        %> the centroid(xing xin) and diameter of eleP.
        
    %-- get the value of Coeffs func at Gauss Points
    valuefunc_fric_g2 = fric_g2(phyGpointsX, phyGpointsY); % [Npoints x 1]
    
    Elem_uh1 = lastStepUh((elemP-1)*NTbases_u+1:elemP*NTbases_u);
    Elem_uh2 = lastStepUh(Nelems*NTbases_u+(elemP-1)*NTbases_u+1:Nelems*NTbases_u+elemP*NTbases_u);
    
    %-- subPart -- get the bases on quad
    [u_elemP_Pb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_u);
        %> u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, [Npoints x NTbases_u]
    

    %<<-- End Part I ---------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
    
    %% Part II, compute the 
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %    
    %-- subPart 3 -- the rhs Dirichlet edges
    
    Newlambda = Newlambda + elemP_tau(1) * valuefunc_fric_g2.*(u_elemP_Pb * Elem_uh1) ...
        + elemP_tau(2) * valuefunc_fric_g2.*(u_elemP_Pb * Elem_uh2);
    
    %<<-- End Part II --------------------------------------------------------------------------------------
end % for FrictionEdgeIndex_2

end % function
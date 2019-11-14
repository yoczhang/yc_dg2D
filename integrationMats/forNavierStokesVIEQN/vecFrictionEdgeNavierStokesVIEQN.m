function [rhsFricEdge_1, rhsFricEdge_2] = ...
    vecFrictionEdgeNavierStokesVIEQN(fric_g1, fric_g2, lastStepLambda, ...
    meshInfo, formulaGauss1D, degreek_u)
%
%
%
%
%   YcZhang 24/10/2017
%
%   Last modified 24/10/2017
%
%


% mesh information: interior edges 
Nelems = meshInfo.Nelems;
FrictionEdges_1 = meshInfo.FrictionEdgeIndex_1; % here we need the FrictionEdgeIndex_1 information.
FrictionEdges_2 = meshInfo.FrictionEdgeIndex_2; % here we need the FrictionEdgeIndex_1 information.

NTbases_u = (degreek_u+1)*(degreek_u+2)/2;
rhsFricEdge_1 = zeros(Nelems*NTbases_u,1);
rhsFricEdge_2 = zeros(Nelems*NTbases_u,1);

for CurrEdge = 1:length(FrictionEdges_1)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    FrictionEindx = FrictionEdges_1(CurrEdge); 
    ePoint1 = meshInfo.node(meshInfo.edge(FrictionEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(FrictionEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = meshInfo.areaEdge(FrictionEindx);
    
    phyGpointsX = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = areaEdge * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
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
    %-- get the rhs Dirchlet vector
    [rhsRow, ~] = getROWandCOL(elemP, 0, 1, NTbases_u);
    
    r_temp_1 = elemP_tau(1) * u_elemP_Pb' *  (lastStepLambda.*valuefunc_fric_g1.*phyGweights);
        %> [NTbases_test x Npoints] * [Npoints x 1] = [NTbases_test x 1].
    r_temp_2 = elemP_tau(2) * u_elemP_Pb' *  (lastStepLambda.*valuefunc_fric_g1.*phyGweights);
        %> [NTbases_test x Npoints] * [Npoints x 1] = [NTbases_test x 1].
    rhsFricEdge_1(rhsRow,1) = rhsFricEdge_1(rhsRow,1) + r_temp_1;
    rhsFricEdge_2(rhsRow,1) = rhsFricEdge_2(rhsRow,1) + r_temp_2;
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
    areaEdge = meshInfo.areaEdge(FrictionEindx);
    
    phyGpointsX = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = areaEdge * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
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
    %-- get the rhs Dirchlet vector
    [rhsRow, ~] = getROWandCOL(elemP, 0, 1, NTbases_u);
    
    r_temp_1 = elemP_tau(1) * u_elemP_Pb' *  (lastStepLambda.*valuefunc_fric_g2.*phyGweights);
        %> [NTbases_test x Npoints] * [Npoints x 1] = [NTbases_test x 1].
    r_temp_2 = elemP_tau(2) * u_elemP_Pb' *  (lastStepLambda.*valuefunc_fric_g2.*phyGweights);
        %> [NTbases_test x Npoints] * [Npoints x 1] = [NTbases_test x 1].
    rhsFricEdge_1(rhsRow,1) = rhsFricEdge_1(rhsRow,1) + r_temp_1;
    rhsFricEdge_2(rhsRow,1) = rhsFricEdge_2(rhsRow,1) + r_temp_2;
    %<<-- End Part II --------------------------------------------------------------------------------------
end % for FrictionEdgeIndex_2


end % function



%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%% get the Row and Col
function [ROW, COL] = getROWandCOL(elemP, elemM, NTbases_trial, NTbases_test)
%
%   elemM, if elemM==0, standsfor the boundary is the boundary edge.
%

if elemM==0 % stands for the edge is the boundary edge.
    pp_Row = ((elemP-1)*NTbases_test+1 : elemP*NTbases_test)' ... 
        * ones(1,NTbases_trial); % [NTbases_test x NTbases_trial]
    pp_Col = ones(NTbases_test, 1) ...
        *((elemP-1)*NTbases_trial+1 : elemP*NTbases_trial); % [NTbases_test x NTbases_trial]
    
    ROW = pp_Row;
    COL = pp_Col;
else
    pp_Row = ((elemP-1)*NTbases_test+1 : elemP*NTbases_test)' ... 
        * ones(1,NTbases_trial); % [NTbases_test x NTbases_trial]
    pp_Col = ones(NTbases_test, 1) ...
        *((elemP-1)*NTbases_trial+1 : elemP*NTbases_trial); % [NTbases_test x NTbases_trial]

    pm_Row = pp_Row; % [NTbases_test x NTbases_trial]
    pm_Col = ones(NTbases_test, 1) ...
        *((elemM-1)*NTbases_trial+1 : elemM*NTbases_trial); % [NTbases_test x NTbases_trial]

    mp_Row = ((elemM-1)*NTbases_test+1 : elemM*NTbases_test)' ...
        *ones(1,NTbases_trial); % [NTbases_test x NTbases_trial]
    mp_Col = pp_Col; % [NTbases_test x NTbases_trial]

    mm_Row = mp_Row; % [NTbases_test x NTbases_trial]
    mm_Col = pm_Col; % [NTbases_test x NTbases_trial]
    
    ROW = [pp_Row, pm_Row;
        mp_Row, mm_Row];
    COL = [pp_Col, pm_Col;
        mp_Col, mm_Col];
end % if elemM==0

end % function getROWandCOL
function C_vecInflowEdge = ...
    vecInflowEdgeAdapTransport(t, Coeffs_func, meshInfo, pde, formulaGauss1D, degreek_test)
%
%   %--------------------------------------------------
%       Just copy from vecInflowEdgeTransport.m
%   %--------------------------------------------------
%
%   In this function we compute the integral {.}\cdot[.] mat on INTERIOR edges.
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
%       meshInfo, the mesh information.
%       option, particularly, we need the \epsilon and the penalty parameter.
%       formulaGauss1D, the 1d Gauss quadrature formula, size: a matrix, [Npoints x 2],
%               the first column is the all Gauss-Points on [0,1],
%               the second is the weights of all Gauss-Points.
%       degreek_trial, the trial polynomial degree k.
%       degreek_test, the test polynomial degree k.
%
%   output:
%       C_vecInflowEdge, convection term inflow edges, [Nelems*NTbases_test x NTelems*NTbases_trial], 
%  
%
%   YcZhang 7/9/2017
%
%   Last modified 7/9/2017
%
%

% mesh information: inflow edges 
Nelems = meshInfo.Nelems;
inflowEdgeIndex = meshInfo.inflowEdgeIndex; % here we need the inflowEdgeIndex information.

vector_u_1 = Coeffs_func{3};
vector_u_2 = Coeffs_func{4};


% setting
NTbases_test = (degreek_test+1)*(degreek_test+2)/2;

C_vecInflowEdge = zeros(Nelems*NTbases_test,1);

for ii = 1:length(inflowEdgeIndex)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    inflowEindx = inflowEdgeIndex(ii); 
    ePoint1 = meshInfo.node(meshInfo.edge(inflowEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(inflowEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = meshInfo.areaEdge(inflowEindx);
    
    phyGpointsX = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = areaEdge * formulaGauss1D(:,2);
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
    elem1 = meshInfo.edge2elem(inflowEindx,1);
    elem2 = meshInfo.edge2elem(inflowEindx,2);
    local_e1 = meshInfo.edge2elem(inflowEindx,3);
    local_e2 = meshInfo.edge2elem(inflowEindx,4);
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
    
    elemP_xT = meshInfo.centroidElem(elemP,1);  
    elemP_yT = meshInfo.centroidElem(elemP,2);
    elemP_hT = meshInfo.hElem(elemP);
        %> the centroid(xing xin) and diameter of eleP.
        
    %-- setting elem_flag and get the value of vector_u
    elemP_flag = settingElemFlagAdapTransport(meshInfo,elemP); % here 
    valueP_vector_u_1 = vector_u_1(phyGpointsX, phyGpointsY, elemP_flag); 
        %> [Npoints x 1], this value_vector_u is a vector-valued.
    valueP_vector_u_2 = vector_u_2(phyGpointsX, phyGpointsY, elemP_flag); 
        %> [Npoints x 1], this value_vector_u is a vector-valued.
        
    value_f_inflow = pde.f_inflow(t, phyGpointsX, phyGpointsY);
        
    %-- subPart -- get the bases on quad    
    [elemP_testPb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_test);
        %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]
        
    %<<-- End Part I ---------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
    
    %% Part II, compute the inflow edges intgrations.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %
    %-- the inflow edges
    % get the position in the big mat    
    [rhsInflowRow, ~] = getROWandCOL(elemP, 0, 1, NTbases_test);
    
    vec_temp = elemP_nu(1) * elemP_testPb' * bsxfun(@times,valueP_vector_u_1.*phyGweights,value_f_inflow) ...
        + elemP_nu(2) * elemP_testPb' * bsxfun(@times,valueP_vector_u_2.*phyGweights,value_f_inflow);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    C_vecInflowEdge(rhsInflowRow,1) = C_vecInflowEdge(rhsInflowRow,1) + vec_temp;
    
end % for ii




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

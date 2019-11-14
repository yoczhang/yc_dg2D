function Kinter_nolinearTerm = ...
    matInterEdgeAverJumpNavierStokes2_nolinearTerm(lastStepUh, ...
    meshInfo, formulaGauss1D, degreek_u)
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   Just copy form matInterEdgeAverJumpNavierStokes1_nolinearTerm
%   %--------------------------------------------
%       dgNavierStokes2 using the Tensor-Stokes, 
%       the dgNavierStokes1 using the Saclar-Stokes.
%   %---------------------------------------------
%
%   In this function we compute the integral {.}\cdot[.] mat on INTERIOR edges.
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases_u denote the number of LOCAL velocity u bases on each K of Th.
%               NTbases_q denote the number of LOCAL pressure p bases on each K of Th.
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
%       degreek_u, the velocity polynomial degree k.
%       degreek_q, the pressure polynomial degree k.
%           %> here we consider that, in the Stokes equation, there will have VELOCITY u and PRESSURE p, 
%           %> so there have \int_e uv ds, then Trial_k=u_degree_k, Test_k=u_degree_k. 
%           %> And also have \int_e pv ds, then Trial_k=p_degree_k, Test_k=u_degree_k. 
%
%   output:
%       matInterEdge_uv, [Nelems*NTbases_u x NTelems*NTbases_u], 
%           \sum_e \int_e {(func1,func2)*grad u}\cdot[v] +  \epsilon\int_e {(func1,func2)*grad v}\cdot[u] + penalty \int_e [u]\cdot[v]. 
%       (matInterEdge_qv1, matInterEdge_qv2), [Nelems*NTbases_q x NTelems*NTbases_u],
%           \sum_e \int_e {p}\cdot[v]
%
%
%   YcZhang 17/10/2017
%
%   Last modified 17/10/2017
%
%

% mesh information: interior edges 
Nelems = meshInfo.Nelems;
interEdges = meshInfo.interEdgeIndex; % here we need the interEdges information.

% setting
NTbases_u = (degreek_u+1)*(degreek_u+2)/2;
dof_u1 = NTbases_u*Nelems;

Kinter_nolinearTerm = sparse(Nelems*NTbases_u,Nelems*NTbases_u);

for ii = 1:length(interEdges)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    interEindx = interEdges(ii); 
    ePoint1 = meshInfo.node(meshInfo.edge(interEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(interEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    eMidpoint = (ePoint1+ePoint2)/2;
    
    phyGpointsX = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = meshInfo.areaEdge(interEindx) * formulaGauss1D(:,2); 
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
        
    %-- subPart -- get the bases on quad
    [u_elemP_Pb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_u);
        %> u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, [Npoints x NTbases_u]
    [u_elemM_Pb, ~, ~] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_u);
        %> u_eleM_Pb, u_eleM_Pbx, u_eleM_Pby, [Npoints x NTbases_u]
    [eMid_u_elemP_Pb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, eMidpoint(1), eMidpoint(2), degreek_u);
        %> u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, [1 x NTbases_u]
        
    lastStepElemP_uh1 =  lastStepUh((elemP-1)*NTbases_u+1:elemP*NTbases_u);
    lastStepElemP_uh2 =  lastStepUh(dof_u1+(elemP-1)*NTbases_u+1:dof_u1+elemP*NTbases_u);
    lastStepElemM_uh1 =  lastStepUh((elemM-1)*NTbases_u+1:elemM*NTbases_u);
    lastStepElemM_uh2 =  lastStepUh(dof_u1+(elemM-1)*NTbases_u+1:dof_u1+elemM*NTbases_u);
    
    valueLastStepElemP_uh1 = u_elemP_Pb*lastStepElemP_uh1; % [Npoints x 1]
    valueLastStepElemP_uh2 = u_elemP_Pb*lastStepElemP_uh2; % [Npoints x 1]
    valueLastStepElemM_uh1 = u_elemM_Pb*lastStepElemM_uh1; % [Npoints x 1]
    valueLastStepElemM_uh2 = u_elemM_Pb*lastStepElemM_uh2; % [Npoints x 1]
    
    eMid_valueLastStepElemP_uh1 = eMid_u_elemP_Pb*lastStepElemP_uh1; % [Npoints x 1]
    eMid_valueLastStepElemP_uh2 = eMid_u_elemP_Pb*lastStepElemP_uh2; % [Npoints x 1]
    %<<-- End Part I ---------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
    
    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %
    %-- subPart 1 -- the AverJump of velotivy trial u and velotivy test v.
    % get the position in the big mat
    [uv_Row, uv_Col] = getROWandCOL(elemP, elemM, NTbases_u, NTbases_u);
        
    %-- subPart 3 -- get the mat of nolinear term 
    averLastStep_uh1 = 0.5*(valueLastStepElemP_uh1 + valueLastStepElemM_uh1); % [Npoints x 1]
    averLastStep_uh2 = 0.5*(valueLastStepElemP_uh2 + valueLastStepElemM_uh2); % [Npoints x 1]
    jumpLastStep_uh = valueLastStepElemP_uh1*elemP_nu(1) + valueLastStepElemP_uh2*elemP_nu(2) ...
        + valueLastStepElemM_uh1*elemM_nu(1) + valueLastStepElemM_uh2*elemM_nu(2);
    
    [K_nolinear_temp1, K_nolinear_temp2]  = getMatOnInterEdge_nolinear( phyGweights, ...
        eMid_valueLastStepElemP_uh1, eMid_valueLastStepElemP_uh2, ...
        averLastStep_uh1, averLastStep_uh2, jumpLastStep_uh, ...
        elemP_nu, elemM_nu, ...
        u_elemP_Pb, ...
        u_elemM_Pb);
    
    Kinter_nolinearTerm = Kinter_nolinearTerm ...
        + sparse(uv_Row(:), uv_Col(:), ...
        -K_nolinear_temp1(:)-0.5*K_nolinear_temp2(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    
    
    %<<-- End Part II --------------------------------------------------------------------------------------
end % for ii


end % function matInterEdgeAverJump







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


function [K_nolinear_temp1, K_nolinear_temp2] = getMatOnInterEdge_nolinear( phyGaussweights, ...
    eMid_valueLastStepElemP_uh1, eMid_valueLastStepElemP_uh2, ...
    averLastStep_uh1, averLastStep_uh2, jumpLastStep_uh, ...
    elemP_nu, elemM_nu, ...
    u_elemP_Pb, ...
    u_elemM_Pb)
%
%
%
%
%
%   YcZhang 21/10/2017
%   Last modified 21/10/2017
%

NTbases = size(u_elemP_Pb,2);

mat0 = zeros(NTbases, NTbases);

if eMid_valueLastStepElemP_uh1*elemP_nu(1) + eMid_valueLastStepElemP_uh2*elemP_nu(2) <= 0
    K1_pp = u_elemP_Pb' ...
        * bsxfun(@times, ...
        (elemP_nu(1)*averLastStep_uh1+ elemP_nu(2)*averLastStep_uh2).*phyGaussweights, ...
        u_elemP_Pb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    K1_pm = u_elemP_Pb' ...
        * bsxfun(@times, ...
        (elemM_nu(1)*averLastStep_uh1+ elemM_nu(2)*averLastStep_uh2).*phyGaussweights, ...
        u_elemM_Pb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    K_nolinear_temp1 = [K1_pp, K1_pm; mat0, mat0];
else
    K1_mp = u_elemM_Pb' ...
        * bsxfun(@times, ...
        (elemP_nu(1)*averLastStep_uh1+ elemP_nu(2)*averLastStep_uh2).*phyGaussweights, ...
        u_elemP_Pb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    K1_mm = u_elemM_Pb' ...
        * bsxfun(@times, ...
        (elemM_nu(1)*averLastStep_uh1+ elemM_nu(2)*averLastStep_uh2).*phyGaussweights, ...
        u_elemM_Pb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    K_nolinear_temp1 = [mat0, mat0; K1_mp, K1_mm];
end % if 

K2_pp = 0.5*u_elemP_Pb' * bsxfun(@times, jumpLastStep_uh.*phyGaussweights, u_elemP_Pb);
K2_mm = 0.5*u_elemM_Pb' * bsxfun(@times, jumpLastStep_uh.*phyGaussweights, u_elemM_Pb);
K_nolinear_temp2 = [K2_pp, mat0; mat0, K2_mm];    


end % function K_nolinear_temp = getMatOnInterEdge_nolinear

%--------------------------------------------------------------------------------------------------------------------%
%----------- end sub function ----------------------------------------------------------------------------------%

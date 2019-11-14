function [Krhs_nolinear_gD1v,Krhs_nolinear_gD2v] = ...
    vecDirichletEdgeAverJumpNavierStokes1_nolinearTerm(func_gD1, func_gD2, lastStepUh, ...
    meshInfo, formulaGauss1D, degreek_u)
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       Using the Saclar-Stokes, 
%       the dgNavierStokes2 using the Tensor-Stokes.
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
%       matDirichletEdge_uv, [Nelems*NTbases_u x NTelems*NTbases_u], 
%           \sum_e \int_e {(func1,func2)*grad u}\cdot[v] +  \epsilon\int_e {(func1,func2)*grad v}\cdot[u] + penalty \int_e [u]\cdot[v]. 
%       (matDirichletEdge_qv1, matDirichletEdge_qv2), [Nelems*NTbases_q x NTelems*NTbases_u],
%           \sum_e \int_e {p}\cdot[v]
%       (rhs_gD1,rhs_gD2), [Nelems*NTbases_u x 1], \sum_Diri_e \int_e gD v 
%       rhs_div, [Nelems*NTbases_q x 1]
%
%
%   YcZhang 15/8/2017
%
%   Last modified 15/8/2017
%
%

% mesh information: interior edges 
Nelems = meshInfo.Nelems;
DirichletEdges = meshInfo.DirichletEdgeIndex; % here we need the DirichletEdges information.

% setting
NTbases_u = (degreek_u+1)*(degreek_u+2)/2;
dof_u1 = NTbases_u*Nelems;

Krhs_nolinear_gD1v = zeros(Nelems*NTbases_u,1);
Krhs_nolinear_gD2v = zeros(Nelems*NTbases_u,1);

for ii = 1:length(DirichletEdges)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    DirichletEindx = DirichletEdges(ii); 
    ePoint1 = meshInfo.node(meshInfo.edge(DirichletEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(DirichletEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    eMidpoint = (ePoint1+ePoint2)/2;
    
    phyGpointsX = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = meshInfo.areaEdge(DirichletEindx) * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
    elem1 = meshInfo.edge2elem(DirichletEindx,1);
    elem2 = meshInfo.edge2elem(DirichletEindx,2);
    local_e1 = meshInfo.edge2elem(DirichletEindx,3);
    local_e2 = meshInfo.edge2elem(DirichletEindx,4);
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
        
    %-- get the value of Coeffs func at Gauss Points
    valuefunc_gD1 = func_gD1(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valuefunc_gD2 = func_gD2(phyGpointsX, phyGpointsY); % [Npoints x 1]
    
    %-- subPart -- get the bases on quad
    [u_elemP_Pb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_u);
        %> u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, [Npoints x NTbases_u]
        
    [eMid_u_elemP_Pb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, eMidpoint(1), eMidpoint(2), degreek_u);
        %> u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, [1 x NTbases_u]
       
    lastStepElemP_uh1 =  lastStepUh((elemP-1)*NTbases_u+1:elemP*NTbases_u);
    lastStepElemP_uh2 =  lastStepUh(dof_u1+(elemP-1)*NTbases_u+1:dof_u1+elemP*NTbases_u);

    valueLastStepElemP_uh1 = u_elemP_Pb*lastStepElemP_uh1; % [Npoints x 1]
    valueLastStepElemP_uh2 = u_elemP_Pb*lastStepElemP_uh2; % [Npoints x 1]
    
    eMid_valueLastStepElemP_uh1 = eMid_u_elemP_Pb*lastStepElemP_uh1; % [1 x 1]
    eMid_valueLastStepElemP_uh2 = eMid_u_elemP_Pb*lastStepElemP_uh2; % [1 x 1]
    

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
    [rhsDirRow, ~] = getROWandCOL(elemP, 0, 1, NTbases_u);

    %-- subPart 5 -- get the mat of Dirichlet nolinear term 
    averLastStep_uh1 = valueLastStepElemP_uh1; % [Npoints x 1]
    averLastStep_uh2 = valueLastStepElemP_uh2; % [Npoints x 1]
    jumpLastStep_uh = valueLastStepElemP_uh1*elemP_nu(1) + valueLastStepElemP_uh2*elemP_nu(2);
    
    [K_nolinear_gD1_temp1,K_nolinear_gD2_temp1, ...
        K_nolinear_gD1_temp2,K_nolinear_gD2_temp2]  = getVecOnDirichletEdge_nolinear( phyGweights, ...
        valuefunc_gD1, valuefunc_gD2, ...
        eMid_valueLastStepElemP_uh1, eMid_valueLastStepElemP_uh2, ...
        averLastStep_uh1, averLastStep_uh2, jumpLastStep_uh, ...
        elemP_nu, ...
        u_elemP_Pb);
    
    Krhs_nolinear_gD1v(rhsDirRow,1) = Krhs_nolinear_gD1v(rhsDirRow,1) ...
        +( -K_nolinear_gD1_temp1(:)-0.5*K_nolinear_gD1_temp2(:) );
    Krhs_nolinear_gD2v(rhsDirRow,1) = Krhs_nolinear_gD2v(rhsDirRow,1) ...
        +( -K_nolinear_gD2_temp1(:)-0.5*K_nolinear_gD2_temp2(:) );
    
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



function [K_nolinear_gD1_temp1,K_nolinear_gD2_temp1, ...
    K_nolinear_gD1_temp2,K_nolinear_gD2_temp2] = getVecOnDirichletEdge_nolinear( phyGaussweights, ...
    valuefunc_gD1, valuefunc_gD2, ...
    eMid_valueLastStepElemP_uh1, eMid_valueLastStepElemP_uh2, ...
    averLastStep_uh1, averLastStep_uh2, jumpLastStep_uh, ...
    elemP_nu, ...
    u_elemP_Pb)
%
%
%
%
%
%   YcZhang 21/10/2017
%   Last modified 21/10/2017
%


if eMid_valueLastStepElemP_uh1*elemP_nu(1) + eMid_valueLastStepElemP_uh2*elemP_nu(2) <= 0
    K1_gD1_pp = u_elemP_Pb' ...
        * bsxfun(@times, ...
        (elemP_nu(1)*averLastStep_uh1+ elemP_nu(2)*averLastStep_uh2).*phyGaussweights, ...
        valuefunc_gD1);
    K1_gD2_pp = u_elemP_Pb' ...
        * bsxfun(@times, ...
        (elemP_nu(1)*averLastStep_uh1+ elemP_nu(2)*averLastStep_uh2).*phyGaussweights, ...
        valuefunc_gD2);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    K_nolinear_gD1_temp1 = K1_gD1_pp;
    K_nolinear_gD2_temp1 = K1_gD2_pp;
else
    NTbases = size(u_elemP_Pb,2);
    vec0 = zeros(NTbases, 1);
    K_nolinear_gD1_temp1 = vec0;
    K_nolinear_gD2_temp1 = vec0;
end % if 

K2_gD1_pp = u_elemP_Pb' * bsxfun(@times, jumpLastStep_uh.*phyGaussweights, valuefunc_gD1);
K2_gD2_pp = u_elemP_Pb' * bsxfun(@times, jumpLastStep_uh.*phyGaussweights, valuefunc_gD2);
K_nolinear_gD1_temp2 = K2_gD1_pp;    
K_nolinear_gD2_temp2 = K2_gD2_pp; 

end % function K_nolinear_temp = getMatOnInterEdge_nolinear

%--------------------------------------------------------------------------------------------------------------------%
%----------- end sub function ----------------------------------------------------------------------------------%

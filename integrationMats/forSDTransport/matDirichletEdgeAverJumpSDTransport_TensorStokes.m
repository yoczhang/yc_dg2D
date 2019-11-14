function [matDirichletEdge_u1v1,matDirichletEdge_u1v2,matDirichletEdge_u2v1,matDirichletEdge_u2v2,matDirichletEdge_qv1,matDirichletEdge_qv2,rhs_gD1,rhs_gD2,rhs_div] = ...
    matDirichletEdgeAverJumpSDTransport_TensorStokes(Coeffs_func, func_gD1, func_gD2, meshInfo, option, formulaGauss1D, degreek_u, degreek_q)
%
%   %---------------------------------------------------------------------
%       Just copy from: matDirichletEdgeAverJumpTensorStokes.m
%   %----------------------------------------------------------------------
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
%       matDirichletEdge_u..v.., [Nelems*NTbases_u x NTelems*NTbases_u], 
%           \sum_e \int_e {(func1,func2)*grad u}\cdot[v] +  \epsilon\int_e {(func1,func2)*grad v}\cdot[u] + penalty \int_e [u]\cdot[v]. 
%       (matDirichletEdge_qv1, matDirichletEdge_qv2), [Nelems*NTbases_q x NTelems*NTbases_u],
%           \sum_e \int_e {p}\cdot[v]
%       (rhs_gD1,rhs_gD2), [Nelems*NTbases_u x 1], \sum_Diri_e \int_e gD v 
%       rhs_div, [Nelems*NTbases_q x 1]
%
%
%   YcZhang 25/8/2017
%
%   Last modified 26/8/2017
%
%

% mesh information: interior edges 
Nelems = meshInfo.Nelems;
DirichletEdges = meshInfo.DirichletEdgeIndex; % here we need the DirichletEdges information.

p_epsilon = option.p_epsilon;
p_sigma = option.p_sigma;
p_beta = option.p_beta;

Coeff_k11 = Coeffs_func{1};
Coeff_k22 = Coeffs_func{2};

% setting
NTbases_u = (degreek_u+1)*(degreek_u+2)/2;
NTbases_q = (degreek_q+1)*(degreek_q+2)/2;

matDirichletEdge_u1v1 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
matDirichletEdge_u1v2 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
matDirichletEdge_u2v1 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
matDirichletEdge_u2v2 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
matDirichletEdge_qv1 = sparse(Nelems*NTbases_u,Nelems*NTbases_q);
matDirichletEdge_qv2 = sparse(Nelems*NTbases_u,Nelems*NTbases_q);
rhs_gD1 = zeros(Nelems*NTbases_u,1);
rhs_gD2 = zeros(Nelems*NTbases_u,1);
rhs_div = zeros(Nelems*NTbases_q,1);


for ii = 1:length(DirichletEdges)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    DirichletEindx = DirichletEdges(ii); 
    ePoint1 = meshInfo.node(meshInfo.edge(DirichletEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(DirichletEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = meshInfo.areaEdge(DirichletEindx);
    
    phyGpointsX = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = meshInfo.areaEdge(DirichletEindx) * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
    ele1 = meshInfo.edge2elem(DirichletEindx,1);
    ele2 = meshInfo.edge2elem(DirichletEindx,2);
    local_e1 = meshInfo.edge2elem(DirichletEindx,3);
    local_e2 = meshInfo.edge2elem(DirichletEindx,4);
        %> [n1, n2]=edge2elem(k,1:2) stands for the n1-th elem and n2-th elem share the k-th edge.
        %> local_e1=edge2elem(k,3) stands for the local edge index in n1-th elem of the k-th edge .
        %> also, local_e2=edge2elem(k,4) stands for the local edge index in n2-th elem of the k-th edge .
    eleP = max(ele1, ele2); % for the boundary edge, there will be one zero of (ele1,ele2), for let the bigger number of element as the element +.
    
    if eleP == ele1
        eIndx_eleP = local_e1;
    else
        eIndx_eleP = local_e2;
    end 
    
    eleP_nu = meshInfo.nuEdge0Elem{eleP}(:,eIndx_eleP); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleP.
    
    eleP_xT = meshInfo.centroidElem(eleP,1);  
    eleP_yT = meshInfo.centroidElem(eleP,2);
    eleP_hT = meshInfo.hElem(eleP);
        %> the centroid(xing xin) and diameter of eleP.
        
    %-- get the value of Coeffs func at Gauss Points
    valueCoeff_k11 = Coeff_k11(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valueCoeff_k22 = Coeff_k22(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valuefunc_gD1 = func_gD1(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valuefunc_gD2 = func_gD2(phyGpointsX, phyGpointsY); % [Npoints x 1]
    
    %-- subPart -- get the bases on quad
    [u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby] = ...
        localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, degreek_u);
        %> u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, [Npoints x NTbases_u]
    
    [q_eleP_Pb, ~, ~] = ...
        localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, degreek_q);
        %> q_eleP_Pb, q_eleP_Pbx, q_eleP_Pby, [Npoints x NTbases_q]
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
    [uv_Row, uv_Col] = getROWandCOL(eleP, 0, NTbases_u, NTbases_u);
    
    % get the matrix
    [uv_A1, uv_A2, uv_A3] = getMatOnDirichletEdge_uv(phyGweights, ...
        valueCoeff_k11, valueCoeff_k22, ...
        eleP_nu, ...
        u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, ...
        u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby);
        %> uv_A1, uv_A2, uv_A3, structure-type, uv_A1.xx [2*NTbases_u  x 2*NTbases_u]
        %> uv_A1.xx=[A1_pp];
        %> the same to uv_A2, uv_A3.
    
    matDirichletEdge_u1v1 = matDirichletEdge_u1v1 ...
        + sparse(uv_Row(:), uv_Col(:), ...
        - 2*uv_A1.A1_u1v1(:) + p_epsilon*2*uv_A2.A2_u1v1(:) + (p_sigma/areaEdge^p_beta)*uv_A3.A3_u1v1(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    
    matDirichletEdge_u1v2 = matDirichletEdge_u1v2 ...
        + sparse(uv_Row(:), uv_Col(:), ...
        - 2*uv_A1.A1_u1v2(:) + p_epsilon*2*uv_A2.A2_u1v2(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    
    matDirichletEdge_u2v1 = matDirichletEdge_u2v1 ...
        + sparse(uv_Row(:), uv_Col(:), ...
        - 2*uv_A1.A1_u2v1(:) + p_epsilon*2*uv_A2.A2_u2v1(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    
    matDirichletEdge_u2v2 = matDirichletEdge_u2v2 ...
        + sparse(uv_Row(:), uv_Col(:), ...
        - 2*uv_A1.A1_u2v2(:) + p_epsilon*2*uv_A2.A2_u2v2(:) + (p_sigma/areaEdge^p_beta)*uv_A3.A3_u2v2(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    
    
    %-- subPart 2 -- the AverJump of pressure trial p and velotivy test v.
    % get the position in the big mat
    [qv_Row, qv_Col] = getROWandCOL(eleP, 0, NTbases_q, NTbases_u);
    
    % get the matrix
    [qv_A1, qv_A2] = getMatOnDirichletEdge_qv(phyGweights, ...
        eleP_nu, ...
        q_eleP_Pb, ...
        u_eleP_Pb);
        %> qv_A1, qv_A2, [NTbases_u x NTbases_q]
        %> qv_A1=[A1_pp, A1_pm;
        %>            A1_mp, A1_mm];
        %> the same to qv_A2.
    
    matDirichletEdge_qv1 = matDirichletEdge_qv1 ...
        + sparse(qv_Row(:), qv_Col(:), qv_A1(:), Nelems*NTbases_u, Nelems*NTbases_q);
    matDirichletEdge_qv2 = matDirichletEdge_qv2 ...
        + sparse(qv_Row(:), qv_Col(:), qv_A2(:), Nelems*NTbases_u, Nelems*NTbases_q);
    
    
    %-- subPart 3 -- the rhs Dirichlet edges
    %-- get the rhs Dirchlet vector
    [rhsDirRow, ~] = getROWandCOL(eleP, 0, 1, NTbases_u);
    
%     r1 = eleP_nu(1) * u_eleP_Pbx' *  (valueCoeff_k11.*valuefunc_gD1.*phyGweights) ...
%         + eleP_nu(2) * u_eleP_Pby' * (valueCoeff_k22.*valuefunc_gD1.*phyGweights);
%         %> [NTbases_test x Npoints] * [Npoints x 1] = [NTbases_test x 1].
    r1 = eleP_nu(1) * u_eleP_Pbx' *  (valueCoeff_k11.*valuefunc_gD1.*phyGweights) ...
        + 0.5 * eleP_nu(2) * u_eleP_Pby' *  (valueCoeff_k11.*valuefunc_gD1.*phyGweights) ...
        + 0.5 * eleP_nu(1) * u_eleP_Pby' *  (valueCoeff_k11.*valuefunc_gD2.*phyGweights);
    r1 = p_epsilon*2*r1 + (p_sigma/areaEdge^p_beta)*u_eleP_Pb' * (valuefunc_gD1.*phyGweights);
    rhs_gD1(rhsDirRow,1) = rhs_gD1(rhsDirRow,1) + r1;
    
%     r2 = eleP_nu(1) * u_eleP_Pbx' *  (valueCoeff_k11.*valuefunc_gD2.*phyGweights) ...
%         + eleP_nu(2) * u_eleP_Pby' * (valueCoeff_k22.*valuefunc_gD2.*phyGweights);
%         %> [NTbases_test x Npoints] * [Npoints x 1] = [NTbases_test x 1].
    r2 = 0.5 * eleP_nu(2) * u_eleP_Pbx' *  (valueCoeff_k11.*valuefunc_gD1.*phyGweights) ...
        + 0.5 * eleP_nu(1) * u_eleP_Pbx' *  (valueCoeff_k11.*valuefunc_gD2.*phyGweights) ...
        + eleP_nu(2) * u_eleP_Pby' *  (valueCoeff_k11.*valuefunc_gD2.*phyGweights);
    r2 = p_epsilon*2*r2 + (p_sigma/areaEdge^p_beta)*u_eleP_Pb' * (valuefunc_gD2.*phyGweights);
    rhs_gD2(rhsDirRow,1) = rhs_gD2(rhsDirRow,1) + r2;
    
    
    %%-- subPart 4 -- the rhs divergence-free Dirichlet term
    [rhsdivDirRow, ~] = getROWandCOL(eleP,0,1,NTbases_q);
    r3 = eleP_nu(1) * q_eleP_Pb' * (valueCoeff_k11.*valuefunc_gD1.*phyGweights) ...
        + eleP_nu(2) * q_eleP_Pb' * (valueCoeff_k22.*valuefunc_gD2.*phyGweights);
    rhs_div(rhsdivDirRow,1) = rhs_div(rhsdivDirRow,1) + r3;
    
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


%% the sub function 
function [A1, A2, A3] = getMatOnDirichletEdge_uv(phyGaussweights, ...
    valueCoeff_1, valueCoeff_2, ...
    eleP_nu, ...
    eleP_trialPb, eleP_trialPbx, eleP_trialPby, ...
    eleP_testPb, eleP_testPbx, eleP_testPby)
%
%   get the Aver{.} and Jump[.] mat on Interior edges.
%
%   int_Agradtest_Jtrial --> \int_e {(func1,func2)\cdot grad_trial}\cdot[test].
%   int_Agradtest_Jtrial --> \int_e {(func1,func2)\cdot grad_test}\cdot[trial].
%   int_Jtrial_Jtest --> \int_e [trial]\cdot[test].
%
%   YcZhang 26/8/2017
%   Last modified 26/8/2017
%

    %-- subPart I, \int_e {(func1,func2)\cdot grad_trial}\cdot[test]
    % {grad_trial}\cdot[test] = 1/2(grad_trial^p + grad_trial^m)\cdot(test^p nu^p + test^m nu^m)
    %AgradtestJtrial_pp = zeros(NTbases_test, NTbases_trial);    
    A1_u1v1_PP = ( eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + 0.5 * eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
    A1_u2v1_PP = ( 0.5 * eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPbx) );
    A1_u1v2_PP = ( 0.5 * eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
    A1_u2v2_PP = ( 0.5 * eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1.A1_u1v1 = A1_u1v1_PP;
    A1.A1_u2v1 = A1_u2v1_PP;
    A1.A1_u1v2 = A1_u1v2_PP;
    A1.A1_u2v2 = A1_u2v2_PP;
        
    
    %-- subPart II, \int_e {(func1,func2)\cdot grad_test}\cdot[trial]
    A2_u1v1_PP = ( eleP_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + 0.5 * eleP_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
    A2_u1v2_PP = ( 0.5 * eleP_nu(2) * eleP_testPbx' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
    A2_u2v1_PP = ( 0.5 * eleP_nu(1) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
    A2_u2v2_PP = ( 0.5 * eleP_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + eleP_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2.A2_u1v1 = A2_u1v1_PP;
    A2.A2_u2v1 = A2_u2v1_PP;
    A2.A2_u1v2 = A2_u1v2_PP;
    A2.A2_u2v2 = A2_u2v2_PP;
    
        
    %-- subPart III, \int_e [trial]\cdot[test]    
    A3_u1v1_PP = eleP_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3.A3_u1v1 = A3_u1v1_PP;
    A3.A3_u2v2 = A3_u1v1_PP;
end % function getMatOnEdge_uv


%% the sub function 
function [int_Atrial_Jtest_1, int_Atrial_Jtest_2] = getMatOnDirichletEdge_qv(phyGaussweights, ...
    eleP_nu, ...
    eleP_trialPb, ...
    eleP_testPb)
%
%   get the Aver{.} and Jump[.] mat on Interior edges.
%
%   int_Agradtest_Jtrial --> \int_e {(func1,func2)\cdot grad_trial}\cdot[test].
%   int_Agradtest_Jtrial --> \int_e {(func1,func2)\cdot grad_test}\cdot[trial].
%   int_Jtrial_Jtest --> \int_e [trial]\cdot[test].
%
%   YcZhang 26/8/2017
%   Last modified 26/8/2017
%

    %-- subPart I, \int_e {(func1,func2)\cdot grad_trial}\cdot[test]
    % {grad_trial}\cdot[test] = 1/2(grad_trial^p + grad_trial^m)\cdot(test^p nu^p + test^m nu^m)
    %AgradtestJtrial_pp = zeros(NTbases_test, NTbases_trial);    
    A1_pp = eleP_nu(1) * eleP_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    int_Atrial_Jtest_1 = A1_pp;
        
    
    %-- subPart II, \int_e {(func1,func2)\cdot grad_test}\cdot[trial]
    A2_pp = eleP_nu(2) * eleP_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    int_Atrial_Jtest_2 = A2_pp;    

end % function getMatOnEdge_qv

%--------------------------------------------------------------------------------------------------------------------%
%----------- end sub function ----------------------------------------------------------------------------------%

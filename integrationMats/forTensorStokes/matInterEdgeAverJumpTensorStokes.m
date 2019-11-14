function [matInterEdge_u1v1, matInterEdge_u1v2, matInterEdge_u2v1, matInterEdge_u2v2, matInterEdge_qv1, matInterEdge_qv2] = ...
    matInterEdgeAverJumpTensorStokes(Coeffs_func, meshInfo, option, formulaGauss1D, degreek_u, degreek_p)
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
%       matInterEdge_uXvX, (X denotes 1 or 2), [Nelems*NTbases_u x NTelems*NTbases_u], 
%           \sum_e \int_e {(func1,func2)*grad u}\cdot[v] +  \epsilon\int_e {(func1,func2)*grad v}\cdot[u] + penalty \int_e [u]\cdot[v]. 
%       (matInterEdge_qv1, matInterEdge_qv2), [Nelems*NTbases_q x NTelems*NTbases_u],
%           \sum_e \int_e {p}\cdot[v]
%
%
%   YcZhang 25/8/2017
%
%   Last modified 26/8/2017
%
%

% mesh information: interior edges 
Nelems = meshInfo.Nelems;
interEdges = meshInfo.interEdgeIndex; % here we need the interEdges information.

p_epsilon = option.p_epsilon;
p_sigma = option.p_sigma;
p_beta = option.p_beta;

Coeff_k11 = Coeffs_func{1};
Coeff_k22 = Coeffs_func{2};

% setting
NTbases_u = (degreek_u+1)*(degreek_u+2)/2;
NTbases_p = (degreek_p+1)*(degreek_p+2)/2;

matInterEdge_u1v1 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
matInterEdge_u1v2 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
matInterEdge_u2v1 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
matInterEdge_u2v2 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
matInterEdge_qv1 = sparse(Nelems*NTbases_u,Nelems*NTbases_p);
matInterEdge_qv2 = sparse(Nelems*NTbases_u,Nelems*NTbases_p);

for ii = 1:length(interEdges)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    interEindx = interEdges(ii); 
    ePoint1 = meshInfo.node(meshInfo.edge(interEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(interEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = meshInfo.areaEdge(interEindx);
    
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
        
    %-- get the value of Coeffs func at Gauss Points
    valueCoeff_k11 = Coeff_k11(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valueCoeff_k22 = Coeff_k22(phyGpointsX, phyGpointsY); % [Npoints x 1]
    
    %-- subPart -- get the bases on quad
    [u_elemP_Pb, u_elemP_Pbx, u_elemP_Pby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_u);
        %> u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, [Npoints x NTbases_u]
    [u_elemM_Pb, u_elemM_Pbx, u_elemM_Pby] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_u);
        %> u_eleM_Pb, u_eleM_Pbx, u_eleM_Pby, [Npoints x NTbases_u]
    
    [p_elemP_Pb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_p);
        %> q_eleP_Pb, q_eleP_Pbx, q_eleP_Pby, [Npoints x NTbases_q]
    [p_elemM_Pb, ~, ~] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_p);
        %> q_eleM_Pb, q_eleM_Pbx, q_eleM_Pby, [Npoints x NTbases_q]
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
    
    % get the matrix
    [uv_A1, uv_A2, uv_A3] = getMatOnInterEdge_uv(phyGweights, ...
        valueCoeff_k11, valueCoeff_k22, ...
        elemP_nu, ...
        u_elemP_Pb, u_elemP_Pbx, u_elemP_Pby, ...
        u_elemP_Pb, u_elemP_Pbx, u_elemP_Pby, ...
        elemM_nu, ...
        u_elemM_Pb, u_elemM_Pbx, u_elemM_Pby, ...
        u_elemM_Pb, u_elemM_Pbx, u_elemM_Pby);
        %> uv_A1, uv_A2, uv_A3, structure-type, uv_A1.xx [2*NTbases_u  x 2*NTbases_u]
        %> uv_A1.xx = [A1_pp, A1_pm; 
        %>            A1_mp, A1_mm];
        %> the same to uv_A2, uv_A3.
    
    matInterEdge_u1v1 = matInterEdge_u1v1 ...
        + sparse(uv_Row(:), uv_Col(:), ...
        - 2*uv_A1.A1_u1v1(:) + p_epsilon*2*uv_A2.A2_u1v1(:) + (p_sigma/areaEdge^p_beta)*uv_A3.A3_u1v1(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    
    matInterEdge_u1v2 = matInterEdge_u1v2 ...
        + sparse(uv_Row(:), uv_Col(:), ...
        - 2*uv_A1.A1_u1v2(:) + p_epsilon*2*uv_A2.A2_u1v2(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    
    matInterEdge_u2v1 = matInterEdge_u2v1 ...
        + sparse(uv_Row(:), uv_Col(:), ...
        - 2*uv_A1.A1_u2v1(:) + p_epsilon*2*uv_A2.A2_u2v1(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    
    matInterEdge_u2v2 = matInterEdge_u2v2 ...
        + sparse(uv_Row(:), uv_Col(:), ...
        - 2*uv_A1.A1_u2v2(:) + p_epsilon*2*uv_A2.A2_u2v2(:) + (p_sigma/areaEdge^p_beta)*uv_A3.A3_u2v2(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    
    
    %-- subPart 2 -- the AverJump of pressure trial p and velotivy test v.
    % get the position in the big mat
    [qv_Row, qv_Col] = getROWandCOL(elemP, elemM, NTbases_p, NTbases_u);
    
    % get the matrix
    [qv_A1, qv_A2] = getMatOnInterEdge_qv(phyGweights, ...
        elemP_nu, ...
        p_elemP_Pb, ...
        u_elemP_Pb, ...
        elemM_nu, ...
        p_elemM_Pb, ...
        u_elemM_Pb);
        %> qv_A1, qv_A2, [2*NTbases_u x 2*NTbases_q]
        %> qv_A1=[A1_pp, A1_pm; 
        %>            A1_mp, A1_mm];
        %> the same to qv_A2.
    
    matInterEdge_qv1 = matInterEdge_qv1 ...
        + sparse(qv_Row(:), qv_Col(:), qv_A1(:), Nelems*NTbases_u, Nelems*NTbases_p);
    matInterEdge_qv2 = matInterEdge_qv2 ...
        + sparse(qv_Row(:), qv_Col(:), qv_A2(:), Nelems*NTbases_u, Nelems*NTbases_p);
    
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
function [A1, A2, A3] = getMatOnInterEdge_uv(phyGaussweights, ...
    valueCoeff_1, valueCoeff_2, ...
    eleP_nu, ...
    eleP_trialPb, eleP_trialPbx, eleP_trialPby, ...
    eleP_testPb, eleP_testPbx, eleP_testPby, ...
    eleM_nu, ...
    eleM_trialPb, eleM_trialPbx, eleM_trialPby, ...
    eleM_testPb, eleM_testPbx, eleM_testPby)
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
    A1_u1v1_PP = 0.5*( eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + 0.5 * eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
    A1_u2v1_PP = 0.5*( 0.5 * eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPbx) );
    A1_u1v2_PP = 0.5*( 0.5 * eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
    A1_u2v2_PP = 0.5*( 0.5 * eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_u1v1_MP = 0.5*( eleM_nu(1) * eleM_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + 0.5 * eleM_nu(2) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
    A1_u2v1_MP = 0.5*( 0.5 * eleM_nu(2) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPbx) );
    A1_u1v2_MP = 0.5*( 0.5 * eleM_nu(1) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
    A1_u2v2_MP = 0.5*( 0.5 * eleM_nu(1) * eleM_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + eleM_nu(2) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_u1v1_PM = 0.5*( eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPbx) ...
        + 0.5 * eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPby) );
    A1_u2v1_PM = 0.5*( 0.5 * eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPbx) );
    A1_u1v2_PM = 0.5*( 0.5 * eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPby) );
    A1_u2v2_PM = 0.5*( 0.5 * eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPbx) ...
        + eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_u1v1_MM = 0.5*( eleM_nu(1) * eleM_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPbx) ...
        + 0.5 * eleM_nu(2) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPby) );
    A1_u2v1_MM = 0.5*( 0.5 * eleM_nu(2) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPbx) );
    A1_u1v2_MM = 0.5*( 0.5 * eleM_nu(1) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPby) );
    A1_u2v2_MM = 0.5*( 0.5 * eleM_nu(1) * eleM_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPbx) ...
        + eleM_nu(2) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1.A1_u1v1 = [A1_u1v1_PP, A1_u1v1_PM;
        A1_u1v1_MP, A1_u1v1_MM];
    A1.A1_u2v1 = [A1_u2v1_PP, A1_u2v1_PM;
        A1_u2v1_MP, A1_u2v1_MM];
    A1.A1_u1v2 = [A1_u1v2_PP, A1_u1v2_PM;
        A1_u1v2_MP, A1_u1v2_MM];
    A1.A1_u2v2 = [A1_u2v2_PP, A1_u2v2_PM;
        A1_u2v2_MP, A1_u2v2_MM];
    
    %-- subPart II, \int_e {(func1,func2)\cdot grad_test}\cdot[trial]    
    A2_u1v1_PP = 0.5*( eleP_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + 0.5 * eleP_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
    A2_u1v2_PP = 0.5*( 0.5 * eleP_nu(2) * eleP_testPbx' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
    A2_u2v1_PP = 0.5*( 0.5 * eleP_nu(1) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
    A2_u2v2_PP = 0.5*( 0.5 * eleP_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + eleP_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_u1v1_PM = 0.5*( eleM_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPb) ...
        + 0.5 * eleM_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
    A2_u1v2_PM = 0.5*( 0.5 * eleM_nu(2) * eleP_testPbx' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
    A2_u2v1_PM = 0.5*( 0.5 * eleM_nu(1) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
    A2_u2v2_PM = 0.5*( 0.5 * eleM_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPb) ...
        + eleM_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_u1v1_MP = 0.5*( eleP_nu(1) * eleM_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + 0.5 * eleP_nu(2) * eleM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
    A2_u1v2_MP = 0.5*( 0.5 * eleP_nu(2) * eleM_testPbx' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
    A2_u2v1_MP = 0.5*( 0.5 * eleP_nu(1) * eleM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
    A2_u2v2_MP = 0.5*( 0.5 * eleP_nu(1) * eleM_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + eleP_nu(2) * eleM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_u1v1_MM = 0.5*( eleM_nu(1) * eleM_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPb) ...
        + 0.5 * eleM_nu(2) * eleM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
    A2_u1v2_MM = 0.5*( 0.5 * eleM_nu(2) * eleM_testPbx' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
    A2_u2v1_MM = 0.5*( 0.5 * eleM_nu(1) * eleM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
    A2_u2v2_MM = 0.5*( 0.5 * eleM_nu(1) * eleM_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPb) ...
        + eleM_nu(2) * eleM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2.A2_u1v1 = [A2_u1v1_PP, A2_u1v1_PM;
        A2_u1v1_MP, A2_u1v1_MM];
    A2.A2_u2v1 = [A2_u2v1_PP, A2_u2v1_PM;
        A2_u2v1_MP, A2_u2v1_MM];
    A2.A2_u1v2 = [A2_u1v2_PP, A2_u1v2_PM;
        A2_u1v2_MP, A2_u1v2_MM];
    A2.A2_u2v2 = [A2_u2v2_PP, A2_u2v2_PM;
        A2_u2v2_MP, A2_u2v2_MM];
        
    %-- subPart III, \int_e [trial]\cdot[test]    
    A3_u1v1_PP = eleP_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb);
    A3_u1v1_PM = -eleP_testPb' * bsxfun(@times,phyGaussweights,eleM_trialPb); 
    A3_u1v1_MP = -eleM_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb);
    A3_u1v1_MM = eleM_testPb' * bsxfun(@times,phyGaussweights,eleM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3.A3_u1v1 = [A3_u1v1_PP, A3_u1v1_PM;
        A3_u1v1_MP, A3_u1v1_MM];
    A3.A3_u2v2 = A3.A3_u1v1;
end % function getMatOnEdge_uv


%% the sub function 
function [int_Atrial_Jtest_1, int_Atrial_Jtest_2] = getMatOnInterEdge_qv(phyGaussweights, ...
    eleP_nu, ...
    eleP_trialPb, ...
    eleP_testPb, ...
    eleM_nu, ...
    eleM_trialPb, ...
    eleM_testPb)
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
    A1_pp = 0.5*( eleP_nu(1) * eleP_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A1_pm = 0.5*( eleP_nu(1) * eleP_testPb' * bsxfun(@times,phyGaussweights,eleM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_mp = 0.5*( eleM_nu(1) * eleM_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_mm = 0.5*( eleM_nu(1) * eleM_testPb' * bsxfun(@times,phyGaussweights,eleM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    int_Atrial_Jtest_1 = [A1_pp, A1_pm;
        A1_mp, A1_mm];
        
    
    %-- subPart II, \int_e {(func1,func2)\cdot grad_test}\cdot[trial]
    A2_pp = 0.5*( eleP_nu(2) * eleP_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A2_pm = 0.5*( eleP_nu(2) * eleP_testPb' * bsxfun(@times,phyGaussweights,eleM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_mp = 0.5*( eleM_nu(2) * eleM_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_mm = 0.5*( eleM_nu(2) * eleM_testPb' * bsxfun(@times,phyGaussweights,eleM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
       
    int_Atrial_Jtest_2 = [A2_pp, A2_pm;
        A2_mp, A2_mm];    

end % function getMatOnEdge_qv

%--------------------------------------------------------------------------------------------------------------------%
%----------- end sub function ----------------------------------------------------------------------------------%

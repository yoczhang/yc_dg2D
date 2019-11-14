function [matInterEdge_uv, matInterEdge_qv1, matInterEdge_qv2, matInterEdge_pq] = ...
    InterEdgeAverJumpMDGTimeDarcy(Coeffs_func, meshInfo, option, formulaGauss1D, degreek_u, degreek_q)
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
%   YcZhang 8/27/2017 -- MM/DD/YYYY
%
%   Last modified 8/27/2017 -- MM/DD/YYYY
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
NTbases_q = (degreek_q+1)*(degreek_q+2)/2;

matInterEdge_uv = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
matInterEdge_qv1 = sparse(Nelems*NTbases_u,Nelems*NTbases_q);
matInterEdge_qv2 = sparse(Nelems*NTbases_u,Nelems*NTbases_q);
matInterEdge_pq = sparse(Nelems*NTbases_q,Nelems*NTbases_q);

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
    
    [q_elemP_Pb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_q);
        %> q_eleP_Pb, q_eleP_Pbx, q_eleP_Pby, [Npoints x NTbases_q]
    [q_elemM_Pb, ~, ~] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_q);
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
        %> uv_A1, uv_A2, uv_A3, [2*NTbases_u  x 2*NTbases_u]
        %> uv_A1=[A1_pp, A1_pm; 
        %>            A1_mp, A1_mm];
        %> the same to uv_A2, uv_A3.
    
    matInterEdge_uv = matInterEdge_uv ...
        + sparse(uv_Row(:), uv_Col(:), ...
        - 0*uv_A1(:) + 0*p_epsilon*uv_A2(:) + (p_sigma/areaEdge^p_beta)*uv_A3(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    
    
    %-- subPart 2 -- the AverJump of pressure trial p and velotivy test v.
    % get the position in the big mat
    [qv_Row, qv_Col] = getROWandCOL(elemP, elemM, NTbases_q, NTbases_u);
    [pq_Row, pq_Col] = getROWandCOL(elemP, elemM, NTbases_q, NTbases_q);
    
    % get the matrix
    [qv_A1, qv_A2, pq_A3] = getMatOnInterEdge_qv(phyGweights, ...
        elemP_nu, ...
        q_elemP_Pb, ...
        u_elemP_Pb, ...
        elemM_nu, ...
        q_elemM_Pb, ...
        u_elemM_Pb);
        %> qv_A1, qv_A2, [2*NTbases_u x 2*NTbases_q]
        %> qv_A1=[A1_pp, A1_pm; 
        %>            A1_mp, A1_mm];
        %> the same to qv_A2.
    
    matInterEdge_qv1 = matInterEdge_qv1 ...
        + sparse(qv_Row(:), qv_Col(:), qv_A1(:), Nelems*NTbases_u, Nelems*NTbases_q);
    matInterEdge_qv2 = matInterEdge_qv2 ...
        + sparse(qv_Row(:), qv_Col(:), qv_A2(:), Nelems*NTbases_u, Nelems*NTbases_q);
    matInterEdge_pq = matInterEdge_pq ...
        + sparse(pq_Row(:), pq_Col(:), pq_A3(:), Nelems*NTbases_q, Nelems*NTbases_q);
    
    
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
function [int_Agradtrial_Jtest, int_Agradtest_Jtrial, int_Jtrial_Jtest] = getMatOnInterEdge_uv(phyGaussweights, ...
    valueCoeff_1, valueCoeff_2, ...
    elemP_nu, ...
    elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
    elemP_testPb, elemP_testPbx, elemP_testPby, ...
    elemM_nu, ...
    elemM_trialPb, elemM_trialPbx, elemM_trialPby, ...
    elemM_testPb, elemM_testPbx, elemM_testPby)
%
%   get the Aver{.} and Jump[.] mat on Interior edges.
%
%   int_Agradtest_Jtrial --> \int_e {(func1,func2)\cdot grad_trial}\cdot[test].
%   int_Agradtest_Jtrial --> \int_e {(func1,func2)\cdot grad_test}\cdot[trial].
%   int_Jtrial_Jtest --> \int_e [trial]\cdot[test].
%
%   YcZhang 16/8/2017
%   Last modified 16/8/2017
%

    %-- subPart I, \int_e {(func1,func2)\cdot grad_trial}\cdot[test]
    % {grad_trial}\cdot[test] = 1/2(grad_trial^p + grad_trial^m)\cdot(test^p nu^p + test^m nu^m)
    %AgradtestJtrial_pp = zeros(NTbases_test, NTbases_trial);    
    A1_pp = 0.5*( elemP_nu(1) * elemP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,elemP_trialPbx) ...
        + elemP_nu(2) * elemP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,elemP_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A1_pm = 0.5*( elemP_nu(1) * elemP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,elemM_trialPbx) ...
        + elemP_nu(2) * elemP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,elemM_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_mp = 0.5*( elemM_nu(1) * elemM_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,elemP_trialPbx) ...
        + elemM_nu(2) * elemM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,elemP_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_mm = 0.5*( elemM_nu(1) * elemM_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,elemM_trialPbx) ...
        + elemM_nu(2) * elemM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,elemM_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    int_Agradtrial_Jtest = [A1_pp, A1_pm;
        A1_mp, A1_mm];
        
    
    %-- subPart II, \int_e {(func1,func2)\cdot grad_test}\cdot[trial]
    A2_pp = 0.5*( elemP_nu(1) * elemP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,elemP_trialPb) ...
        + elemP_nu(2) * elemP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,elemP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A2_pm = 0.5*( elemM_nu(1) * elemP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,elemM_trialPb) ...
        + elemM_nu(2) * elemP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,elemM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_mp = 0.5*( elemP_nu(1) * elemM_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,elemP_trialPb) ...
        + elemP_nu(2) * elemM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,elemP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_mm = 0.5*( elemM_nu(1) * elemM_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,elemM_trialPb) ...
        + elemM_nu(2) * elemM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,elemM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
       
    int_Agradtest_Jtrial = [A2_pp, A2_pm;
        A2_mp, A2_mm];    
    
        
    %-- subPart III, \int_e [trial]\cdot[test]    
    A3_pp = elemP_testPb' * bsxfun(@times,phyGaussweights,elemP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A3_pm = -elemP_testPb' * bsxfun(@times,phyGaussweights,elemM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3_mp = -elemM_testPb' * bsxfun(@times,phyGaussweights,elemP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3_mm = elemM_testPb' * bsxfun(@times,phyGaussweights,elemM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    int_Jtrial_Jtest = [A3_pp, A3_pm;
        A3_mp, A3_mm];
end % function getMatOnEdge_uv


%% the sub function 
function [int_Atrial_Jtest_1, int_Atrial_Jtest_2, int_Jtrial_Jtest_3] = getMatOnInterEdge_qv(phyGaussweights, ...
    elemP_nu, ...
    elemP_trialPb, ...
    elemP_testPb, ...
    elemM_nu, ...
    elemM_trialPb, ...
    elemM_testPb)
%
%   get the Aver{.} and Jump[.] mat on Interior edges.
%
%   int_Agradtest_Jtrial --> \int_e {(func1,func2)\cdot grad_trial}\cdot[test].
%   int_Agradtest_Jtrial --> \int_e {(func1,func2)\cdot grad_test}\cdot[trial].
%   int_Jtrial_Jtest --> \int_e [trial]\cdot[test].
%

    %-- subPart I, \int_e {(func1,func2)\cdot grad_trial}\cdot[test]
    % {grad_trial}\cdot[test] = 1/2(grad_trial^p + grad_trial^m)\cdot(test^p nu^p + test^m nu^m)
    %AgradtestJtrial_pp = zeros(NTbases_test, NTbases_trial);    
    A1_pp = 0.5*( elemP_nu(1) * elemP_testPb' * bsxfun(@times,phyGaussweights,elemP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A1_pm = 0.5*( elemP_nu(1) * elemP_testPb' * bsxfun(@times,phyGaussweights,elemM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_mp = 0.5*( elemM_nu(1) * elemM_testPb' * bsxfun(@times,phyGaussweights,elemP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_mm = 0.5*( elemM_nu(1) * elemM_testPb' * bsxfun(@times,phyGaussweights,elemM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    int_Atrial_Jtest_1 = [A1_pp, A1_pm;
        A1_mp, A1_mm];
        
    
    %-- subPart II, \int_e {(func1,func2)\cdot grad_test}\cdot[trial]
    A2_pp = 0.5*( elemP_nu(2) * elemP_testPb' * bsxfun(@times,phyGaussweights,elemP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A2_pm = 0.5*( elemP_nu(2) * elemP_testPb' * bsxfun(@times,phyGaussweights,elemM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_mp = 0.5*( elemM_nu(2) * elemM_testPb' * bsxfun(@times,phyGaussweights,elemP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_mm = 0.5*( elemM_nu(2) * elemM_testPb' * bsxfun(@times,phyGaussweights,elemM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
       
    int_Atrial_Jtest_2 = [A2_pp, A2_pm;
        A2_mp, A2_mm];    
    
    %-- subPart III, \int_e [trial]\cdot[test]    
    A3_pp = elemP_trialPb' * bsxfun(@times,phyGaussweights,elemP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A3_pm = -elemP_trialPb' * bsxfun(@times,phyGaussweights,elemM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3_mp = -elemM_trialPb' * bsxfun(@times,phyGaussweights,elemP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3_mm = elemM_trialPb' * bsxfun(@times,phyGaussweights,elemM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    int_Jtrial_Jtest_3 = [A3_pp, A3_pm;
        A3_mp, A3_mm];

end % function getMatOnEdge_qv

%--------------------------------------------------------------------------------------------------------------------%
%----------- end sub function ----------------------------------------------------------------------------------%

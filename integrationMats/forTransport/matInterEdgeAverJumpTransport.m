function [matInterEdge, matInterEdge_upwind_term] = ...
    matInterEdgeAverJumpTransport(Coeffs_func, meshInfo, option,formulaGauss1D,degreek_trial,degreek_test)
%
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
%       matInterEdge, [Nelems*NTbases_test x NTelems*NTbases_trial], 
%           \sum_e \int_e {(func1,func2)*grad u}\cdot[v] +  \epsilon\int_e {(func1,func2)*grad v}\cdot[u] + penalty \int_e [u]\cdot[v]. 
%       
%
%   YcZhang 5/8/2017
%
%   Last modified 16/8/2017
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
vector_u_1 = Coeffs_func{3};
vector_u_2 = Coeffs_func{4};

% setting
NTbases_trial = (degreek_trial+1)*(degreek_trial+2)/2;
NTbases_test = (degreek_test+1)*(degreek_test+2)/2;

matInterEdge = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);
matInterEdge_upwind_term = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);


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
    
    ele1 = meshInfo.edge2elem(interEindx,1);
    ele2 = meshInfo.edge2elem(interEindx,2);
    local_e1 = meshInfo.edge2elem(interEindx,3);
    local_e2 = meshInfo.edge2elem(interEindx,4);
        %> [n1, n2]=edge2elem(k,1:2) stands for the n1-th elem and n2-th elem share the k-th edge.
        %> local_e1=edge2elem(k,3) stands for the local edge index in n1-th elem of the k-th edge .
        %> also, local_e2=edge2elem(k,4) stands for the local edge index in n2-th elem of the k-th edge .
    eleP = min(ele1, ele2); % let the smaller number of element as the element +.
    eleM = max(ele1, ele2); % let the larger number of element as the element -.
    
    if eleP == ele1
        eIndx_eleP = local_e1;
        eIndx_eleM = local_e2;
    else
        eIndx_eleP = local_e2;
        eIndx_eleM = local_e1;
    end 
    
    eleP_nu = meshInfo.nuEdge0Elem{eleP}(:,eIndx_eleP); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleP.
    eleM_nu = meshInfo.nuEdge0Elem{eleM}(:,eIndx_eleM); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleM.
    
    eleP_xT = meshInfo.centroidElem(eleP,1);  
    eleP_yT = meshInfo.centroidElem(eleP,2);
    eleP_hT = meshInfo.hElem(eleP);
        %> the centroid(xing xin) and diameter of eleP.
    eleM_xT = meshInfo.centroidElem(eleM,1);  
    eleM_yT = meshInfo.centroidElem(eleM,2);
    eleM_hT = meshInfo.hElem(eleM);
        %> the centroid(xing xin) and diameter of eleM.
        
    %-- get the value of Coeffs func at Gauss Points
    valueCoeff_k11 = Coeff_k11(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valueCoeff_k22 = Coeff_k22(phyGpointsX, phyGpointsY); % [Npoints x 1]
    
    % setting elem_flag and get the value of vector_u
    elem_flag = settingElemFlagTransport(ii); % here 
    value_vector_u_1 = vector_u_1(phyGpointsX, phyGpointsY, elem_flag); 
        %> [Npoints x 1], this value_vector_u is a vector-valued.
    value_vector_u_2 = vector_u_2(phyGpointsX, phyGpointsY, elem_flag); 
        %> [Npoints x 1], this value_vector_u is a vector-valued.
    
    %-- subPart -- get the bases on quad
    [eleP_trialPb, eleP_trialPbx, eleP_trialPby] = ...
        localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, degreek_trial);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
    [eleM_trialPb, eleM_trialPbx, eleM_trialPby] = ...
        localBases2D(eleM_xT, eleM_yT, eleM_hT, phyGpointsX, phyGpointsY, degreek_trial);
        %> eleMtrialPb, eleMtrialPbx, eleMtrialPby, [Npoints x NTg1bases_trial]
    
    [eleP_testPb, eleP_testPbx, eleP_testPby] = ...
        localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, degreek_test);
        %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]
    [eleM_testPb, eleM_testPbx, eleM_testPby] = ...
        localBases2D(eleM_xT, eleM_yT, eleM_hT, phyGpointsX, phyGpointsY, degreek_test);
        %> eleMtestPb, eleMtestPbx, eleMtestPby, [Npoints x NTg1bases_test]
        
    %<<-- End Part I ---------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
    
    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %
    %-- subPart I -- the AverJump of velotivy trial u and velotivy test v 
    % get the position in the big mat    
    [Row, Col] = getROWandCOL(eleP, eleM, NTbases_trial, NTbases_test);
    
    % get the matrix
    [uv_A1, uv_A2, uv_A3] = getMatOnInterEdge(phyGweights, ...
        valueCoeff_k11, valueCoeff_k22, ...
        eleP_nu, ...
        eleP_trialPb, eleP_trialPbx, eleP_trialPby, ...
        eleP_testPb, eleP_testPbx, eleP_testPby, ...
        eleM_nu, ...
        eleM_trialPb, eleM_trialPbx, eleM_trialPby, ...
        eleM_testPb, eleM_testPbx, eleM_testPby);
        %> uv_A1, uv_A2, uv_A3, [NTbases_test x NTbases_trial]
        %> uv_A1=[A1_pp, A1_pm; 
        %>            A1_mp, A1_mm];
        %> the same to uv_A2, uv_A3.
        
    upwind_term = getMatOnInterEdge_upwindterm(phyGweights, ...
        value_vector_u_1, ...
        value_vector_u_2, ...
        eleP_nu, ...
        eleP_trialPb, ...
        eleP_testPb, ...
        eleM_nu, ...
        eleM_trialPb, ...
        eleM_testPb);
    
    matInterEdge = matInterEdge ...
        + sparse(Row(:), Col(:), ...
        - uv_A1(:) + p_epsilon*uv_A2(:) + (p_sigma/areaEdge^p_beta)*uv_A3(:), ...
        Nelems*NTbases_test, Nelems*NTbases_trial);
    
    matInterEdge_upwind_term = matInterEdge_upwind_term ...
        + sparse(Row(:), Col(:), ...
        upwind_term(:), ...
        Nelems*NTbases_test, Nelems*NTbases_trial);
        
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
function [int_Agradtrial_Jtest, int_Agradtest_Jtrial, int_Jtrial_Jtest] = getMatOnInterEdge(phyGaussweights, ...
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
    %-- subPart I, \int_e {(func1,func2)\cdot grad_trial}\cdot[test]
    % {grad_trial}\cdot[test] = 1/2(grad_trial^p + grad_trial^m)\cdot(test^p nu^p + test^m nu^m)
    %AgradtestJtrial_pp = zeros(NTbases_test, NTbases_trial);    
    A1_pp = 0.5*( eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A1_pm = 0.5*( eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPbx) ...
        + eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_mp = 0.5*( eleM_nu(1) * eleM_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + eleM_nu(2) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_mm = 0.5*( eleM_nu(1) * eleM_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPbx) ...
        + eleM_nu(2) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    int_Agradtrial_Jtest = [A1_pp, A1_pm;
        A1_mp, A1_mm];
        
    
    %-- subPart II, \int_e {(func1,func2)\cdot grad_test}\cdot[trial]
    A2_pp = 0.5*( eleP_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + eleP_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A2_pm = 0.5*( eleM_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPb) ...
        + eleM_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_mp = 0.5*( eleP_nu(1) * eleM_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + eleP_nu(2) * eleM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_mm = 0.5*( eleM_nu(1) * eleM_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPb) ...
        + eleM_nu(2) * eleM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
       
    int_Agradtest_Jtrial = [A2_pp, A2_pm;
        A2_mp, A2_mm];    
    
        
    %-- subPart III, \int_e [trial]\cdot[test]    
    A3_pp = eleP_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A3_pm = -eleP_testPb' * bsxfun(@times,phyGaussweights,eleM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3_mp = -eleM_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3_mm = eleM_testPb' * bsxfun(@times,phyGaussweights,eleM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    int_Jtrial_Jtest = [A3_pp, A3_pm;
        A3_mp, A3_mm];
end % function getMatOnEdge


%% the sub function 
function upwind_term = getMatOnInterEdge_upwindterm(phyGaussweights, ...
    vector_u_1, ...
    vector_u_2, ...
    eleP_nu, ...
    eleP_trialPb, ...
    eleP_testPb, ...
    eleM_nu, ...
    eleM_trialPb, ...
    eleM_testPb)
%
%   get the vector_u\cdot n c^{up}testv on interior edges using upwind scheme.
%
%
    NTbases_trial = size(eleP_trialPb,2);
    NTbases_test = size(eleP_testPb,2);
    A_pp = zeros(NTbases_test,NTbases_trial);
    A_pm = zeros(NTbases_test,NTbases_trial);
    A_mp = zeros(NTbases_test,NTbases_trial);
    A_mm = zeros(NTbases_test,NTbases_trial);
    
% %     if (vector_u(1)*eleP_nu(1)+vector_u(2)*eleP_nu(2)) >= 0
% %         A_pp = eleP_nu(1) * eleP_testPb' * bsxfun(@times,vector_u_1.*phyGaussweights,eleP_trialPb) ...
% %             + eleP_nu(2) * eleP_testPb' * bsxfun(@times,vector_u_2.*phyGaussweights,eleP_trialPb);
% %     else
% %         A_pm = eleP_nu(1) * eleP_testPb' * bsxfun(@times,vector_u_1.*phyGaussweights,eleM_trialPb) ...
% %             + eleP_nu(2) * eleP_testPb' * bsxfun(@times,vector_u_2.*phyGaussweights,eleM_trialPb);
% %     end
% %     if (vector_u(1)*eleM_nu(1)+vector_u(2)*eleM_nu(2)) >= 0
% %         A_mm = eleM_nu(1) * eleM_testPb' * bsxfun(@times,vector_u_1.*phyGaussweights,eleM_trialPb) ...
% %             + eleM_nu(2) * eleM_testPb' * bsxfun(@times,vector_u_2.*phyGaussweights,eleM_trialPb);
% %     else
% %         A_mp = eleM_nu(1) * eleM_testPb' * bsxfun(@times,vector_u_1.*phyGaussweights,eleP_trialPb) ...
% %             + eleM_nu(2) * eleM_testPb' * bsxfun(@times,vector_u_2.*phyGaussweights,eleP_trialPb);
% %     end
    
    
    if (vector_u_1(1)*eleP_nu(1)+vector_u_2(1)*eleP_nu(2)) >= 0
        A_pp = eleP_nu(1) * eleP_testPb' * bsxfun(@times,vector_u_1.*phyGaussweights,eleP_trialPb) ...
            + eleP_nu(2) * eleP_testPb' * bsxfun(@times,vector_u_2.*phyGaussweights,eleP_trialPb);
        A_mp = eleM_nu(1) * eleM_testPb' * bsxfun(@times,vector_u_1.*phyGaussweights,eleP_trialPb) ...
            + eleM_nu(2) * eleM_testPb' * bsxfun(@times,vector_u_2.*phyGaussweights,eleP_trialPb);
    else
        A_pm = eleP_nu(1) * eleP_testPb' * bsxfun(@times,vector_u_1.*phyGaussweights,eleM_trialPb) ...
            + eleP_nu(2) * eleP_testPb' * bsxfun(@times,vector_u_2.*phyGaussweights,eleM_trialPb);
        A_mm = eleM_nu(1) * eleM_testPb' * bsxfun(@times,vector_u_1.*phyGaussweights,eleM_trialPb) ...
            + eleM_nu(2) * eleM_testPb' * bsxfun(@times,vector_u_2.*phyGaussweights,eleM_trialPb);
    end

    upwind_term = [A_pp, A_pm;
        A_mp, A_mm];
    
end % function getMatOnEdge

%--------------------------------------------------------------------------------------------------------------------%
%----------- end sub function ----------------------------------------------------------------------------------%

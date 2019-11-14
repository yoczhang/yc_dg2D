function [D_matDirichletEdge, C_matDirichletEdge] = ...
    matDirichletEdgeAverJumpTransport(Coeffs_func, meshInfo, option, formulaGauss1D, degreek_trial, degreek_test)
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
%       D_matDirichletEdge, diffusion term Dirichlet edges, [Nelems*NTbases_test x NTelems*NTbases_trial], 
%           \sum_e \int_e {(func1,func2)*grad u}\cdot[v] +  \epsilon\int_e {(func1,func2)*grad v}\cdot[u] + penalty \int_e [u]\cdot[v]. 
%       C_matDirichletEdge, convection term Dirichlet edges, [Nelems*NTbases_test x NTelems*NTbases_trial], 
%  
%
%   YcZhang 4/9/2017
%
%   Last modified 4/9/2017
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
vector_u_1 = Coeffs_func{3};
vector_u_2 = Coeffs_func{4};

% setting
NTbases_trial = (degreek_trial+1)*(degreek_trial+2)/2;
NTbases_test = (degreek_test+1)*(degreek_test+2)/2;

D_matDirichletEdge = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);
C_matDirichletEdge = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);


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
    phyGweights = areaEdge * formulaGauss1D(:,2);
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

    %-- setting elem_flag and get the value of vector_u
    elem_flag = settingElemFlagTransport(ii); % here 
    value_vector_u_1 = vector_u_1(phyGpointsX, phyGpointsY, elem_flag); 
        %> [Npoints x 1], this value_vector_u is a vector-valued.
    value_vector_u_2 = vector_u_2(phyGpointsX, phyGpointsY, elem_flag); 
        %> [Npoints x 1], this value_vector_u is a vector-valued.
    
    
    %-- subPart -- get the bases on quad
    [eleP_trialPb, eleP_trialPbx, eleP_trialPby] = ...
        localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, degreek_trial);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
    
    [eleP_testPb, eleP_testPbx, eleP_testPby] = ...
        localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, degreek_test);
        %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]
        
    %<<-- End Part I ---------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
    
    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %
    %-- subPart 1 -- the interior edges
    % get the position in the big mat    
    [Row, Col] = getROWandCOL(eleP, 0, NTbases_trial, NTbases_test);
    
    % get the matrix
    [uv_A1, uv_A2, uv_A3] = getMatOnDirichletEdge(phyGweights, ...
        valueCoeff_k11, valueCoeff_k22, ...
        eleP_nu, ...
        eleP_trialPb, eleP_trialPbx, eleP_trialPby, ...
        eleP_testPb, eleP_testPbx, eleP_testPby);
        %> uv_A1, uv_A2, uv_A3, [NTbases_test x NTbases_trial]
        %> uv_A1=[A1_pp];
        %> the same to uv_A2, uv_A3.
    
    D_matDirichletEdge = D_matDirichletEdge ...
        + sparse(Row(:), Col(:), ...
        - uv_A1(:) + p_epsilon*uv_A2(:) + (p_sigma/areaEdge^p_beta)*uv_A3(:), ...
        Nelems*NTbases_test, Nelems*NTbases_trial);
    
    C_mat_temp = eleP_nu(1) * eleP_testPb' * bsxfun(@times,value_vector_u_1.*phyGweights,eleP_trialPb) ...
        + eleP_nu(2) * eleP_testPb' * bsxfun(@times,value_vector_u_2.*phyGweights,eleP_trialPb);
    C_matDirichletEdge = C_matDirichletEdge ...
        + sparse(Row(:), Col(:), ...
        C_mat_temp(:), ...
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
function [int_Agradtrial_Jtest, int_Agradtest_Jtrial, int_Jtrial_Jtest] = getMatOnDirichletEdge(phyGaussweights, ...
    valueCoeff_1, valueCoeff_2, ...
    eleP_nu, ...
    eleP_trialPb, eleP_trialPbx, eleP_trialPby, ...
    eleP_testPb, eleP_testPbx, eleP_testPby)
%
%   get the Aver{.} and Jump[.] mat on Dirichlet edges.
%
%   int_Agradtest_Jtrial --> \int_e {(func1,func2)\cdot grad_trial}\cdot[test].
%   int_Agradtest_Jtrial --> \int_e {(func1,func2)\cdot grad_test}\cdot[trial].
%   int_Jtrial_Jtest --> \int_e [trial]\cdot[test].
%
    %-- subPart I, \int_e {(func1,func2)\cdot grad_trial}\cdot[test]
    % {grad_trial}\cdot[test] = 1/2(grad_trial^p + grad_trial^m)\cdot(test^p nu^p + test^m nu^m)
    %AgradtestJtrial_pp = zeros(NTbases_test, NTbases_trial);    
    A1_pp = eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    int_Agradtrial_Jtest = A1_pp;
        
    
    %-- subPart II, \int_e {(func1,func2)\cdot grad_test}\cdot[trial]
    A2_pp = eleP_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + eleP_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    int_Agradtest_Jtrial = A2_pp;    
    
        
    %-- subPart III, \int_e [trial]\cdot[test]    
    A3_pp = eleP_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].   
    int_Jtrial_Jtest = A3_pp;
end % function getMatOnBoundEdge

%--------------------------------------------------------------------------------------------------------------------%
%----------- end sub function ----------------------------------------------------------------------------------%

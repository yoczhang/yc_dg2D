function [matInterEdge, matDiffusiveAdapEdge, matConvectionUpwindEdge] = ...
    matInterEdgeAverJumpAdapTransport_noAdap1(Coeffs_func, meshInfo, option,formulaGauss1D,degreek_trial,degreek_test)
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
%   YcZhang 23/9/2017
%
%   Last modified 23/9/2017
%
%

% mesh information: interior edges 
Nelems = meshInfo.Nelems;
interEdges = meshInfo.interEdgeIndex; % here we need the interEdges information.

p_epsilon = option.p_epsilon;
% p_epsilon1 = option.p_epsilon1; % this is using to the adaptive-transport
p_sigma = option.p_sigma;
% p_sigma1 = option.p_sigma1; % this is using to the adaptive-transport
p_beta = option.p_beta;

vector_u_1 = Coeffs_func{3};
vector_u_2 = Coeffs_func{4};

% setting
NTbases_trial = (degreek_trial+1)*(degreek_trial+2)/2;
NTbases_test = (degreek_test+1)*(degreek_test+2)/2;

matInterEdge = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);
matConvectionUpwindEdge = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);
matDiffusiveAdapEdge = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);

%-------------- for test -------------------------------
H_DC = min( meshInfo.DiffusivityCoeffs(:));
%---------------------------------------------------------
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
        
    %-- get the Diffusivity coefficients of elemP and elemM
    elemP_DC = meshInfo.DiffusivityCoeffs(elemP);
    elemM_DC = meshInfo.DiffusivityCoeffs(elemM);
    
    % setting elem_flag and get the value of vector_u
    elemP_flag = settingElemFlagAdapTransport(meshInfo,elemP); % here 
    valueP_vector_u_1 = vector_u_1(phyGpointsX, phyGpointsY, elemP_flag); 
        %> [Npoints x 1], this value_vector_u is a vector-valued.
    valueP_vector_u_2 = vector_u_2(phyGpointsX, phyGpointsY, elemP_flag); 
        %> [Npoints x 1], this value_vector_u is a vector-valued.
        
    elemM_flag = settingElemFlagAdapTransport(meshInfo,elemM); % here 
    valueM_vector_u_1 = vector_u_1(phyGpointsX, phyGpointsY, elemM_flag); 
        %> [Npoints x 1], this value_vector_u is a vector-valued.
    valueM_vector_u_2 = vector_u_2(phyGpointsX, phyGpointsY, elemM_flag); 
        %> [Npoints x 1], this value_vector_u is a vector-valued.
        
    if sum(abs(valueP_vector_u_1-valueM_vector_u_1))+sum(abs(valueP_vector_u_2-valueM_vector_u_2))>0
        error('in matInterEdgeAverJumpAdapTransport_noAdap1.m, valueP_vector_u ~= valueM_vector_u')
    end
    
    %-- subPart -- get the bases on quad
    [elemP_trialPb, elemP_trialPbx, elemP_trialPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_trial);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
    [elemM_trialPb, elemM_trialPbx, elemM_trialPby] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_trial);
        %> eleMtrialPb, eleMtrialPbx, eleMtrialPby, [Npoints x NTg1bases_trial]
    
    [elemP_testPb, elemP_testPbx, elemP_testPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_test);
        %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]
    [elemM_testPb, elemM_testPbx, elemM_testPby] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_test);
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
    [Row, Col] = getROWandCOL(elemP, elemM, NTbases_trial, NTbases_test);

    % get the matrix on all the interior edges, 
    elemP_value_DC = elemP_DC * ones(length(phyGweights),1);
	elemM_value_DC = elemM_DC * ones(length(phyGweights),1);
    [uv_A1, uv_A2, uv_A3] = ...
        getMatOnInterEdge(phyGweights, ...
        elemP_value_DC, elemM_value_DC, ...
        elemP_nu, ...
        elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
        elemP_testPb, elemP_testPbx, elemP_testPby, ...
        elemM_nu, ...
        elemM_trialPb, elemM_trialPbx, elemM_trialPby, ...
        elemM_testPb, elemM_testPbx, elemM_testPby);
            %> uv_A1, uv_A2, uv_A3, [NTbases_test x NTbases_trial]
            %> uv_A1=[A1_pp, A1_pm; 
            %>            A1_mp, A1_mm];
            %> the same to uv_A2, uv_A3.
       
%     %-------------- for test -------------------------------
%     if elemP_DC == H_DC && elemM_DC == H_DC
%         matInterEdge = matInterEdge ...
%             + sparse(Row(:), Col(:), ...
%             - uv_A1(:) + p_epsilon*uv_A2(:) + (0.0/areaEdge^p_beta)*uv_A3(:), ...
%             Nelems*NTbases_test, Nelems*NTbases_trial);
%     else
%         matInterEdge = matInterEdge ...
%             + sparse(Row(:), Col(:), ...
%             - uv_A1(:) + p_epsilon*uv_A2(:) + (p_sigma/areaEdge^p_beta)*uv_A3(:), ...
%             Nelems*NTbases_test, Nelems*NTbases_trial);
%     end
%     %--------------------------------------------------------
    
    matInterEdge = matInterEdge ...
        + sparse(Row(:), Col(:), ...
        - uv_A1(:) + p_epsilon*uv_A2(:) + (p_sigma/areaEdge^p_beta)*uv_A3(:), ...
        Nelems*NTbases_test, Nelems*NTbases_trial);

    ConvectionUpwind_mat = getMatOnConvectionUpwindEdge(phyGweights, ...
        valueP_vector_u_1, ...
        valueP_vector_u_2, ...
        elemP_nu, ...
        elemP_trialPb, ...
        elemP_testPb, ...
        elemM_nu, ...
        elemM_trialPb, ...
        elemM_testPb);
    
    matConvectionUpwindEdge = matConvectionUpwindEdge ...
        + sparse(Row(:), Col(:), ...
        ConvectionUpwind_mat(:), ...
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
function [Adap_A1_normal, Adap_A2_normal, Adap_A3_normal] = ...
	getMatOnInterEdge(phyGaussweights, ...
	elemP_value_DC, elemM_value_DC, ...
	elemP_nu, ...
	elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
	elemP_testPb, elemP_testPbx, elemP_testPby, ...
	elemM_nu, ...
	elemM_trialPb, elemM_trialPbx, elemM_trialPby, ...
	elemM_testPb, elemM_testPbx, elemM_testPby)
%
%   %---------------------------------------------------------
%       get the mat on Diffusive adaptive edges, using the NORMAL scheme.
%   %---------------------------------------------------------
%
%   input:
%       elemP_value_DC, elemM_value_DC, [Npoints x 1], the Diffusivity
%       coefficients on elemP and elemM.
%
%
%
%   YcZhang 23/9/2017
%
%   Last modified 23/9/2017
%
%

A1_pp = 0.5*( elemP_nu(1) * elemP_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemP_trialPbx) ...
    + elemP_nu(2) * elemP_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemP_trialPby) );
    %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].

A1_pm = 0.5*( elemP_nu(1) * elemP_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemM_trialPbx) ...
    + elemP_nu(2) * elemP_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemM_trialPby) );
    %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].

A1_mp = 0.5*( elemM_nu(1) * elemM_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemP_trialPbx) ...
    + elemM_nu(2) * elemM_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemP_trialPby) );
    %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].

A1_mm = 0.5*( elemM_nu(1) * elemM_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemM_trialPbx) ...
    + elemM_nu(2) * elemM_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemM_trialPby) );
    %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].

Adap_A1_normal = [A1_pp, A1_pm;
    A1_mp, A1_mm];


%-- subPart II, \int_e {(func1,func2)\cdot grad_test}\cdot[trial]
A2_pp = 0.5*( elemP_nu(1) * (elemP_value_DC.*elemP_testPbx)' * bsxfun(@times,phyGaussweights,elemP_trialPb) ...
    + elemP_nu(2) * (elemP_value_DC.*elemP_testPby)' * bsxfun(@times,phyGaussweights,elemP_trialPb) );
    %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].

A2_pm = 0.5*( elemM_nu(1) * (elemP_value_DC.*elemP_testPbx)' * bsxfun(@times,phyGaussweights,elemM_trialPb) ...
    + elemM_nu(2) * (elemP_value_DC.*elemP_testPby)' * bsxfun(@times,phyGaussweights,elemM_trialPb) );
    %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].

A2_mp = 0.5*( elemP_nu(1) * (elemM_value_DC.*elemM_testPbx)' * bsxfun(@times,phyGaussweights,elemP_trialPb) ...
    + elemP_nu(2) * (elemM_value_DC.*elemM_testPby)' * bsxfun(@times,phyGaussweights,elemP_trialPb) );
    %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].

A2_mm = 0.5*( elemM_nu(1) * (elemM_value_DC.*elemM_testPbx)' * bsxfun(@times,phyGaussweights,elemM_trialPb) ...
    + elemM_nu(2) * (elemM_value_DC.*elemM_testPby)' * bsxfun(@times,phyGaussweights,elemM_trialPb) );
    %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
   
Adap_A2_normal = [A2_pp, A2_pm;
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

Adap_A3_normal = [A3_pp, A3_pm;
    A3_mp, A3_mm];
        
        
end % function getMatOnDiffusiveUpwindEdge_normal(.)


%% the sub function 
function upwind_term = getMatOnConvectionUpwindEdge(phyGaussweights, ...
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
    % first we need to allocate the matrix, because in the last some matrix
    % may not be assigned values.
    NTbases_trial = size(eleP_trialPb,2);
    NTbases_test = size(eleP_testPb,2);
    A_pp = zeros(NTbases_test,NTbases_trial);
    A_pm = zeros(NTbases_test,NTbases_trial);
    A_mp = zeros(NTbases_test,NTbases_trial);
    A_mm = zeros(NTbases_test,NTbases_trial);
    
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

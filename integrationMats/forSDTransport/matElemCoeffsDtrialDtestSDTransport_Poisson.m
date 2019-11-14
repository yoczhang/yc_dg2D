function [matElemDuDv, vecRhsF] = ...
    matElemCoeffsDtrialDtestSDTransport_Poisson(Coeffs_func, func_f, meshInfo, formulaGauss2D, degreek_trial, degreek_test)
%
%   %---------------------------------------------------------------------
%       Just copy from: matElemCoeffsDtrialDtestPoisson.m
%   %----------------------------------------------------------------------
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
%                   Coeffs{1} is the (\int_T Coeffs{1}uv dx),
%                   Coeffs{2}, and Coeffs{3} are the K=[k11, k12; k21, k22], 
%                   let k12=0, k21=0, the Coeffs{2}=k11, Coeffs{3}= k22. 
%                   Coeffs{2} is the (\int_T Coeffs{2} u_x v dx). 
%                   Coeffs{3} is the (\int_T Coeffs{3} u_y v dx).
%       rhs_f_func, the cell-type, here
%                   rhs_f_func{i} is the first rhs f func.
%       meshInfos, the mesh information.
%       formulaGauss2D, the 1d Gauss quadrature formula, size: a matrix, [Npoints x 3],
%               the first column is the x-coordinates of all Gauss-Points,
%               the second column is the y-coordinates of all Gauss-Points,
%               the third is the weights of all Gauss-Points.
%       g1Trial_k, the group1 trial polynomial degree k.
%       g1Test_k, the group1 test polynomial degree k.
%           %> here we consider that, in the Stokes equation, there will have VELOCITY u and PRESSURE p, 
%           %> so there have \int_e uv ds, then g1Trial_k=u_degree_k, g1Test_k=u_degree_k. 
%           %> And also have \int_e pv ds, then g2Trial_k=p_degree_k, g2Test_k=u_degree_k. 
%
%   output:
%       matElemDuDv, [Nelems*NTg1bases_test x NTelems*NTg1bases_trial], \sum_e \int_e {grad u}\cdot[v]. 
%       rhs_fh, [Nelems*NTg1bases_test x 1].
%
%
%   YcZhang 11/8/2017
%
%   Last modified 16/8/2017
%
%

Nelems = meshInfo.Nelems;

Coeff_k11 = Coeffs_func{1};
Coeff_k22 = Coeffs_func{2};
Coeff_1 = Coeffs_func{3};


% group1 setting
NTbases_trial = (degreek_trial+1)*(degreek_trial+2)/2;
NTbases_test = (degreek_test+1)*(degreek_test+2)/2;

matElemDuDv = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);
vecRhsF = zeros(Nelems*NTbases_test,1);


% the simple build-in function to get the mat [NTg1bases_test x NTg1bases_trial]
getMatOnElem = @(funcValue, trialValue, testValue, phyGweights) ...
    testValue' * bsxfun(@times, phyGweights.*funcValue, trialValue);
    %> input:
    %>      funcValue, [Npoints x 1], the value the coeffients function at Gauss points.
    %>      trialValue, [Npoints x NTbases_trial], the trial bases values at Gauss points. 
    %>      testValue, [Npoints x NTbases_test], the test bases values at Gauss points. 
    %
    %> output:
    %>       [NTbases_test x NTbases_trial], the mat of (\int_T uv  dx) or
    %>                                                       (\int_T u_x v_y  dx) or ...

for ii = 1:Nelems
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    singleElem = meshInfo.elem{ii,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges on singleElem.
    
    concavePointElem_ii = meshInfo.concavePointElem(ii); % [1 x 1], the concavePointElem of ii-the elem.
    coordv = meshInfo.node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    centroidElem = meshInfo.centroidElem(ii,:); % [1 x 2], the centroid(xing xin) (x-coord, y-coord) of element.
    
    coordTri0Elem = getcoordTri0Elem(singleNE, concavePointElem_ii, centroidElem, coordv); 
        %> if singleNE==3, coordTri0Elem, [3 x 2]. else coordTri0Elem, [3*singleNE x 2].
    
    ele_xT = centroidElem(1);
    ele_yT = centroidElem(2); 
    ele_hT = meshInfo.hElem(ii); 
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
    %-- get the Diffusivity coefficients of elemP and elemM
    elem_DC = meshInfo.DiffusivityCoeffs(ii);
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %-- the mat and vec setting
    TrialTest_temp = zeros(NTbases_test, NTbases_trial); % acutally, in the Poisson, this var will not be used.
    DxTrialDxTest_temp = zeros(NTbases_test, NTbases_trial);
    DyTrialDyTest_temp = zeros(NTbases_test, NTbases_trial);
    rhs_f_temp = zeros(NTbases_test,1);
    
    %-- the integration on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D);

        % get the bases values on quad points
        [trialPb, trialPbx, trialPby] = localBases2D(ele_xT, ele_yT, ele_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_trial);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases_trial].
        [testPb, testPbx, testPby] = localBases2D(ele_xT, ele_yT, ele_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_test);
            %> testPb, testPbx, testPby, [Npoints x NTbases_test].
        
        %--- the funcValue may be chosen by case.
        valueCoeff_1 = Coeff_1(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k11 = Coeff_k11(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k22 = Coeff_k22(phyGpoints(:,1), phyGpoints(:,2));
        valuefunc_f = func_f(phyGpoints(:,1), phyGpoints(:,2));
        elem_value_DC = elem_DC * ones(length(phyGweights),1);
        
        %-- MatElem
        TrialTest_temp = TrialTest_temp + getMatOnElem(elem_value_DC.*valueCoeff_1, trialPb, testPb, phyGweights);
        DxTrialDxTest_temp = DxTrialDxTest_temp + getMatOnElem(elem_value_DC.*valueCoeff_k11, trialPbx, testPbx, phyGweights);
        DyTrialDyTest_temp = DyTrialDyTest_temp + getMatOnElem(elem_value_DC.*valueCoeff_k22, trialPby, testPby, phyGweights);
        
        %-- VecF   
        rhs_f_temp = rhs_f_temp + testPb' * ( phyGweights.*valuefunc_f );
    end % for nt
    
    %-- add MatElem to the position in the big mat
    Row_elem = ((ii-1)*NTbases_test+1 : ii*NTbases_test)' ...
        * ones(1,NTbases_trial); % [NTg1bases_test x NTg1bases_trial]
    Col_elem = ones(NTbases_test, 1) ...
        *((ii-1)*NTbases_trial+1 : ii*NTbases_trial); % [NTg1bases_test x NTg1bases_trial]
    
    matElemDuDv = matElemDuDv ...
        + sparse(Row_elem(:), Col_elem(:), ...
        DxTrialDxTest_temp(:) + DyTrialDyTest_temp(:), ...
        Nelems*NTbases_test, Nelems*NTbases_trial);
    
    %-- value_rhs_f -- add to the position in the big vector
    Row_vec = ((ii-1)*NTbases_test+1 : ii*NTbases_test);
    vecRhsF(Row_vec,1) = vecRhsF(Row_vec,1) + rhs_f_temp;
    
    %<<-- End Part II --------------------------------------------------------------------------------------
    
end % for ii

end % function matElemCoeffsTrialTest





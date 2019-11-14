function [matElem_uv, matElem_DuDv, matElem_uvx, vecRhsF] = ...
    matElemCoeffsDtrialDtestAdapTransport_fortestSDTransport1(Coeffs_func, func_f, meshInfo, formulaGauss2D, degreek_trial, degreek_test)

%   YcZhang 11/8/2017
%
%   Last modified 16/8/2017
%
%

twoDomain2oneDomain_elem = ...
    [1, 1;
    2, 2;
    3, 5;
    4, 6;
    5, 3;
    6, 4;
    7, 7;
    8, 8];

Nelems = meshInfo.Nelems;

Coeff_k11 = Coeffs_func{1};
Coeff_k22 = Coeffs_func{2};
Coeff_1 = @(x,y) 1+0.*x;
vector_u_1 = Coeffs_func{3};
vector_u_2 = Coeffs_func{4};

% group1 setting
NTbases_trial = (degreek_trial+1)*(degreek_trial+2)/2;
NTbases_test = (degreek_test+1)*(degreek_test+2)/2;

matElem_DuDv = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);
matElem_uv = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);
matElem_uvx = sparse(Nelems*NTbases_test,Nelems*NTbases_trial);
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

for ii_twoDomain = 1:Nelems
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    
    %-------------- for test --------------------------------------
    % using twoDomain2oneDomain_elem to change the elemIndex.
    ii = twoDomain2oneDomain_elem(ii_twoDomain,2);
    %---------------------------------------------------------------
    
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
    
    %-- get the Diffusivity coefficients of elemP and elemM
    elem_DC = meshInfo.DiffusivityCoeffs(ii);
    
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
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
    TrialDxTest_temp = zeros(NTbases_test, NTbases_trial);
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
        
        % setting elem_flag and get the value of vector_u
        elem_flag = settingElemFlagAdapTransport(meshInfo,ii); % here 
        value_vector_u_1 = vector_u_1(phyGpoints(:,1), phyGpoints(:,2), elem_flag); % this value_vector_u is a vector-valued.
        value_vector_u_2 = vector_u_2(phyGpoints(:,1), phyGpoints(:,2), elem_flag); % this value_vector_u is a vector-valued.
        
        %--- the funcValue may be chosen by case.
        valueCoeff_1 = Coeff_1(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k11 = Coeff_k11(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k22 = Coeff_k22(phyGpoints(:,1), phyGpoints(:,2));
        valuefunc_f = func_f(0,phyGpoints(:,1), phyGpoints(:,2),elem_flag);
        elem_value_DC = elem_DC * ones(length(phyGweights),1);
        
        
        %-- MatElem
        TrialTest_temp = TrialTest_temp + getMatOnElem(valueCoeff_1, trialPb, testPb, phyGweights);
        DxTrialDxTest_temp = DxTrialDxTest_temp + getMatOnElem(elem_value_DC.*valueCoeff_k11, trialPbx, testPbx, phyGweights);
        DyTrialDyTest_temp = DyTrialDyTest_temp + getMatOnElem(elem_value_DC.*valueCoeff_k22, trialPby, testPby, phyGweights);
        TrialDxTest_temp = TrialDxTest_temp + getMatOnElem(value_vector_u_1, trialPb, testPbx, phyGweights) ...
            + getMatOnElem(value_vector_u_2, trialPb, testPby, phyGweights);
        
        %-- VecF   
        rhs_f_temp = rhs_f_temp + testPb' * ( phyGweights.*valuefunc_f );
    end % for nt
    
    %-- add MatElem to the position in the big mat
    Row_elem = ((ii_twoDomain-1)*NTbases_test+1 : ii_twoDomain*NTbases_test)' ...
        * ones(1,NTbases_trial); % [NTg1bases_test x NTg1bases_trial]
    Col_elem = ones(NTbases_test, 1) ...
        *((ii_twoDomain-1)*NTbases_trial+1 : ii_twoDomain*NTbases_trial); % [NTg1bases_test x NTg1bases_trial]
    
    matElem_uv = matElem_uv ...
        + sparse(Row_elem(:), Col_elem(:), ...
        TrialTest_temp(:), ...
        Nelems*NTbases_test, Nelems*NTbases_trial);
    
    matElem_DuDv = matElem_DuDv ...
        + sparse(Row_elem(:), Col_elem(:), ...
        DxTrialDxTest_temp(:) + DyTrialDyTest_temp(:), ...
        Nelems*NTbases_test, Nelems*NTbases_trial);
    
    matElem_uvx = matElem_uvx ...
        + sparse(Row_elem(:), Col_elem(:), ...
        TrialDxTest_temp(:), ...
        Nelems*NTbases_test, Nelems*NTbases_trial);
    
    %-- value_rhs_f -- add to the position in the big vector
    Row_vec = ((ii_twoDomain-1)*NTbases_test+1 : ii_twoDomain*NTbases_test);
    vecRhsF(Row_vec,1) = vecRhsF(Row_vec,1) + rhs_f_temp;
    
    %<<-- End Part II --------------------------------------------------------------------------------------
    
end % for ii

end % function matElemCoeffsTrialTest



%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function coordTri0Elem = getcoordTri0Elem(singleNE, beginP_n, baryElem, coordv)
%
%
%   input:
%       beginP_n, the n-th point of singleElem, and from the n-th point to
%       construct the little triangles.
%

m = @(x) mod(x,singleNE)+(x==singleNE)*singleNE;
n = beginP_n;

if n == 1
    if singleNE == 3
        coordTri0Elem = coordv; % [3 x 2]
    elseif singleNE == 4
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 5
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(1,:)];
    elseif singleNE == 6
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 7
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 8
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(8,:);
            baryElem; coordv(8,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 9
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(8,:);
            baryElem; coordv(8,:); coordv(9,:);
            baryElem; coordv(9,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 10
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(8,:);
            baryElem; coordv(8,:); coordv(9,:);
            baryElem; coordv(9,:); coordv(10,:);
            baryElem; coordv(10,:); coordv(1,:)]; % [3*singleNE x 2]
    end % if singleNE == 3
    
else
    if singleNE == 3
        coordTri0Elem = coordv; % [3 x 2]
    elseif singleNE == 4
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:)];
    elseif singleNE == 5
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:)];
    elseif singleNE == 6
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:)];
    elseif singleNE == 7
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:)];  
    elseif singleNE == 8
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:); ...
            coordv(m(n),:); coordv(m(n+6),:); coordv(m(n+7),:)];  
    elseif singleNE == 9
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:); ...
            coordv(m(n),:); coordv(m(n+6),:); coordv(m(n+7),:); ...
            coordv(m(n),:); coordv(m(n+7),:); coordv(m(n+8),:)];
    elseif singleNE == 10
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:); ...
            coordv(m(n),:); coordv(m(n+6),:); coordv(m(n+7),:); ...
            coordv(m(n),:); coordv(m(n+7),:); coordv(m(n+8),:); ...
            coordv(m(n),:); coordv(m(n+8),:); coordv(m(n+9),:)];
    end % if singleNE == 3
    
end % if n==1

end % function getcoordTri0Elem
%%<<-- End sub function 1 ---------------------------------------------------------------

%%>> -- Begin sub function 2 -------------------------------------------------------------
function [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D)
%
%   output:
%       phyGpoints, [Npoints x 2]
%       phyGweights, [Npoints x 1]
x1=coordTri_nt(1,1);
y1=coordTri_nt(1,2);
x2=coordTri_nt(2,1);
y2=coordTri_nt(2,2);
x3=coordTri_nt(3,1);
y3=coordTri_nt(3,2);
JacobiTri=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));

phyGweights = JacobiTri * formulaGauss2D(:,3);
phyGpoints(:,1)=x1+(x2-x1)*formulaGauss2D(:,1)+(x3-x1)*formulaGauss2D(:,2);
phyGpoints(:,2)=y1+(y2-y1)*formulaGauss2D(:,1)+(y3-y1)*formulaGauss2D(:,2);
end % function getGaussLocalTri
%%<<-- End sub function 2 ---------------------------------------------------------------

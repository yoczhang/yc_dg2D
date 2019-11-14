function [g1_matElemCoeffsTrialTest, g1_matElemCoeffsDxTrialDxTest, g1_matElemCoeffsDyTrialDyTest, rhs_fh, varargout] = ...
    matElemCoeffsDtrialDtestPoisson_backup(Coeffs_func, rhs_f_func, meshInfo, formulaGauss2D, g1Trial_k, g1Test_k, varargin)
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
%       varargin, if length(varargin)~=0, the varargin must include g2Trial_k, g2Test_k, (i.e. 
%                     group2 trial polynomian degree k, and group2 test polynominal degree k).
%
%   output:
%       g1_matElemCoeffsTrialTest, [Nelems*NTg1bases_test x NTelems*NTg1bases_trial], \sum_e \int_e {grad u}\cdot[v]. 
%       g1_matElemCoeffsDxTrialDxTest, [Nelems*NTg1bases_test x NTelems*NTg1bases_trial], \sum_e \int_e {grad v}\cdot[u]. 
%       g1_matElemCoeffsDyTrialDyTest, [Nelems*NTg1bases_test x NTelems*NTg1bases_trial], \sum_e \int_e [u]\cdot[v]. 
%       rhs_fh, cell-type.
%       varargout, cell-type, matlab built-in var.
%
%
%   YcZhang 11/8/2017
%
%   Last modified 11/8/2017
%
%

% mesh information: interior edges 
Nelems = meshInfo.Nelems;

Coeff_k11 = Coeffs_func{1};
Coeff_k22 = Coeffs_func{2};
Coeff_1 = Coeffs_func{3};


% group1 setting
NTg1bases_trial = (g1Trial_k+1)*(g1Trial_k+2)/2;
NTg1bases_test = (g1Test_k+1)*(g1Test_k+2)/2;

g1_matElemCoeffsTrialTest = sparse(Nelems*NTg1bases_test,Nelems*NTg1bases_trial);
g1_matElemCoeffsDxTrialDxTest = sparse(Nelems*NTg1bases_test,Nelems*NTg1bases_trial);
%g1_matElemCoeffsDxTrialDyTest = saprse(Nelems*NTg1bases_test,Nelems*NTg1bases_trial);
%g1_matElemCoeffsDyTrialDxTest = saprse(Nelems*NTg1bases_test,Nelems*NTg1bases_trial);
g1_matElemCoeffsDyTrialDyTest = sparse(Nelems*NTg1bases_test,Nelems*NTg1bases_trial);

% group 2 setting
have_group2 = 0;
if length(varargin)==2
    have_group2 = 1;
    g2Trial_k = varargin{1};
    g2Test_k = varargin{2};   
    
    NTg2bases_trial = (g2Trial_k+1)*(g2Trial_k+2)/2;
    NTg2bases_test = (g2Test_k+1)*(g2Test_k+2)/2;
    
    g2_matElemCoeffsTrialTest = saprse(Nelems*NTg2bases_test,Nelems*NTg2bases_trial);
    g2_matElemCoeffsDxTrialDxTest = saprse(Nelems*NTg2bases_test,Nelems*NTg2bases_trial);
    g2_matElemCoeffsDyTrialDyTest = saprse(Nelems*NTg2bases_test,Nelems*NTg2bases_trial);
end

% rhs f function setting
rhs_fh = cell(length(rhs_f_func),1);
for ii = 1:length(rhs_f_func)
    rhs_fh{ii} = zeros(Nelems*NTg1bases_test,1);
end % for ii

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
    
    coordv = meshInfo.node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    baryElem = meshInfo.baryElem(ii,:); % [1 x 2], the barycenter (x-coord, y-coord) of element.
    
    coordTri0Elem = getcoordTri0Elem(singleNE, baryElem, coordv); 
        %> if singleNE==3, coordTri0Elem, [3 x 2]. else coordTri0Elem, [3*singleNE x 2].
    
    ele_xT = baryElem(1);
    ele_yT = baryElem(2); 
    ele_hT = meshInfo.diameters(ii); 
    
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %-- group 1
    g1_TrialTest_temp = zeros(NTg1bases_test, NTg1bases_trial);
    g1_DxTrialDxTest_temp = zeros(NTg1bases_test, NTg1bases_trial);
    %g1_DxTrialDyTest_temp = zeros(NTg1bases_test, NTg1bases_trial);
    %g1_DyTrialDxTest_temp = zeros(NTg1bases_test, NTg1bases_trial);
    g1_DyTrialDyTest_temp = zeros(NTg1bases_test, NTg1bases_trial);
    
    %-- group 2
    if have_group2
        g2_TrialTest_temp = zeros(NTg2bases_test, NTg2bases_trial);
        g2_DxTrialDxTest_temp = zeros(NTg2bases_test, NTg2bases_trial);
        g2_DyTrialDyTest_temp = zeros(NTg2bases_test, NTg2bases_trial);
    end % have_group2
    
    %-- rhs f function 
    value_rhs_f = cell(length(rhs_f_func),1);
    for jj = 1:length(rhs_f_func)
        value_rhs_f{jj} = zeros(NTg1bases_test,1);
    end % for jj
    
    
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D);

        % get the bases values on quad points
        [g1trialPb, g1trialPbx, g1trialPby] = localBases2D(ele_xT, ele_yT, ele_hT, phyGpoints(:,1), phyGpoints(:,2), g1Trial_k);
            %> g1trialPb, g1trialPbx, g1trialPby, [Npoints x NTg1bases_trial].
        [g1testPb, g1testPbx, g1testPby] = localBases2D(ele_xT, ele_yT, ele_hT, phyGpoints(:,1), phyGpoints(:,2), g1Test_k);
            %> g1testPb, g1testPbx, g1testPby, [Npoints x NTg1bases_test].
        
        %--- the funcValue may be chosen by case.
        valueCoeff_1 = Coeff_1(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k11 = Coeff_k11(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k22 = Coeff_k22(phyGpoints(:,1), phyGpoints(:,2));
        
        %-- group 1
        g1_TrialTest_temp = g1_TrialTest_temp + getMatOnElem(valueCoeff_1, g1trialPb, g1testPb, phyGweights);
        g1_DxTrialDxTest_temp = g1_DxTrialDxTest_temp + getMatOnElem(valueCoeff_k11, g1trialPbx, g1testPbx, phyGweights);
        g1_DyTrialDyTest_temp = g1_DyTrialDyTest_temp + getMatOnElem(valueCoeff_k22, g1trialPby, g1testPby, phyGweights);
        
        %-- group 2
        if have_group2
            [g2trialPb, g2trialPbx, g2trialPby] = localBases2D(ele_xT, ele_yT, ele_hT, phyGpoints(:,1), phyGpoints(:,2), g1Trial_k);
                %> g2testPb, g2testPbx, g2testPby, [Npoints x NTg2bases_trial].
            [g2testPb, g2testPbx, g2testPby] = localBases2D(ele_xT, ele_yT, ele_hT, phyGpoints(:,1), phyGpoints(:,2), g1Test_k);
                %> g2testPb, g2testPbx, g2testPby, [Npoints x NTg2bases_test].
                
            g2_TrialTest_temp = g2_TrialTest_temp + getMatOnElem(valueCoeff_1, g2trialPb, g2testPb, phyGweights);
            g2_DxTrialDxTest_temp = g2_DxTrialDxTest_temp + getMatOnElem(valueCoeff_k11, g2trialPbx, g2testPbx, phyGweights);
            g2_DyTrialDyTest_temp = g2_DyTrialDyTest_temp + getMatOnElem(valueCoeff_k22, g2trialPby, g2testPby, phyGweights);
        end % if have_group2   
        
        %-- rhs f integration        
        for jj = 1:length(rhs_f_func)
            value_rhs_f{jj} = value_rhs_f{jj} ...
                + g1testPb' * ( phyGweights.*rhs_f_func{jj}(phyGpoints(:,1), phyGpoints(:,2)) );
        end % for jj
        
    end % for nt
    
    %-- gruop 1 -- add to the position in the big mat
    g1_Row = ((ii-1)*NTg1bases_test+1 : ii*NTg1bases_test)' ...
        * ones(1,NTg1bases_trial); % [NTg1bases_test x NTg1bases_trial]
    g1_Col = ones(NTg1bases_test, 1) ...
        *((ii-1)*NTg1bases_trial+1 : ii*NTg1bases_trial); % [NTg1bases_test x NTg1bases_trial]
    
    g1_matElemCoeffsTrialTest = g1_matElemCoeffsTrialTest ...
        + sparse(g1_Row(:), g1_Col(:), g1_TrialTest_temp(:), Nelems*NTg1bases_test, Nelems*NTg1bases_trial);
    g1_matElemCoeffsDxTrialDxTest = g1_matElemCoeffsDxTrialDxTest ...
        + sparse(g1_Row(:), g1_Col(:), g1_DxTrialDxTest_temp(:), Nelems*NTg1bases_test, Nelems*NTg1bases_trial);
    g1_matElemCoeffsDyTrialDyTest = g1_matElemCoeffsDyTrialDyTest ...
        + sparse(g1_Row(:), g1_Col(:), g1_DyTrialDyTest_temp(:), Nelems*NTg1bases_test, Nelems*NTg1bases_trial);
    
    %-- gruop 2 -- add to the position in the big mat
    if have_group2
        g2_Row = ((ii-1)*NTg2bases_test+1 : ii*NTg2bases_test)' ...
            * ones(1,NTg2bases_trial); % [NTg2bases_test x NTg2bases_trial]
        g2_Col = ones(NTg2bases_test, 1) ...
            *((ii-1)*NTg2bases_trial+1 : ii*NTg2bases_trial); % [NTg2bases_test x NTg2bases_trial]

        g2_matElemCoeffsTrialTest = g2_matElemCoeffsTrialTest ...
            + sparse(g2_Row(:), g2_Col(:), g2_TrialTest_temp(:), Nelems*NTg2bases_test, Nelems*NTg2bases_trial);
        g2_matElemCoeffsDxTrialDxTest = g2_matElemCoeffsDxTrialDxTest ...
            + sparse(g2_Row(:), g2_Col(:), g2_DxTrialDxTest_temp(:), Nelems*NTg2bases_test, Nelems*NTg2bases_trial);
        g2_matElemCoeffsDyTrialDyTest = g2_matElemCoeffsDyTrialDyTest ...
            + sparse(g2_Row(:), g2_Col(:), g2_DyTrialDyTest_temp(:), Nelems*NTg2bases_test, Nelems*NTg2bases_trial);
    end % if have_group2
    
    %-- value_rhs_f -- add to the position in the big vector
    for jj = 1:length(rhs_f_func)
        rhs_fh{jj}((ii-1)*NTg1bases_test+1 : ii*NTg1bases_test) = rhs_fh{jj}((ii-1)*NTg1bases_test+1 : ii*NTg1bases_test) ...
            + value_rhs_f{jj}(:);
    end % for jj
    
    %<<-- End Part II --------------------------------------------------------------------------------------
    
end % for ii

if have_group2
    varargout{1} = g2_matElemCoeffsTrialTest;
    varargout{2} = g2_matElemCoeffsDxTrialDxTest;
    varargout{3} = g2_matElemCoeffsDyTrialDyTest;
end % if have_group2


end % function matElemCoeffsTrialTest



%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%

%%>>-- Begin sub function 1 ---------------------------------------------------------
function coordTri0Elem = getcoordTri0Elem(singleNE, baryElem, coordv)

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
    end 

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

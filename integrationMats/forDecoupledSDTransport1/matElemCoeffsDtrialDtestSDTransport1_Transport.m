function [C_G, C_GxxPlusGyy, C_u1G0xPlusu2G0y, vecRhsF] = ...
    matElemCoeffsDtrialDtestSDTransport1_Transport(pde, Coeffs_func, DarcymeshInfo, StokesmeshInfo, ...
    lastStep_DarcyUh, lastStep_StokesUh, formulaGauss2D, ...
    degreek_Darcy, degreek_Stokes, degreek_trialC, degreek_testC)
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       dgDecoupledSDTransport1, 
%       for the transport eqn, have NO inflow and outflow boundaryEdges,
%       i.e. for the transport eqn, all the boundaryEdges are Dirichlet Edges.
%   %---------------------------------------------   
%
%   In this function we compute the integral uv, u_xx\times v_xx + u_yy\times v_yy mat on element T.
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
%       matElem_uv, [Nelems*NTg1bases_test x NTelems*NTg1bases_trial], \sum_T \int_T  uv dx.
%       matElemDuDv, [Nelems*NTg1bases_test x NTelems*NTg1bases_trial], \sum_T \int_T u_xx v_xx + u_yy v_yy. 
%       rhs_fh, [Nelems*NTg1bases_test x 1].
%
%
%   YcZhang 11/8/2017
%
%   Last modified 16/8/2017
%
%

if isempty(lastStep_DarcyUh)
    error('in matElemCoeffsDtrialDtestSDTransport1_Transport.m, lastStep_DarcyUh is []')
end
if isempty(lastStep_StokesUh)
    error('in matElemCoeffsDtrialDtestSDTransport1_Transport.m, lastStep_StokesUh is []')
end

c_darcy = -pde.K;

S_Nelems = StokesmeshInfo.Nelems;
D_Nelems = DarcymeshInfo.Nelems;
S_Nbases = (degreek_Stokes+1)*(degreek_Stokes+2)/2;
D_Nbases = (degreek_Darcy+1)*(degreek_Darcy+2)/2;

dof_S_u1 = S_Nelems*S_Nbases;
% dof_S_p = S_Nelems*((degreek_Stokes-1+1)*(degreek_Stokes-1+2)/2);
% dof_D_u = D_Nelems*D_Nbases;
% dof_SD = (dof_D_u+2*dof_S_u1+dof_S_p);

Nelems = S_Nelems + D_Nelems;
% Nelems = meshInfo.Nelems;

Coeff_k11 = Coeffs_func{1};
Coeff_k22 = Coeffs_func{2};
Coeff_1 = @(x,y) 1+0.*x;
% vector_u_1 = Coeffs_func{3};
% vector_u_2 = Coeffs_func{4};

% group1 setting
C_Nbases_trial = (degreek_trialC+1)*(degreek_trialC+2)/2;
C_Nbases_test = (degreek_testC+1)*(degreek_testC+2)/2;

C_GxxPlusGyy = sparse(Nelems*C_Nbases_test,Nelems*C_Nbases_trial);
C_G = sparse(Nelems*C_Nbases_test,Nelems*C_Nbases_trial);
C_u1G0xPlusu2G0y = sparse(Nelems*C_Nbases_test,Nelems*C_Nbases_trial);
vecRhsF = zeros(Nelems*C_Nbases_test,1);


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

for CurrElem = 1:D_Nelems
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    singleElem = DarcymeshInfo.elem{CurrElem,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges on singleElem.
    
    concavePointElem_ii = DarcymeshInfo.concavePointElem(CurrElem); % [1 x 1], the concavePointElem of ii-the elem.
    coordv = DarcymeshInfo.node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    centroidElem = DarcymeshInfo.centroidElem(CurrElem,:); % [1 x 2], the centroid(xing xin) (x-coord, y-coord) of element.
    
    coordTri0Elem = getcoordTri0Elem(singleNE, concavePointElem_ii, centroidElem, coordv); 
        %> if singleNE==3, coordTri0Elem, [3 x 2]. else coordTri0Elem, [3*singleNE x 2].
    
    elem_xT = centroidElem(1);
    elem_yT = centroidElem(2); 
    elem_hT = DarcymeshInfo.hElem(CurrElem); 
    
    %-- get the Diffusivity coefficients of elemP and elemM
    elem_Darcy_DC = DarcymeshInfo.DiffusivityCoeffs(CurrElem);
    
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %-- the mat and vec setting
    C_G_temp = zeros(C_Nbases_test, C_Nbases_trial); % acutally, in the Poisson, this var will not be used.
    C_Gxx_temp = zeros(C_Nbases_test, C_Nbases_trial);
    C_Gyy_temp = zeros(C_Nbases_test, C_Nbases_trial);
    C_u1G0xPlusu2G0y_temp = zeros(C_Nbases_test, C_Nbases_trial);
    rhs_f_temp = zeros(C_Nbases_test,1);
    
    %-- the integration on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D);

        % get the bases values on quad points
        [trialPb_C, trialPbx_C, trialPby_C] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_trialC);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases_trial].
        [testPb_C, testPbx_C, testPby_C] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_testC);
            %> testPb, testPbx, testPby, [Npoints x NTbases_test].
        
        %--- the funcValue may be chosen by case.
        valueCoeff_1 = Coeff_1(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k11 = Coeff_k11(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k22 = Coeff_k22(phyGpoints(:,1), phyGpoints(:,2));
        valuefunc_f = pde.c0(0,phyGpoints(:,1), phyGpoints(:,2));
        elem_valueDarcy_DC = elem_Darcy_DC * ones(length(phyGweights),1);
        
        [~, trialPbx_D, trialPby_D] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_Darcy);
            %> [Npoints x Nbases]
            
        % in the Darcy domain
        elemUh = lastStep_DarcyUh((CurrElem-1)*D_Nbases+1:CurrElem*D_Nbases);
            
        value_vector_u_1 = c_darcy.*(trialPbx_D*elemUh);
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
        value_vector_u_2 = c_darcy.*(trialPby_D*elemUh);
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
                
        %------------------------ yc test ------------------------
% %         value_vector_u_1 = c_darcy.*(pde.phix(phyGpoints(:,1), phyGpoints(:,2)));
% %         value_vector_u_2 = c_darcy.*(pde.phiy(phyGpoints(:,1), phyGpoints(:,2)));
        %-----------------------------------------------------------
                
        C_u1G0xPlusu2G0y_temp = C_u1G0xPlusu2G0y_temp ...
            + getMatOnElem(value_vector_u_1, trialPb_C, testPbx_C, phyGweights) ...
            + getMatOnElem(value_vector_u_2, trialPb_C, testPby_C, phyGweights);
        
        
        
        %-- MatElem
        C_G_temp = C_G_temp + getMatOnElem(valueCoeff_1, trialPb_C, testPb_C, phyGweights);
        C_Gxx_temp = C_Gxx_temp + getMatOnElem(elem_valueDarcy_DC.*valueCoeff_k11, trialPbx_C, testPbx_C, phyGweights);
        C_Gyy_temp = C_Gyy_temp + getMatOnElem(elem_valueDarcy_DC.*valueCoeff_k22, trialPby_C, testPby_C, phyGweights);
        
        %-- VecF   
        rhs_f_temp = rhs_f_temp + testPb_C' * ( phyGweights.*valuefunc_f );
    end % for nt
    
    %-- add MatElem to the position in the big mat
    C_Row_elem = ((CurrElem-1)*C_Nbases_test+1 : (CurrElem)*C_Nbases_test)' ...
        * ones(1,C_Nbases_trial); % [NTg1bases_test x NTg1bases_trial]
    C_Col_elem = ones(C_Nbases_test, 1) ...
        *((CurrElem-1)*C_Nbases_trial+1 : (CurrElem)*C_Nbases_trial); % [NTg1bases_test x NTg1bases_trial]
    
    C_G = C_G ...
        + sparse(C_Row_elem(:), C_Col_elem(:), ...
        C_G_temp(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    
    C_GxxPlusGyy = C_GxxPlusGyy ...
        + sparse(C_Row_elem(:), C_Col_elem(:), ...
        C_Gxx_temp(:) + C_Gyy_temp(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    
    C_u1G0xPlusu2G0y = C_u1G0xPlusu2G0y ...
        + sparse(C_Row_elem(:), C_Col_elem(:), ...
        C_u1G0xPlusu2G0y_temp(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    
    %-- value_rhs_f -- add to the position in the big vector
    Row_vec = ((CurrElem-1)*C_Nbases_test+1 : (CurrElem)*C_Nbases_test);
    
    vecRhsF(Row_vec,1) = vecRhsF(Row_vec,1) + rhs_f_temp;
    
    %<<-- End Part II --------------------------------------------------------------------------------------
    
end % for ii
    
for CurrElem = 1:S_Nelems
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    singleElem = StokesmeshInfo.elem{CurrElem,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges on singleElem.
    
    concavePointElem_ii = StokesmeshInfo.concavePointElem(CurrElem); % [1 x 1], the concavePointElem of ii-the elem.
    coordv = StokesmeshInfo.node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    centroidElem = StokesmeshInfo.centroidElem(CurrElem,:); % [1 x 2], the centroid(xing xin) (x-coord, y-coord) of element.
    
    coordTri0Elem = getcoordTri0Elem(singleNE, concavePointElem_ii, centroidElem, coordv); 
        %> if singleNE==3, coordTri0Elem, [3 x 2]. else coordTri0Elem, [3*singleNE x 2].
    
    elem_xT = centroidElem(1);
    elem_yT = centroidElem(2); 
    elem_hT = StokesmeshInfo.hElem(CurrElem); 
    
    %-- get the Diffusivity coefficients of elemP and elemM
    elem_Stokes_DC = StokesmeshInfo.DiffusivityCoeffs(CurrElem);
    
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %-- the mat and vec setting
    C_G_temp = zeros(C_Nbases_test, C_Nbases_trial); % acutally, in the Poisson, this var will not be used.
    C_Gxx_temp = zeros(C_Nbases_test, C_Nbases_trial);
    C_Gyy_temp = zeros(C_Nbases_test, C_Nbases_trial);
    C_u1G0xPlusu2G0y_temp = zeros(C_Nbases_test, C_Nbases_trial);
    rhs_f_temp = zeros(C_Nbases_test,1);
    
    %-- the integration on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D);

        % get the bases values on quad points
        [trialPb_C, trialPbx_C, trialPby_C] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_trialC);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases_trial].
        [testPb_C, testPbx_C, testPby_C] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_testC);
            %> testPb, testPbx, testPby, [Npoints x NTbases_test].
        
        %--- the funcValue may be chosen by case.
        valueCoeff_1 = Coeff_1(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k11 = Coeff_k11(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k22 = Coeff_k22(phyGpoints(:,1), phyGpoints(:,2));
        valuefunc_f = pde.c0(0, phyGpoints(:,1), phyGpoints(:,2));
        elem_valueStokes_DC = elem_Stokes_DC * ones(length(phyGweights),1);
        
        % setting elem_flag and get the value of vector_u
        %elem_flag = settingElemFlagTransport(CurrElem); % here 
        %value_vector_u_1 = vector_u_1(phyGpoints(:,1), phyGpoints(:,2), elem_flag); % this value_vector_u is a vector-valued.
        %value_vector_u_2 = vector_u_2(phyGpoints(:,1), phyGpoints(:,2), elem_flag); % this value_vector_u is a vector-valued.
        
        
        [trialPb_Stokes, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_Stokes);
            %> [Npoints x Nbases]

        % in the Stokes domain
        elemUh1 = lastStep_StokesUh((CurrElem-1)*S_Nbases+1:CurrElem*S_Nbases);
        elemU2h = lastStep_StokesUh(...
            dof_S_u1+(CurrElem-1)*S_Nbases+1:dof_S_u1+CurrElem*S_Nbases);
        value_vector_u_1 = trialPb_Stokes*elemUh1; 
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
        value_vector_u_2 = trialPb_Stokes*elemU2h;
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
                
        %------------------------ yc test ------------------------
% %         value_vector_u_1 = (pde.u1(phyGpoints(:,1), phyGpoints(:,2)));
% %         value_vector_u_2 = (pde.u2(phyGpoints(:,1), phyGpoints(:,2)));
        %-----------------------------------------------------------
                           
        C_u1G0xPlusu2G0y_temp = C_u1G0xPlusu2G0y_temp ...
            + getMatOnElem(value_vector_u_1, trialPb_C, testPbx_C, phyGweights) ...
            + getMatOnElem(value_vector_u_2, trialPb_C, testPby_C, phyGweights);

        
        
        %-- MatElem
        C_G_temp = C_G_temp + getMatOnElem(valueCoeff_1, trialPb_C, testPb_C, phyGweights);
        C_Gxx_temp = C_Gxx_temp + getMatOnElem(elem_valueStokes_DC.*valueCoeff_k11, trialPbx_C, testPbx_C, phyGweights);
        C_Gyy_temp = C_Gyy_temp + getMatOnElem(elem_valueStokes_DC.*valueCoeff_k22, trialPby_C, testPby_C, phyGweights);

        %-- VecF          
        rhs_f_temp = rhs_f_temp + testPb_C' * ( phyGweights.*valuefunc_f );
    end % for nt
    
    %-- add MatElem to the position in the big mat
    C_Row_elem = ((D_Nelems+CurrElem-1)*C_Nbases_test+1 : (D_Nelems+CurrElem)*C_Nbases_test)' ...
        * ones(1,C_Nbases_trial); % [NTg1bases_test x NTg1bases_trial]
    C_Col_elem = ones(C_Nbases_test, 1) ...
        * ((D_Nelems+CurrElem-1)*C_Nbases_trial+1 : (D_Nelems+CurrElem)*C_Nbases_trial); % [NTg1bases_test x NTg1bases_trial]
    
    C_G = C_G ...
        + sparse(C_Row_elem(:), C_Col_elem(:), ...
        C_G_temp(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    
    C_GxxPlusGyy = C_GxxPlusGyy ...
        + sparse(C_Row_elem(:), C_Col_elem(:), ...
        C_Gxx_temp(:) + C_Gyy_temp(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    
    C_u1G0xPlusu2G0y = C_u1G0xPlusu2G0y ...
        + sparse(C_Row_elem(:), C_Col_elem(:), ...
        C_u1G0xPlusu2G0y_temp(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    
    %-- value_rhs_f -- add to the position in the big vector
    Row_vec = ((D_Nelems+CurrElem-1)*C_Nbases_test+1 : (D_Nelems+CurrElem)*C_Nbases_test);
    
    vecRhsF(Row_vec,1) = vecRhsF(Row_vec,1) + rhs_f_temp;
    
    %<<-- End Part II --------------------------------------------------------------------------------------
   
end % for CurrElem

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

if n == 0
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
    
end % if n==0

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

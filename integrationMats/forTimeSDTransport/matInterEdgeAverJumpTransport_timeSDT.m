function [C_diff_K_noAdapInterE, C_diff_K_AdapInterE, C_conv_K_InterE] = ...
    matInterEdgeAverJumpTransport_timeSDT(...
    pde, DarcymeshInfo, StokesmeshInfo, interfacemeshInfo, ...
    lastStep_DarcyUh, lastStep_StokesUh, option, formulaGauss1D, ...
    degreek_Darcy, degreek_Stokes, degreek_trialC, degreek_testC)
%
%   %-----------------------------------------
%       Just copy from matInterEdgeAverJumpSDTransport1_Transport.m
%   %-----------------------------------------
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
%   YcZhang 3/11/2017
%
%   Last modified 3/11/2017
%
%

if isempty(lastStep_DarcyUh)
    error('in matInterEdgeAverJumpSDTransport1_Transport.m, lastStep_DarcyUh is []')
end
if isempty(lastStep_StokesUh)
    error('in matInterEdgeAverJumpSDTransport1_Transport.m, lastStep_StokesUh is []')
end

c_darcy = -pde.K;

% mesh information: interior edges 
S_Nelems = StokesmeshInfo.Nelems;
D_Nelems = DarcymeshInfo.Nelems;
S_Nbases = (degreek_Stokes+1)*(degreek_Stokes+2)/2;
D_Nbases = (degreek_Darcy+1)*(degreek_Darcy+2)/2;

dof_S_u1 = S_Nelems*S_Nbases;
% dof_S_p = S_Nelems*((degreek_Stokes-1+1)*(degreek_Stokes-1+2)/2);
% dof_D_u = D_Nelems*D_Nbases;
% dof_SD = (dof_D_u+2*dof_S_u1+dof_S_p);

Nelems = S_Nelems + D_Nelems;
%Nelems = meshInfo.Nelems;
Darcy_interEdges = DarcymeshInfo.interEdgeIndex; % here we need the interEdges information.
Stokes_interEdges = StokesmeshInfo.interEdgeIndex;

p_epsilon_c = option.p_epsilon_c;
p_epsilon1 = option.p_epsilon1; % this is using to the adaptive-transport
p_sigma_c = option.p_sigma_c;
% p_sigma1 = option.p_sigma1; % this is using to the adaptive-transport
p_beta_c = option.p_beta_c;

% Coeff_k11 = Coeffs_func{1};
% Coeff_k22 = Coeffs_func{2};
% vector_u_1 = Coeffs_func{3};
% vector_u_2 = Coeffs_func{4};

% setting
C_Nbases_trial = (degreek_trialC+1)*(degreek_trialC+2)/2;
C_Nbases_test = (degreek_testC+1)*(degreek_testC+2)/2;

C_diff_K_noAdapInterE = sparse(Nelems*C_Nbases_test,Nelems*C_Nbases_trial);
C_diff_K_AdapInterE = sparse(Nelems*C_Nbases_test,Nelems*C_Nbases_trial);
C_conv_K_InterE = sparse(Nelems*C_Nbases_test,Nelems*C_Nbases_trial);

for CurrEdge = 1:length(Darcy_interEdges)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    interEindx = Darcy_interEdges(CurrEdge); 
    ePoint1 = DarcymeshInfo.node(DarcymeshInfo.edge(interEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = DarcymeshInfo.node(DarcymeshInfo.edge(interEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = DarcymeshInfo.areaEdge(interEindx);
    
    phyGpointsX = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = DarcymeshInfo.areaEdge(interEindx) * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
    elem1 = DarcymeshInfo.edge2elem(interEindx,1);
    elem2 = DarcymeshInfo.edge2elem(interEindx,2);
    local_e1 = DarcymeshInfo.edge2elem(interEindx,3);
    local_e2 = DarcymeshInfo.edge2elem(interEindx,4);
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
    
    elemP_nu = DarcymeshInfo.nuEdge0Elem{elemP}(:,eIndx_elemP); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleP.
    elemM_nu = DarcymeshInfo.nuEdge0Elem{elemM}(:,eIndx_elemM); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleM.
    
    elemP_xT = DarcymeshInfo.centroidElem(elemP,1);  
    elemP_yT = DarcymeshInfo.centroidElem(elemP,2);
    elemP_hT = DarcymeshInfo.hElem(elemP);
        %> the centroid(xing xin) and diameter of eleP.
    elemM_xT = DarcymeshInfo.centroidElem(elemM,1);  
    elemM_yT = DarcymeshInfo.centroidElem(elemM,2);
    elemM_hT = DarcymeshInfo.hElem(elemM);
        %> the centroid(xing xin) and diameter of eleM.
        
    %-- get the Diffusivity coefficients of elemP and elemM
    elemP_Darcy_DC = DarcymeshInfo.DiffusivityCoeffs(elemP);
    elemM_Darcy_DC = DarcymeshInfo.DiffusivityCoeffs(elemM);
    
    %-- subPart -- get the bases on quad
    [elemP_trialPb, elemP_trialPbx, elemP_trialPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_trialC);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
    [elemM_trialPb, elemM_trialPbx, elemM_trialPby] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_trialC);
        %> eleMtrialPb, eleMtrialPbx, eleMtrialPby, [Npoints x NTg1bases_trial]
    
    [elemP_testPb, elemP_testPbx, elemP_testPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_testC);
        %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]
    [elemM_testPb, elemM_testPbx, elemM_testPby] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_testC);
        %> eleMtestPb, eleMtestPbx, eleMtestPby, [Npoints x NTg1bases_test]
        
    %<<-- End Part I ---------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
    
    %% Part II, compute the 
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %
    %-- subPart I -- the AverJump of velotivy trial u and velotivy test v 
    % get the position in the big mat    
    [C_DarcyRow, C_DarcyCol] = getROWandCOL(elemP, elemM, C_Nbases_trial, C_Nbases_test);
    
    %------------------ convection term ----------------------------------
    [~, elemPtrialPbx_D, elemPtrialPby_D] = localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_Darcy);

	% in the Darcy domain
	elemPUh = lastStep_DarcyUh((elemP-1)*D_Nbases+1:elemP*D_Nbases);
	value_vector_u_1 = c_darcy.*(elemPtrialPbx_D*elemPUh);
        %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
	value_vector_u_2 = c_darcy.*(elemPtrialPby_D*elemPUh);
        %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
       
	%------------------------ yc test ------------------------
% % 	value_vector_u_1 = c_darcy.*(pde.phix(phyGpointsX, phyGpointsY));
% % 	value_vector_u_2 = c_darcy.*(pde.phiy(phyGpointsX, phyGpointsY));
	%-----------------------------------------------------------
            
	ConvectionUpwind_mat = getMatOnConvectionUpwindEdge(phyGweights, ...
        value_vector_u_1, ...
        value_vector_u_2, ...
        elemP_nu, ...
        elemP_trialPb, ...
        elemP_testPb, ...
        elemM_nu, ...
        elemM_trialPb, ...
        elemM_testPb);
        
	C_conv_K_InterE = C_conv_K_InterE ...
        + sparse(C_DarcyRow(:), C_DarcyCol(:), ...
        ConvectionUpwind_mat(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);

    %---------------------------------------------------------------------------------------------------

    %------------------------------------ diffusive term --------------------------------------------
    % get the matrix on all the interior edges, 
    elemP_valueDarcy_DC = elemP_Darcy_DC * ones(length(phyGweights),1);
	elemM_valueDarcy_DC = elemM_Darcy_DC * ones(length(phyGweights),1);
    [uv_A1, uv_A2, uv_A3] = ...
        getMatOnInterEdge(phyGweights, ...
        elemP_valueDarcy_DC, elemM_valueDarcy_DC, ...
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
            
    
    if elemP_Darcy_DC ~= elemM_Darcy_DC
        % the following H_ and P_ notations according to the paper:
        % discontinuous galerkin methods for convection-diffusion equations
        % for varying and vanishing diffusivity. ---- J.Proft and B.Riviere.
        % H_ denotes: Hyperbolic solution behavior (small diffusion conefficient),
        % P_ denotes: Parabolic solution behavior (larger diffusion conefficient).
        H_DC = min(elemP_Darcy_DC,elemM_Darcy_DC);
        P_DC = max(elemP_Darcy_DC,elemM_Darcy_DC);
        
        if H_DC == elemP_Darcy_DC
            H_elem = elemP;
            H_elem_nu = elemP_nu;
            
            %P_elem = elemM;
            %P_elem_nu = elemM_nu;
        else
            H_elem = elemM;
            H_elem_nu = elemM_nu;
            
            %P_elem = elemP;
            %P_elem_nu = elemP_nu;
        end 
    end 
    
    if elemP_Darcy_DC ~= elemM_Darcy_DC && (value_vector_u_1(1)*H_elem_nu(1)+value_vector_u_2(1)*H_elem_nu(2) > 0)
        theta_DC = H_DC/P_DC;
        %theta_DC=0;
        
        [Adap_A1_up, Adap_A2_up, Adap_A3_up] = ...
            getMatOnDiffusiveAdapEdgeHP_up(phyGweights, ...
            H_elem, elemP, ...
            elemP_valueDarcy_DC, elemM_valueDarcy_DC, ...
            elemP_nu, ...
            elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
            elemP_testPb, elemP_testPbx, elemP_testPby, ...
            elemM_nu, ...
            elemM_trialPb, elemM_trialPbx, elemM_trialPby, ...
            elemM_testPb, elemM_testPbx, elemM_testPby);
        
        C_diff_K_AdapInterE = C_diff_K_AdapInterE ...
            + sparse(C_DarcyRow(:), C_DarcyCol(:), ...
            -(1-theta_DC)*Adap_A1_up(:) - theta_DC*uv_A1(:) ...
            + p_epsilon1*( (1-theta_DC)*Adap_A2_up(:) + theta_DC*uv_A2(:) ) ...
            + (1-theta_DC)*(1/areaEdge^p_beta_c)*Adap_A3_up(:) + theta_DC*(p_sigma_c/areaEdge^p_beta_c)*uv_A3(:), ...
            Nelems*C_Nbases_test, Nelems*C_Nbases_trial);

    else
        C_diff_K_noAdapInterE = C_diff_K_noAdapInterE ...
            + sparse(C_DarcyRow(:), C_DarcyCol(:), ...
            - uv_A1(:) + p_epsilon_c*uv_A2(:) + (p_sigma_c/areaEdge^p_beta_c)*uv_A3(:), ...
            Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
        
    end % if
    %---------------------------------------------------------------------------------------------------

    %<<-- End Part II --------------------------------------------------------------------------------------
end % for CurrEdge


for CurrEdge = 1:interfacemeshInfo.Ninterface
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    ePoint1 = interfacemeshInfo.interfaceNode(interfacemeshInfo.interfaceEdge(CurrEdge,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = interfacemeshInfo.interfaceNode(interfacemeshInfo.interfaceEdge(CurrEdge,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = sqrt( sum((ePoint1-ePoint2).^2) );
    
    phyGpointsX = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = areaEdge * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
    elemS = interfacemeshInfo.interface2StokesElem(CurrEdge,1);
    eIndx_elemS = interfacemeshInfo.interface2StokesElem(CurrEdge,2);
        %> i.e., [n, e] = interface2StokesElem(ii,:), n standsfor the ii-the interface belongs to n-th elem, 
        %> and e standsfor the ii-th interface is the e-th edge of n-th elem. 
    elemD = interfacemeshInfo.interface2DarcyElem(CurrEdge,1);
    eIndx_elemD = interfacemeshInfo.interface2DarcyElem(CurrEdge,2);
    
    elemS_nu = StokesmeshInfo.nuEdge0Elem{elemS}(:,eIndx_elemS); 
    elemD_nu = DarcymeshInfo.nuEdge0Elem{elemD}(:,eIndx_elemD);
    
    elemS_xT = StokesmeshInfo.centroidElem(elemS,1);  
    elemS_yT = StokesmeshInfo.centroidElem(elemS,2);
    elemS_hT = StokesmeshInfo.hElem(elemS);
        %> the centroid(xing xin) and diameter of elemS.
        
    elemD_xT = DarcymeshInfo.centroidElem(elemD,1);  
    elemD_yT = DarcymeshInfo.centroidElem(elemD,2);
    elemD_hT = DarcymeshInfo.hElem(elemD);
        %> the centroid(xing xin) and diameter of elemD.
        
    %-- get the Diffusivity coefficients of elemP and elemM
    elemS_DC = StokesmeshInfo.DiffusivityCoeffs(elemS);
    elemD_DC = DarcymeshInfo.DiffusivityCoeffs(elemD);
        
    %-------------- elemP and elemM assignment --------------------- 
    % because the Darcy elems is the small elem, and Stokes is the bigger
    % elem, so:
    elemP = elemD; elemM = elemS; % here because elemM is the based on elemD.
    elemP_nu = elemD_nu; elemM_nu = elemS_nu;
    elemP_xT = elemD_xT; elemP_yT = elemD_yT; elemP_hT = elemD_hT; 
    elemM_xT = elemS_xT; elemM_yT = elemS_yT; elemM_hT = elemS_hT; 
    elemP_Darcy_DC = elemD_DC; elemM_Stokes_DC = elemS_DC;
    %--------------------------------------------------------------------------

    %-- subPart -- get the bases on quad
    [elemP_trialPb, elemP_trialPbx, elemP_trialPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_trialC);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
    [elemM_trialPb, elemM_trialPbx, elemM_trialPby] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_trialC);
        %> eleMtrialPb, eleMtrialPbx, eleMtrialPby, [Npoints x NTg1bases_trial]
    
    [elemP_testPb, elemP_testPbx, elemP_testPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_testC);
        %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]
    [elemM_testPb, elemM_testPbx, elemM_testPby] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_testC);
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
    [C_interfaceRow, C_interfaceCol] = getROWandCOL(elemP, D_Nelems+elemM, C_Nbases_trial, C_Nbases_test);
    
    %---------------------- convection term ------------------------------------------------------
    [elemMtrialPb_S, ~, ~] = localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_Stokes);

	% in the Stokes domain
	elemPUh1 = lastStep_StokesUh((elemS-1)*S_Nbases+1:elemS*S_Nbases);
	elemPUh2 = lastStep_StokesUh(dof_S_u1+(elemS-1)*S_Nbases+1:dof_S_u1+elemS*S_Nbases);
        
	% here we use the Stokes domain velocity as the input velocity.
	value_vector_u_1 = elemMtrialPb_S*elemPUh1;
        %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
	value_vector_u_2 = elemMtrialPb_S*elemPUh2;
        %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
            
	%------------------------ yc test ------------------------
% % 	value_vector_u_1 = (pde.u1(phyGpointsX, phyGpointsY));
% % 	value_vector_u_2 = (pde.u2(phyGpointsX, phyGpointsY));
	%-----------------------------------------------------------

	ConvectionUpwind_mat = getMatOnConvectionUpwindEdge(phyGweights, ...
        value_vector_u_1, ...
        value_vector_u_2, ...
        elemP_nu, ...
        elemP_trialPb, ...
        elemP_testPb, ...
        elemM_nu, ...
        elemM_trialPb, ...
        elemM_testPb);
        
	C_conv_K_InterE = C_conv_K_InterE ...
        + sparse(C_interfaceRow(:), C_interfaceCol(:), ...
        ConvectionUpwind_mat(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);

    %---------------------------------------------------------------------------------------------------

    %------------------------------------ diffusive term ------------------------------------------
    % get the matrix on all the interior edges, 
    elemP_valueDarcy_DC = elemP_Darcy_DC * ones(length(phyGweights),1);
	elemM_valueStokes_DC = elemM_Stokes_DC * ones(length(phyGweights),1);
    [uv_A1, uv_A2, uv_A3] = ...
        getMatOnInterEdge(phyGweights, ...
        elemP_valueDarcy_DC, elemM_valueStokes_DC, ...
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
            
    if elemP_Darcy_DC ~= elemM_Stokes_DC
        % the following H_ and P_ notations according to the paper:
        % discontinuous galerkin methods for convection-diffusion equations
        % for varying and vanishing diffusivity. ---- J.Proft and B.Riviere.
        % H_ denotes: Hyperbolic solution behavior (small diffusion conefficient),
        % P_ denotes: Parabolic solution behavior (larger diffusion conefficient).
        H_DC = min(elemP_Darcy_DC,elemM_Stokes_DC);
        P_DC = max(elemP_Darcy_DC,elemM_Stokes_DC);
        
        if H_DC == elemP_Darcy_DC
            H_elem = elemP;
            H_elem_nu = elemP_nu;
        else
            H_elem = D_Nelems+elemM;
            H_elem_nu = elemM_nu;
        end 
    end 
    
    if elemP_Darcy_DC ~= elemM_Stokes_DC && (value_vector_u_1(1)*H_elem_nu(1)+value_vector_u_2(1)*H_elem_nu(2) > 0)
        theta_DC = H_DC/P_DC;
        
        [Adap_A1_up, Adap_A2_up, Adap_A3_up] = ...
            getMatOnDiffusiveAdapEdgeHP_up(phyGweights, ...
            H_elem, elemP, ...
            elemP_valueDarcy_DC, elemM_valueStokes_DC, ...
            elemP_nu, ...
            elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
            elemP_testPb, elemP_testPbx, elemP_testPby, ...
            elemM_nu, ...
            elemM_trialPb, elemM_trialPbx, elemM_trialPby, ...
            elemM_testPb, elemM_testPbx, elemM_testPby);
        
        C_diff_K_AdapInterE = C_diff_K_AdapInterE ...
            + sparse(C_interfaceRow(:), C_interfaceCol(:), ...
            -(1-theta_DC)*Adap_A1_up(:) - theta_DC*uv_A1(:) ...
            + p_epsilon1*( (1-theta_DC)*Adap_A2_up(:) + theta_DC*uv_A2(:) ) ...
            + (1-theta_DC)*(1/areaEdge^p_beta_c)*Adap_A3_up(:) + theta_DC*(p_sigma_c/areaEdge^p_beta_c)*uv_A3(:), ...
            Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    else
        C_diff_K_noAdapInterE = C_diff_K_noAdapInterE ...
            + sparse(C_interfaceRow(:), C_interfaceCol(:), ...
            - uv_A1(:) + p_epsilon_c*uv_A2(:) + (p_sigma_c/areaEdge^p_beta_c)*uv_A3(:), ...
            Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    end % if
    %--------------------------------------------------------------------------------------------------
    

    %<<-- End Part II --------------------------------------------------------------------------------------
end % for CurrEdge


for CurrEdge = 1:length(Stokes_interEdges)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    interEindx = Stokes_interEdges(CurrEdge); 
    ePoint1 = StokesmeshInfo.node(StokesmeshInfo.edge(interEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = StokesmeshInfo.node(StokesmeshInfo.edge(interEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = StokesmeshInfo.areaEdge(interEindx);
    
    phyGpointsX = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = StokesmeshInfo.areaEdge(interEindx) * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
    elem1 = StokesmeshInfo.edge2elem(interEindx,1);
    elem2 = StokesmeshInfo.edge2elem(interEindx,2);
    local_e1 = StokesmeshInfo.edge2elem(interEindx,3);
    local_e2 = StokesmeshInfo.edge2elem(interEindx,4);
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
    
    elemP_nu = StokesmeshInfo.nuEdge0Elem{elemP}(:,eIndx_elemP); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleP.
    elemM_nu = StokesmeshInfo.nuEdge0Elem{elemM}(:,eIndx_elemM); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleM.
    
    elemP_xT = StokesmeshInfo.centroidElem(elemP,1);  
    elemP_yT = StokesmeshInfo.centroidElem(elemP,2);
    elemP_hT = StokesmeshInfo.hElem(elemP);
        %> the centroid(xing xin) and diameter of eleP.
    elemM_xT = StokesmeshInfo.centroidElem(elemM,1);  
    elemM_yT = StokesmeshInfo.centroidElem(elemM,2);
    elemM_hT = StokesmeshInfo.hElem(elemM);
        %> the centroid(xing xin) and diameter of eleM.
        
    %-- get the Diffusivity coefficients of elemP and elemM
    elemP_Stokes_DC = StokesmeshInfo.DiffusivityCoeffs(elemP);
    elemM_Stokes_DC = StokesmeshInfo.DiffusivityCoeffs(elemM);
        
%     %-- get the value of Coeffs func at Gauss Points
%     valueCoeff_k11 = Coeff_k11(phyGpointsX, phyGpointsY); % [Npoints x 1]
%     valueCoeff_k22 = Coeff_k22(phyGpointsX, phyGpointsY); % [Npoints x 1]
    
%     % setting elem_flag and get the value of vector_u
%     elem_flag = settingElemFlagAdapTransport(ii); % here 
%     value_vector_u_1 = vector_u_1(phyGpointsX, phyGpointsY, elem_flag); 
%         %> [Npoints x 1], this value_vector_u is a vector-valued.
%     value_vector_u_2 = vector_u_2(phyGpointsX, phyGpointsY, elem_flag); 
%         %> [Npoints x 1], this value_vector_u is a vector-valued.
    
    %-- subPart -- get the bases on quad
    [elemP_trialPb, elemP_trialPbx, elemP_trialPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_trialC);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
    [elemM_trialPb, elemM_trialPbx, elemM_trialPby] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_trialC);
        %> eleMtrialPb, eleMtrialPbx, eleMtrialPby, [Npoints x NTg1bases_trial]
    
    [elemP_testPb, elemP_testPbx, elemP_testPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_testC);
        %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]
    [elemM_testPb, elemM_testPbx, elemM_testPby] = ...
        localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_testC);
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
    [C_StokesRow, C_StokesCol] = getROWandCOL(D_Nelems+elemP, D_Nelems+elemM, C_Nbases_trial, C_Nbases_test);
    
    %---------------------- convection term ------------------------------------------------
    [elemPtrialPb_S, ~, ~] = localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_Stokes);

	% in the Stokes domain
	elemPUh1 = lastStep_StokesUh((elemP-1)*S_Nbases+1:elemP*S_Nbases);
	elemPUh2 = lastStep_StokesUh(dof_S_u1+(elemP-1)*S_Nbases+1:dof_S_u1+elemP*S_Nbases);
        
	value_vector_u_1 = elemPtrialPb_S*elemPUh1;
        %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
	value_vector_u_2 = elemPtrialPb_S*elemPUh2;
        %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
                
	%------------------------ yc test ------------------------
% % 	value_vector_u_1 = (pde.u1(phyGpointsX, phyGpointsY));
% % 	value_vector_u_2 = (pde.u2(phyGpointsX, phyGpointsY));
	%-----------------------------------------------------------
            
	ConvectionUpwind_mat = getMatOnConvectionUpwindEdge(phyGweights, ...
        value_vector_u_1, ...
        value_vector_u_2, ...
        elemP_nu, ...
        elemP_trialPb, ...
        elemP_testPb, ...
        elemM_nu, ...
        elemM_trialPb, ...
        elemM_testPb);
        
	C_conv_K_InterE = C_conv_K_InterE ...
        + sparse(C_StokesRow(:), C_StokesCol(:), ...
        ConvectionUpwind_mat(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);

    %---------------------------------------------------------------------------------------------

    %-------------------------- diffusive term ------------------------------------------------
    % get the matrix on all the interior edges, 
    elemP_valueStokes_DC = elemP_Stokes_DC * ones(length(phyGweights),1);
	elemM_valueStokes_DC = elemM_Stokes_DC * ones(length(phyGweights),1);
    [uv_A1, uv_A2, uv_A3] = ...
        getMatOnInterEdge(phyGweights, ...
        elemP_valueStokes_DC, elemM_valueStokes_DC, ...
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
            
    if elemP_Stokes_DC ~= elemM_Stokes_DC
        % the following H_ and P_ notations according to the paper:
        % discontinuous galerkin methods for convection-diffusion equations
        % for varying and vanishing diffusivity. ---- J.Proft and B.Riviere.
        % H_ denotes: Hyperbolic solution behavior (small diffusion conefficient),
        % P_ denotes: Parabolic solution behavior (larger diffusion conefficient).
        H_DC = min(elemP_Stokes_DC,elemM_Stokes_DC);
        P_DC = max(elemP_Stokes_DC,elemM_Stokes_DC);
        
        if H_DC == elemP_Stokes_DC
            H_elem = elemP;
            H_elem_nu = elemP_nu;
            
            %P_elem = elemM;
            %P_elem_nu = elemM_nu;
        else
            H_elem = elemM;
            H_elem_nu = elemM_nu;
            
            %P_elem = elemP;
            %P_elem_nu = elemP_nu;
        end 
    end 
    
    if elemP_Stokes_DC ~= elemM_Stokes_DC && (value_vector_u_1(1)*H_elem_nu(1)+value_vector_u_2(1)*H_elem_nu(2) > 0)
        theta_DC = H_DC/P_DC;
        %theta_DC=0;
        
        [Adap_A1_up, Adap_A2_up, Adap_A3_up] = ...
            getMatOnDiffusiveAdapEdgeHP_up(phyGweights, ...
            H_elem, elemP, ...
            elemP_valueStokes_DC, elemM_valueStokes_DC, ...
            elemP_nu, ...
            elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
            elemP_testPb, elemP_testPbx, elemP_testPby, ...
            elemM_nu, ...
            elemM_trialPb, elemM_trialPbx, elemM_trialPby, ...
            elemM_testPb, elemM_testPbx, elemM_testPby);
        
        C_diff_K_AdapInterE = C_diff_K_AdapInterE ...
            + sparse(C_StokesRow(:), C_StokesCol(:), ...
            -(1-theta_DC)*Adap_A1_up(:) - theta_DC*uv_A1(:) ...
            + p_epsilon1*( (1-theta_DC)*Adap_A2_up(:) + theta_DC*uv_A2(:) ) ...
            + (1-theta_DC)*(1/areaEdge^p_beta_c)*Adap_A3_up(:) + theta_DC*(p_sigma_c/areaEdge^p_beta_c)*uv_A3(:), ...
            Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    else
        C_diff_K_noAdapInterE = C_diff_K_noAdapInterE ...
            + sparse(C_StokesRow(:), C_StokesCol(:), ...
            - uv_A1(:) + p_epsilon_c*uv_A2(:) + (p_sigma_c/areaEdge^p_beta_c)*uv_A3(:), ...
            Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    end % if
    %---------------------------------------------------------------------------------------------
    


    %<<-- End Part II --------------------------------------------------------------------------------------
end % for CurrEdge

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
function [A1, A2, A3] = getMatOnDiffusiveAdapEdgeHP_up(phyGaussweights, ...
    H_elem, elemP, ...
    elemP_value_DC, elemM_value_DC, ...
    elemP_nu, ...
    elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
    elemP_testPb, elemP_testPbx, elemP_testPby, ...
    elemM_nu, ...
    elemM_trialPb, elemM_trialPbx, elemM_trialPby, ...
    elemM_testPb, elemM_testPbx, elemM_testPby)
%
%   %---------------------------------------------------------
%       get the mat on Diffusive adaptive edges, using the UPWIND scheme.
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


%-- first we need to allocate the matrix, because in the last some matrix
%-- may not be assigned values.
NTbases_trial = size(elemP_trialPbx,2);
NTbases_test = size(elemP_testPbx,2);
A1_pp = zeros(NTbases_test,NTbases_trial);
A1_pm = zeros(NTbases_test,NTbases_trial);
A1_mp = zeros(NTbases_test,NTbases_trial);
A1_mm = zeros(NTbases_test,NTbases_trial);

A2_pp = zeros(NTbases_test,NTbases_trial);
A2_pm = zeros(NTbases_test,NTbases_trial);
A2_mp = zeros(NTbases_test,NTbases_trial);
A2_mm = zeros(NTbases_test,NTbases_trial);

if H_elem == elemP
    A1_pp = elemP_nu(1) * elemP_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemP_trialPbx) ...
        + elemP_nu(2) * elemP_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemP_trialPby);
    A1_mp = elemM_nu(1) * elemM_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemP_trialPbx) ...
        + elemM_nu(2) * elemM_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemP_trialPby);
    
    A2_pp = elemP_nu(1) * (elemP_value_DC.*elemP_testPbx)' * bsxfun(@times,phyGaussweights,elemP_trialPb) ...
        + elemP_nu(2) * (elemP_value_DC.*elemP_testPby)' * bsxfun(@times,phyGaussweights,elemP_trialPb);
    A2_pm = elemM_nu(1) * (elemP_value_DC.*elemP_testPbx)' * bsxfun(@times,phyGaussweights,elemM_trialPb) ...
        + elemM_nu(2) * (elemP_value_DC.*elemP_testPby)' * bsxfun(@times,phyGaussweights,elemM_trialPb);
    
    %-- subPart III, \int_e [trial]\cdot[test]    
    A3_pp = elemP_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemP_trialPb);
    A3_pm = -elemP_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemM_trialPb);   
    A3_mp = -elemM_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemP_trialPb);
    A3_mm = elemM_testPb' * bsxfun(@times,elemP_value_DC.*phyGaussweights,elemM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].

else
    A1_pm = elemP_nu(1) * elemP_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemM_trialPbx) ...
        + elemP_nu(2) * elemP_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemM_trialPby);
	A1_mm = elemM_nu(1) * elemM_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemM_trialPbx) ...
        + elemM_nu(2) * elemM_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemM_trialPby);
    
    A2_mp = elemP_nu(1) * (elemM_value_DC.*elemM_testPbx)' * bsxfun(@times,phyGaussweights,elemP_trialPb) ...
        + elemP_nu(2) * (elemM_value_DC.*elemM_testPby)' * bsxfun(@times,phyGaussweights,elemP_trialPb);
	A2_mm = elemM_nu(1) * (elemM_value_DC.*elemM_testPbx)' * bsxfun(@times,phyGaussweights,elemM_trialPb) ...
        + elemM_nu(2) * (elemM_value_DC.*elemM_testPby)' * bsxfun(@times,phyGaussweights,elemM_trialPb);
    
    %-- subPart III, \int_e [trial]\cdot[test]    
    A3_pp = elemP_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemP_trialPb);
    A3_pm = -elemP_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemM_trialPb);   
    A3_mp = -elemM_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemP_trialPb);
    A3_mm = elemM_testPb' * bsxfun(@times,elemM_value_DC.*phyGaussweights,elemM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
end

% % %-- subPart III, \int_e [trial]\cdot[test]    
% % A3_pp = elemP_testPb' * bsxfun(@times,phyGaussweights,elemP_trialPb);
% % 	%> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
% %     
% % A3_pm = -elemP_testPb' * bsxfun(@times,phyGaussweights,elemM_trialPb);
% % 	%> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
% %         
% % A3_mp = -elemM_testPb' * bsxfun(@times,phyGaussweights,elemP_trialPb);
% % 	%> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
% %         
% % A3_mm = elemM_testPb' * bsxfun(@times,phyGaussweights,elemM_trialPb);
% % 	%> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].


A1 = [A1_pp, A1_pm;
    A1_mp, A1_mm];

A2 = [A2_pp, A2_pm;
    A2_mp, A2_mm];

A3 = [A3_pp, A3_pm;
    A3_mp, A3_mm];

end % function
        


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

for ii = 1:interfacemeshInfo.Ninterface
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    ePoint1 = interfacemeshInfo.interfaceNode(interfacemeshInfo.interfaceEdge(ii,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = interfacemeshInfo.interfaceNode(interfacemeshInfo.interfaceEdge(ii,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = sqrt( sum((ePoint1-ePoint2).^2) );
    
    phyGpointsX = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = areaEdge * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
    elemS = interfacemeshInfo.interface2StokesElem(ii,1);
    eIndx_elemS = interfacemeshInfo.interface2StokesElem(ii,2);
        %> i.e., [n, e] = interface2StokesElem(ii,:), n standsfor the ii-the interface belongs to n-th elem, 
        %> and e standsfor the ii-th interface is the e-th edge of n-th elem. 
    elemD = interfacemeshInfo.interface2DarcyElem(ii,1);
    eIndx_elemD = interfacemeshInfo.interface2DarcyElem(ii,2);
    
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
    elemP = elemD; elemM = D_Nelems+elemS; % here because elemM is the based on elemD.
    elemP_nu = elemD_nu; elemM_nu = elemS_nu;
    elemP_xT = elemD_xT; elemP_yT = elemD_yT; elemP_hT = elemD_hT; 
    elemM_xT = elemS_xT; elemM_yT = elemS_yT; elemM_hT = elemS_hT; 
    elemP_DC = elemD_DC; elemM_DC = elemS_DC;
    %--------------------------------------------------------------------------
    
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
    [C_interfaceRow, C_interfaceCol] = getROWandCOL(elemP, elemM, C_Nbases_trial, C_Nbases_test);
    
    %---------------------- convection term ------------------------------------------------------
    [elemPtrialPb_D, elemPtrialPbx_D, elemPtrialPby_D] = localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_Darcy);
    [elemMtrialPb_S, ~, ~] = localBases2D(elemM_xT, elemM_yT, elemM_hT, phyGpointsX, phyGpointsY, degreek_Stokes);
    if ~isempty(lastStep_StokesU)
        % in the Stokes domain
        elemPUh1 = lastStep_StokesU((elemS-1)*S_Nbases:elemS*S_Nbases);
        elemPUh2 = lastStep_StokesU(dof_S_u1+(elemS-1)*S_Nbases:dof_S_u1+elemS*S_Nbases);
        
        % here we use the Stokes domain velocity as the input velocity.
        value_vector_u_1 = elemMtrialPb_S*elemPUh1;
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
        value_vector_u_2 = elemMtrialPb_S*elemPUh2;
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
                
        ConvectionUpwind_mat = getMatOnConvectionUpwindEdge(phyGweights, ...
            value_vector_u_1, ...
            value_vector_u_2, ...
            elemP_nu, ...
            elemP_trialPb, ...
            elemP_testPb, ...
            elemM_nu, ...
            elemM_trialPb, ...
            elemM_testPb);
        
        IC_conv_K = IC_conv_K ...
            + sparse(C_interfaceRow(:), C_interfaceCol(:), ...
            ConvectionUpwind_mat(:), ...
            Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    else
        error('in matInterEdgeAverJumpSDTransport_Transport.m, lastStep_DarcyU is []')
    end
        
	if ~isempty(lastStep_C)
        elemP_Ch = lastStep_C((elemD-1)*C_Nbases_trial+1 : elemD*C_Nbases_trial);
        elemM_Ch = lastStep_C((elemS-1)*C_Nbases_trial+1 : elemS*C_Nbases_trial);
        value_elemPCh = elemP_trialPb*elemP_Ch;
        value_elemMCh = elemM_trialPb*elemM_Ch;
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
        value_elemCh = (value_elemPCh+value_elemMCh)/2;
        
        %---- Darcy domain --------------
        IDC_conv_K_elemPP = elemP_nu(1) * elemP_testPb' * bsxfun(@times,value_elemCh.*phyGweights,elemPtrialPbx_D) ...
            + elemP_nu(2) * elemP_testPb' * bsxfun(@times,value_elemCh.*phyGweights,elemPtrialPby_D);
        [IDC_Row, IDC_COL] = getROWandCOL(elemP, 0, D_Nbases, C_Nbases_test);
        IC_conv_K = IC_conv_K ...
            + sparse(IDC_Row(:), IDC_COL(:), ...
            IDC_conv_K_elemPP(:), ...
            Nelems*C_Nbases_test, dof_SD);
        %-------------------------------------
        
        %---- Stokes domain -------------
        ISC_conv_K_elemMM_u1 = elemM_nu(1) * elemM_testPb' * bsxfun(@times,value_elemCh.*phyGweights,elemMtrialPb_S);
        ISC_conv_K_elemMM_u2 = elemM_nu(2) * elemM_testPb' * bsxfun(@times,value_elemCh.*phyGweights,elemMtrialPb_S);
        [ISC_Row_u1, ISC_COL_u1] = getROWandCOL(elemM, 0, S_Nbases, C_Nbases_test);
        [ISC_Row_u2, ISC_COL_u2] = getROWandCOL(elemM, 0, S_Nbases, C_Nbases_test);
        ISC_COL_u2 = dof_S_u1 + ISC_COL_u2;
        
        IC_conv_K = IC_conv_K ...
            + sparse(ISC_Row_u1(:), ISC_COL_u1(:), ...
            ISC_conv_K_elemMM_u1(:), ...
            Nelems*C_Nbases_test, dof_SD);
        IC_conv_K = IC_conv_K ...
            + sparse(ISC_Row_u2(:), ISC_COL_u2(:), ...
            ISC_conv_K_elemMM_u2(:), ...
            Nelems*C_Nbases_test, dof_SD);
        %-------------------------------------
        
    else
        error('in matElemCoeffsDtrialDtestSDTransport_Transport.m, lastStep_C is []')
	end
    %---------------------------------------------------------------------------------------------------

    %------------------------------------ diffusive term ------------------------------------------
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
            
    if elemP_DC ~= elemM_DC
        % the following H_ and P_ notations according to the paper:
        % discontinuous galerkin methods for convection-diffusion equations
        % for varying and vanishing diffusivity. ---- J.Proft and B.Riviere.
        % H_ denotes: Hyperbolic solution behavior (small diffusion conefficient),
        % P_ denotes: Parabolic solution behavior (larger diffusion conefficient).
        H_DC = min(elemP_DC,elemM_DC);
        P_DC = max(elemP_DC,elemM_DC);
        
        if H_DC == elemP_DC
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
    
    if elemP_DC ~= elemM_DC && (value_vector_u_1(1)*H_elem_nu(1)+value_vector_u_2(1)*H_elem_nu(2) > 0)
        theta_DC = H_DC/P_DC;
        %theta_DC=0;
        
        [Adap_A1_up, Adap_A2_up, Adap_A3_up] = ...
            getMatOnDiffusiveAdapEdgeHP_up(phyGweights, ...
            H_elem, elemP, ...
            elemP_value_DC, elemM_value_DC, ...
            elemP_nu, ...
            elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
            elemP_testPb, elemP_testPbx, elemP_testPby, ...
            elemM_nu, ...
            elemM_trialPb, elemM_trialPbx, elemM_trialPby, ...
            elemM_testPb, elemM_testPbx, elemM_testPby);
        
        C_diff_K_AdapE = C_diff_K_AdapE ...
            + sparse(C_interfaceRow(:), C_interfaceCol(:), ...
            -(1-theta_DC)*Adap_A1_up(:) - theta_DC*uv_A1(:) ...
            + p_epsilon1*( (1-theta_DC)*Adap_A2_up(:) + theta_DC*uv_A2(:) ) ...
            + (1-theta_DC)*(1/areaEdge^p_beta)*Adap_A3_up(:) + theta_DC*(p_sigma/areaEdge^p_beta)*uv_A3(:), ...
            Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    else
        C_diff_K_noAdapE = C_diff_K_noAdapE ...
            + sparse(C_interfaceRow(:), C_interfaceCol(:), ...
            - uv_A1(:) + p_epsilon*uv_A2(:) + (p_sigma/areaEdge^p_beta)*uv_A3(:), ...
            Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    end % if
    %--------------------------------------------------------------------------------------------------
    

    %<<-- End Part II --------------------------------------------------------------------------------------
end % for ii
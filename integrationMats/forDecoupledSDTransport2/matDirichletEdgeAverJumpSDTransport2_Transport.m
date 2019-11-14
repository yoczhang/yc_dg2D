function [diff_DiriK, conv_DiriK] = matDirichletEdgeAverJumpSDTransport2_Transport(...
    pde, DarcymeshInfo, StokesmeshInfo, ...
    lastStep_DarcyUh, lastStep_StokesUh, option, formulaGauss1D, ...
    degreek_Darcy, degreek_Stokes, degreek_trialC, degreek_testC)
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       dgDecoupledSDTransport2, 
%       for the transport eqn, have NO inflow and outflow boundaryEdges.
%   %---------------------------------------------
%
%   %-----------------------------------------------------------------
%       Just copy from matDirichletEdgeAverJumpTransport.m
%   %-----------------------------------------------------------------
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
%   YcZhang 3/11/2017
%
%   Last modified 3/11/2017
%
%

c_darcy = -pde.K/pde.mu;
Coeff_k11 = @(x,y) 1 + 0.*x;
Coeff_k22 = @(x,y) 1 + 0.*x;

% mesh information: Dirichlet edges 
S_Nelems = StokesmeshInfo.Nelems;
D_Nelems = DarcymeshInfo.Nelems;
S_Nbases = (degreek_Stokes+1)*(degreek_Stokes+2)/2;
D_Nbases = (degreek_Darcy+1)*(degreek_Darcy+2)/2;

dof_S_u1 = S_Nelems*S_Nbases;

Nelems = D_Nelems + S_Nelems;

Darcy_DirichletEdges = DarcymeshInfo.DirichletEdgeIndex; % here we need the DirichletEdges information.
Stokes_DirichletEdges = StokesmeshInfo.DirichletEdgeIndex; % here we need the DirichletEdges information.

p_epsilon_c = option.p_epsilon_c;
p_sigma_c = option.p_sigma_c;
p_beta_c = option.p_beta_c;

%--- check the LastStep Uh
if isempty(lastStep_DarcyUh)
    error('in matDirichletEdgeAverJumpSDTransport2_Transport.m, lastStep_DarcyUh is []')
end
if isempty(lastStep_StokesUh)
    error('in matDirichletEdgeAverJumpSDTransport2_Transport.m, lastStep_StokesUh is []')
end

%--- setting Nbases
C_Nbases_trial = (degreek_trialC+1)*(degreek_trialC+2)/2;
C_Nbases_test = (degreek_testC+1)*(degreek_testC+2)/2;

diff_DiriK = sparse(Nelems*C_Nbases_test,Nelems*C_Nbases_trial);
conv_DiriK = sparse(Nelems*C_Nbases_test,Nelems*C_Nbases_trial);


for CurrEdge = 1:length(Darcy_DirichletEdges)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    DirichletEindx = Darcy_DirichletEdges(CurrEdge); 
    ePoint1 = DarcymeshInfo.node(DarcymeshInfo.edge(DirichletEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = DarcymeshInfo.node(DarcymeshInfo.edge(DirichletEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = DarcymeshInfo.areaEdge(DirichletEindx);
    
    phyGpointsX = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = areaEdge * formulaGauss1D(:,2);
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
    elem1 = DarcymeshInfo.edge2elem(DirichletEindx,1);
    elem2 = DarcymeshInfo.edge2elem(DirichletEindx,2);
    local_e1 = DarcymeshInfo.edge2elem(DirichletEindx,3);
    local_e2 = DarcymeshInfo.edge2elem(DirichletEindx,4);
        %> [n1, n2]=edge2elem(k,1:2) stands for the n1-th elem and n2-th elem share the k-th edge.
        %> local_e1=edge2elem(k,3) stands for the local edge index in n1-th elem of the k-th edge .
        %> also, local_e2=edge2elem(k,4) stands for the local edge index in n2-th elem of the k-th edge .
    elemP = max(elem1, elem2); % for the boundary edge, there will be one zero of (ele1,ele2), for let the bigger number of element as the element +.
    
    if elemP == elem1
        eIndx_elemP = local_e1;
    else
        eIndx_elemP = local_e2;
    end 
    
    elemP_nu = DarcymeshInfo.nuEdge0Elem{elemP}(:,eIndx_elemP); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleP.
    
    elemP_xT = DarcymeshInfo.centroidElem(elemP,1);  
    elemP_yT = DarcymeshInfo.centroidElem(elemP,2);
    elemP_hT = DarcymeshInfo.hElem(elemP);
        %> the centroid(xing xin) and diameter of eleP.
        
    %-- get the value of Coeffs func at Gauss Points
    valueCoeff_k11 = Coeff_k11(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valueCoeff_k22 = Coeff_k22(phyGpointsX, phyGpointsY); % [Npoints x 1]  
    
    %-- subPart -- get the bases on quad
    [elemP_trialPb, elemP_trialPbx, elemP_trialPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_trialC);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
    
    [elemP_testPb, elemP_testPbx, elemP_testPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_testC);
        %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]
        
    elemP_Darcy_DC = DarcymeshInfo.DiffusivityCoeffs(elemP);
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
    [Row, Col] = getROWandCOL(elemP, 0, C_Nbases_trial, C_Nbases_test);

    %------------------------------------ convection term -------------------------------------
    [~, elemPtrialPbx_D, elemPtrialPby_D] = localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_Darcy);

	% in the Darcy domain
	elemPUh = lastStep_DarcyUh((elemP-1)*D_Nbases+1:elemP*D_Nbases);
	value_vector_u_1 = c_darcy.*(elemPtrialPbx_D*elemPUh);
        %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
	value_vector_u_2 = c_darcy.*(elemPtrialPby_D*elemPUh);
        %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
            
% % 	%------------------------ yc test ------------------------
% % 	value_vector_u_1 = c_darcy.*(pde.phix(phyGpointsX, phyGpointsY));
% % 	value_vector_u_2 = c_darcy.*(pde.phiy(phyGpointsX, phyGpointsY));
% % 	%-----------------------------------------------------------
                
	conv_mat_temp = elemP_nu(1) * elemP_testPb' * bsxfun(@times,value_vector_u_1.*phyGweights,elemP_trialPb) ...
        + elemP_nu(2) * elemP_testPb' * bsxfun(@times,value_vector_u_2.*phyGweights,elemP_trialPb);
        
	conv_DiriK = conv_DiriK ...
        + sparse(Row(:), Col(:), ...
        conv_mat_temp(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);

    %--------------------------------------------------------------------------------------------------

    %------------------------------------ diffusive term ------------------------------------------
    % get the matrix
    elemP_valueDarcy_DC = elemP_Darcy_DC * ones(length(phyGweights),1);
    [diff_uv_A1, diff_uv_A2, diff_uv_A3] = getMatOnDirichletEdge(phyGweights, ...
        elemP_valueDarcy_DC.*valueCoeff_k11, elemP_valueDarcy_DC.*valueCoeff_k22, ...
        elemP_nu, ...
        elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
        elemP_testPb, elemP_testPbx, elemP_testPby);
        %> uv_A1, uv_A2, uv_A3, [NTbases_test x NTbases_trial]
        %> uv_A1=[A1_pp];
        %> the same to uv_A2, uv_A3.
    
    diff_DiriK = diff_DiriK ...
        + sparse(Row(:), Col(:), ...
        - diff_uv_A1(:) + p_epsilon_c*diff_uv_A2(:) + (p_sigma_c/areaEdge^p_beta_c)*diff_uv_A3(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    %--------------------------------------------------------------------------------------------------

    %<<-- End Part II --------------------------------------------------------------------------------------
end % for CurrEdge

for CurrEdge = 1:length(Stokes_DirichletEdges)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    DirichletEindx = Stokes_DirichletEdges(CurrEdge); 
    ePoint1 = StokesmeshInfo.node(StokesmeshInfo.edge(DirichletEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = StokesmeshInfo.node(StokesmeshInfo.edge(DirichletEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = StokesmeshInfo.areaEdge(DirichletEindx);
    
    phyGpointsX = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = areaEdge * formulaGauss1D(:,2);
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
    elem1 = StokesmeshInfo.edge2elem(DirichletEindx,1);
    elem2 = StokesmeshInfo.edge2elem(DirichletEindx,2);
    local_e1 = StokesmeshInfo.edge2elem(DirichletEindx,3);
    local_e2 = StokesmeshInfo.edge2elem(DirichletEindx,4);
        %> [n1, n2]=edge2elem(k,1:2) stands for the n1-th elem and n2-th elem share the k-th edge.
        %> local_e1=edge2elem(k,3) stands for the local edge index in n1-th elem of the k-th edge .
        %> also, local_e2=edge2elem(k,4) stands for the local edge index in n2-th elem of the k-th edge .
    elemP = max(elem1, elem2); % for the boundary edge, there will be one zero of (ele1,ele2), for let the bigger number of element as the element +.
    
    if elemP == elem1
        eIndx_elemP = local_e1;
    else
        eIndx_elemP = local_e2;
    end 
    
    elemP_nu = StokesmeshInfo.nuEdge0Elem{elemP}(:,eIndx_elemP); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleP.
    
    elemP_xT = StokesmeshInfo.centroidElem(elemP,1);  
    elemP_yT = StokesmeshInfo.centroidElem(elemP,2);
    elemP_hT = StokesmeshInfo.hElem(elemP);
        %> the centroid(xing xin) and diameter of eleP.
        
    %-- get the value of Coeffs func at Gauss Points
    valueCoeff_k11 = Coeff_k11(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valueCoeff_k22 = Coeff_k22(phyGpointsX, phyGpointsY); % [Npoints x 1]  
    
    %-- subPart -- get the bases on quad
    [elemP_trialPb, elemP_trialPbx, elemP_trialPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_trialC);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
    
    [elemP_testPb, elemP_testPbx, elemP_testPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_testC);
        %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]
        
	elemP_Stokes_DC = StokesmeshInfo.DiffusivityCoeffs(elemP);
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
    [Row, Col] = getROWandCOL(D_Nelems+elemP, 0, C_Nbases_trial, C_Nbases_test);

    %------------------------------------ convection term -------------------------------------
    [elemPtrialPb_S, ~, ~] = localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_Darcy);

	% in the Stokes domain
	elemPUh1 = lastStep_StokesUh((elemP-1)*S_Nbases+1:elemP*S_Nbases);
	elemPUh2 = lastStep_StokesUh(dof_S_u1+(elemP-1)*S_Nbases+1:dof_S_u1+elemP*S_Nbases);
        
	value_vector_u_1 = elemPtrialPb_S*elemPUh1;
        %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
	value_vector_u_2 = elemPtrialPb_S*elemPUh2;
        %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]

% % 	%------------------------ yc test ------------------------
% % 	value_vector_u_1 = (pde.u1(phyGpointsX, phyGpointsY));
% % 	value_vector_u_2 = (pde.u2(phyGpointsX, phyGpointsY));
% % 	%-----------------------------------------------------------
                
	conv_mat_temp = elemP_nu(1) * elemP_testPb' * bsxfun(@times,value_vector_u_1.*phyGweights,elemP_trialPb) ...
        + elemP_nu(2) * elemP_testPb' * bsxfun(@times,value_vector_u_2.*phyGweights,elemP_trialPb);
        
	conv_DiriK = conv_DiriK ...
        + sparse(Row(:), Col(:), ...
        conv_mat_temp(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);

    %--------------------------------------------------------------------------------------------------

    %------------------------------------ diffusive term ------------------------------------------
    % get the matrix
    elemP_valueStokes_DC = elemP_Stokes_DC * ones(length(phyGweights),1);
    [diff_uv_A1, diff_uv_A2, diff_uv_A3] = getMatOnDirichletEdge(phyGweights, ...
        elemP_valueStokes_DC.*valueCoeff_k11, elemP_valueStokes_DC.*valueCoeff_k22, ...
        elemP_nu, ...
        elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
        elemP_testPb, elemP_testPbx, elemP_testPby);
        %> uv_A1, uv_A2, uv_A3, [NTbases_test x NTbases_trial]
        %> uv_A1=[A1_pp];
        %> the same to uv_A2, uv_A3.
    
    diff_DiriK = diff_DiriK ...
        + sparse(Row(:), Col(:), ...
        - diff_uv_A1(:) + p_epsilon_c*diff_uv_A2(:) + (p_sigma_c/areaEdge^p_beta_c)*diff_uv_A3(:), ...
        Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
    %--------------------------------------------------------------------------------------------------

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
function [int_Agradtrial_Jtest, int_Agradtest_Jtrial, int_Jtrial_Jtest] = getMatOnDirichletEdge(phyGaussweights, ...
    valueCoeff_1, valueCoeff_2, ...
    elemP_nu, ...
    elemP_trialPb, elemP_trialPbx, elemP_trialPby, ...
    elemP_testPb, elemP_testPbx, elemP_testPby)
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
    A1_pp = elemP_nu(1) * elemP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,elemP_trialPbx) ...
        + elemP_nu(2) * elemP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,elemP_trialPby);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    int_Agradtrial_Jtest = A1_pp;
        
    
    %-- subPart II, \int_e {(func1,func2)\cdot grad_test}\cdot[trial]
    A2_pp = elemP_nu(1) * elemP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,elemP_trialPb) ...
        + elemP_nu(2) * elemP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,elemP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    int_Agradtest_Jtrial = A2_pp;    
    
        
    %-- subPart III, \int_e [trial]\cdot[test]    
    A3_pp = elemP_testPb' * bsxfun(@times,phyGaussweights,elemP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].   
    int_Jtrial_Jtest = A3_pp;
end % function getMatOnBoundEdge

%--------------------------------------------------------------------------------------------------------------------%
%----------- end sub function ----------------------------------------------------------------------------------%

function [interface_uphi_n1, interface_uphi_n2, interface_phiu_n1, interface_phiu_n2] = ...
    matSDTransport1_InterfaceStokesDarcy(Coeffs_func, StokesmeshInfo, DarcymeshInfo, interfacemeshInfo, formulaGauss1D, degreek_u, degreek_phi)
%
%   %---------------------------------------------------------------------
%       Just copy from: matInterfaceStokesDarcy.m
%   %----------------------------------------------------------------------
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       dgDecoupledSDTransport1, 
%       for the transport eqn, have NO inflow and outflow boundaryEdges,
%       i.e. for the transport eqn, all the boundaryEdges are Dirichlet Edges.
%   %---------------------------------------------
%
%   This function is to compute the coupled interface integration of Stoke and Darcy.
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases_trial denote the number of LOCAL trial bases on each K of Th.
%               NTbases_test denote the number of LOCAL test bases on each K of Th.
%   
%
%   input:
%       Coeffs, the cell-type, here 
%                   Coeffs{1}, and Coeffs{2} are the K=[k11, k12; k21, k22], 
%                   Coeffs{3} is the (\int_T Coeffs{1}uv dx),
%                   let k12=0, k21=0, the Coeffs{1}=k11, Coeffs{2}= k22. 
%                   Coeffs{1} is the (\int_T Coeffs{1} u_x v dx). 
%                   Coeffs{2} is the (\int_T Coeffs{2} u_y v dx).
%       StokesmeshInfo, mesh structure-type, the Stokes mesh information.
%       DarcymeshInfo, mesh structure-type, the Darcy mesh information.
%       interfacemeshInfo, mesh structure-type, the interface mesh information.
%       formulaGauss1D, the 1d Gauss quadrature formula, size: a matrix, [Npoints x 2],
%               the first column is the all Gauss-Points on [0,1],
%               the second is the weights of all Gauss-Points.
%       degreek_u, the Stokes u polynomial degree.
%       degreek_phi, the Darcy phi polynomial degree.
%
%   output:
%       interface_uphi, phi is the Darcy bases, is the test function, [Nelems_D*NTbases_phi x Nelems_S*NTbases_u].
%       interface_uphi, u is the Stokes bases, is the test function, [Nelems_S*NTbases_u x Nelems_D*NTbases_phi].
%
%
%   YcZhang 30/8/2017
%
%   Last modified 30/8/2017
%

% mesh information 
Nelems_S = StokesmeshInfo.Nelems;
Nelems_D = DarcymeshInfo.Nelems;

% bases setting
NTbases_u = (degreek_u+1)*(degreek_u+2)/2;
NTbases_phi = (degreek_phi+1)*(degreek_phi+2)/2;

% coefficient setting
Coeff_u = Coeffs_func{1};
Coeff_phi = Coeffs_func{2};

% output matrix setting
interface_uphi_n1 = sparse(Nelems_D*NTbases_phi,Nelems_S*NTbases_u);
interface_uphi_n2 = sparse(Nelems_D*NTbases_phi,Nelems_S*NTbases_u);
interface_phiu_n1 = sparse(Nelems_S*NTbases_u,Nelems_D*NTbases_phi);
interface_phiu_n2 = sparse(Nelems_S*NTbases_u,Nelems_D*NTbases_phi);

for ii = 1:interfacemeshInfo.Ninterface
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    ePoint1 = interfacemeshInfo.interfaceNode(interfacemeshInfo.interfaceEdge(ii,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = interfacemeshInfo.interfaceNode(interfacemeshInfo.interfaceEdge(ii,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = sqrt( sum((ePoint1-ePoint2).^2) );
    
    phyGpointsX = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = areaEdge * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
        
    elemS = interfacemeshInfo.interface2StokesElem(ii,1);
    eIndx_elemS = interfacemeshInfo.interface2StokesElem(ii,2);
        %> i.e., [n, e] = interface2StokesElem(ii,:), n standsfor the ii-the interface belongs to n-th elem, 
        %> and e standsfor the ii-th interface is the e-th edge of n-th elem. 
    elemD = interfacemeshInfo.interface2DarcyElem(ii,1);
    
    elemS_nu = StokesmeshInfo.nuEdge0Elem{elemS}(:,eIndx_elemS); 
    
    elemS_xT = StokesmeshInfo.centroidElem(elemS,1);  
    elemS_yT = StokesmeshInfo.centroidElem(elemS,2);
    elemS_hT = StokesmeshInfo.hElem(elemS);
        %> the centroid(xing xin) and diameter of elemS.
        
    elemD_xT = DarcymeshInfo.centroidElem(elemD,1);  
    elemD_yT = DarcymeshInfo.centroidElem(elemD,2);
    elemD_hT = DarcymeshInfo.hElem(elemD);
        %> the centroid(xing xin) and diameter of elemD.
    
    %-- get the value of Coeffs func at Gauss Points
    valueCoeff_u = Coeff_u(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valueCoeff_phi = Coeff_phi(phyGpointsX, phyGpointsY); % [Npoints x 1]
    
    %-- subPart -- get the bases on quad
    [u_elemS_Pb, ~, ~] = ...
        localBases2D(elemS_xT, elemS_yT, elemS_hT, phyGpointsX, phyGpointsY, degreek_u);
        %> u_elemS_Pb, [Npoints x NTbases_u]
    [phi_elemD_Pb, ~, ~] = ...
        localBases2D(elemD_xT, elemD_yT, elemD_hT, phyGpointsX, phyGpointsY, degreek_phi);
        %> phi_elemD_Pb, [Npoints x NTbases_phi]

    %<<-- End Part I ---------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
    
    %% Part II, 
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %
    %-- subPart 1 -- get the position in the big mat
    uphi_Row = ((elemD-1)*NTbases_phi+1 : elemD*NTbases_phi)' ... 
        * ones(1,NTbases_u); % [NTbases_phi x NTbases_u]
    uphi_Col = ones(NTbases_phi, 1) ...
        *((elemS-1)*NTbases_u+1 : elemS*NTbases_u); % [NTbases_phi x NTbases_u]
    
    phiu_Row = ((elemS-1)*NTbases_u+1 : elemS*NTbases_u)' ... 
        * ones(1,NTbases_phi); % [NTbases_u x NTbases_phi]
    phiu_Col = ones(NTbases_u, 1) ...
        *((elemD-1)*NTbases_phi+1 : elemD*NTbases_phi); % [NTbases_u x NTbases_phi]
    
    %-- subPart 2 -- get the integration mat
    uphin1_temp = elemS_nu(1) * phi_elemD_Pb' * bsxfun(@times,valueCoeff_u.*phyGweights,u_elemS_Pb);
    uphin2_temp = elemS_nu(2) * phi_elemD_Pb' * bsxfun(@times,valueCoeff_u.*phyGweights,u_elemS_Pb);
    %> [NTbases_phi x Npoints] * [Npoints x NTbases_u] = [NTbases_phi x NTbases_u].
    
    phiun1_temp = elemS_nu(1) * u_elemS_Pb' * bsxfun(@times,valueCoeff_phi.*phyGweights,phi_elemD_Pb);
    phiun2_temp = elemS_nu(2) * u_elemS_Pb' * bsxfun(@times,valueCoeff_phi.*phyGweights,phi_elemD_Pb);
    %> [NTbases_u x Npoints] * [Npoints x NTbases_phi] = [NTbases_u x NTbases_phi].
    
    %-- subPart 3 -- add to the big mat
    interface_uphi_n1 = interface_uphi_n1 ...
        + sparse(uphi_Row(:), uphi_Col(:), -uphin1_temp(:), Nelems_D*NTbases_phi, Nelems_S*NTbases_u);
    interface_uphi_n2 = interface_uphi_n2 ...
        + sparse(uphi_Row(:), uphi_Col(:), -uphin2_temp(:), Nelems_D*NTbases_phi, Nelems_S*NTbases_u);
    interface_phiu_n1 = interface_phiu_n1 ...
        + sparse(phiu_Row(:), phiu_Col(:), phiun1_temp(:), Nelems_S*NTbases_u, Nelems_D*NTbases_phi);
    interface_phiu_n2 = interface_phiu_n2 ...
        + sparse(phiu_Row(:), phiu_Col(:), phiun2_temp(:), Nelems_S*NTbases_u, Nelems_D*NTbases_phi);
    
    
end % for ii


end % function
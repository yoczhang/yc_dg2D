function [interface_u1v1_test, interface_u2v1_test, interface_u1v2_test, interface_u2v2_test] = ...
    matInterfaceOfStokes_timeSDT(Coeffs_func, StokesmeshInfo, interfacemeshInfo, formulaGauss1D, degreek_u)
%
%   %-----------------------------------------
%       Just copy from matInterfaceOfStokes_timeSD.m
%   %-----------------------------------------
%
%   %>-----------------------------------------------------------------------------------------------------------------------------------------
%       Here we proviede two methods to get the interface_u1v1, interface_u2v1, ....
%       the 1-th method, we just use the StokesmeshInfo, and it's StokesmeshInfo.interfaceEdgeIndex.
%       the 2-th method, we use both StokesmeshInfo and interfacemeshInfo.
%       In fact, the 1-th and 2-th method give the same results, See the test: aaa = interface_u1v1_test - interface_u1v1 in end.
%       And in the following we comment the 1-th method.
%   %>------------------------------------------------------------------------------------------------------------------------------------------
%
%   This function is to compute the interface integration of Stokes domain.
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
%       interfacemeshInfo, mesh structure-type, the interface mesh information.
%       formulaGauss1D, the 1d Gauss quadrature formula, size: a matrix, [Npoints x 2],
%               the first column is the all Gauss-Points on [0,1],
%               the second is the weights of all Gauss-Points.
%       degreek_u, the trial and test polynomial degree k.
%
%   output:
%       interface_u1v1, interface_u2v1, interface_u1v2, interface_u2v2,
%           [Nelems*NTbases_test x NTelems*NTbases_trial].
%
%
%
%   YcZhang 30/8/2017
%
%   Last modified 3/9/2017
%

% mesh information: interior edges 
Nelems = StokesmeshInfo.Nelems;

% bases setting
NTbases_u = (degreek_u+1)*(degreek_u+2)/2;

% coefficient setting
Coeff1 = Coeffs_func{1};
Coeff2 = Coeffs_func{2};

% % % output matrix setting
% % interface_u1v1 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
% % interface_u2v1 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
% % interface_u1v2 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
% % interface_u2v2 = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
% % 
% % % interface info setting
% % Ninterface = length(StokesmeshInfo.interfaceEdgeIndex);
% % interfaceEdges = StokesmeshInfo.interfaceEdgeIndex;
% 
% % for ii = 1:Ninterface
% %     %% Part I, get the information about 
% %     % 1. physical GaussPoints, 
% %     % 2. different element bases on phy GaussPoints on ii-th edge.
% %     %
% %     %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
% %     interfaceIndex = interfaceEdges(ii);
% %     ePoint1 = StokesmeshInfo.node(StokesmeshInfo.edge(interfaceIndex,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
% %     ePoint2 = StokesmeshInfo.node(StokesmeshInfo.edge(interfaceIndex,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
% %     
% %     phyGpointsX = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
% %         %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
% %     phyGpointsY = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
% %         %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
% %     phyGweights = StokesmeshInfo.areaEdge(interfaceIndex) * formulaGauss1D(:,2); 
% %         %> [Npoints x 1], the physical Gauss weights on ii-th edge.
% %         
% %     ele1 = StokesmeshInfo.edge2elem(interfaceIndex,1);
% %     ele2 = StokesmeshInfo.edge2elem(interfaceIndex,2);
% %     local_e1 = StokesmeshInfo.edge2elem(interfaceIndex,3);
% %     local_e2 = StokesmeshInfo.edge2elem(interfaceIndex,4);
% %         %> [n1, n2]=edge2elem(k,1:2) stands for the n1-th elem and n2-th elem share the k-th edge.
% %         %> local_e1=edge2elem(k,3) stands for the local edge index in n1-th elem of the k-th edge .
% %         %> also, local_e2=edge2elem(k,4) stands for the local edge index in n2-th elem of the k-th edge .
% %     eleP = max(ele1, ele2); % for the boundary edge, there will be one zero of (ele1,ele2), for let the bigger number of element as the element +.
% %     
% %     if eleP == ele1
% %         eIndx_eleP = local_e1;
% %     else
% %         eIndx_eleP = local_e2;
% %     end 
% %     
% %     eleP_nu = StokesmeshInfo.nuEdge0Elem{eleP}(:,eIndx_eleP); 
% %         %> [1 x 2], the outward unit normal vector of ii-th edge on eleP.
% %     eleP_tau = [-eleP_nu(2); eleP_nu(1)];
% %     
% %     eleP_xT = StokesmeshInfo.centroidElem(eleP,1);  
% %     eleP_yT = StokesmeshInfo.centroidElem(eleP,2);
% %     eleP_hT = StokesmeshInfo.hElem(eleP);
% %         %> the centroid(xing xin) and diameter of eleP.
% %         
% %     %-- get the value of Coeffs func at Gauss Points
% %     valueCoeff1 = Coeff1(phyGpointsX, phyGpointsY); % [Npoints x 1]
% %     valueCoeff2 = Coeff2(phyGpointsX, phyGpointsY); % [Npoints x 1]
% %     
% %     %-- subPart -- get the bases on quad
% %     [u_eleP_Pb, ~, ~] = ...
% %         localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, degreek_u);
% %         %> u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, [Npoints x NTbases_u]
% % 
% %     %<<-- End Part I ---------------------------------------------------------------------------------
% %     
% % %% ------------------------------------------ Dividing line -------------------------------------------------- %%
% % %-------------------------------------------------------------------------------------------------------------------%
% %     
% %     %% Part II, 
% %     %
% %     %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
% %     %
% %     %-- subPart 1 -- the AverJump of velotivy trial u and velotivy test v.
% %     % get the position in the big mat
% %     [uv_Row, uv_Col] = getROWandCOL(eleP, 0, NTbases_u, NTbases_u);
% %     
% %     u1v1_temp = eleP_tau(1) * eleP_tau(1) * u_eleP_Pb' * bsxfun(@times,valueCoeff1.*phyGweights,u_eleP_Pb);
% %         %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
% %     u2v1_temp = eleP_tau(2) * eleP_tau(1) * u_eleP_Pb' * bsxfun(@times,valueCoeff2.*phyGweights,u_eleP_Pb);
% %         %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
% %     u1v2_temp = eleP_tau(1) * eleP_tau(2) * u_eleP_Pb' * bsxfun(@times,valueCoeff1.*phyGweights,u_eleP_Pb);
% %         %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
% %     u2v2_temp = eleP_tau(2) * eleP_tau(2) * u_eleP_Pb' * bsxfun(@times,valueCoeff2.*phyGweights,u_eleP_Pb);
% %         %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
% %     
% %     interface_u1v1 = interface_u1v1 + sparse(uv_Row(:), uv_Col(:), u1v1_temp(:), Nelems*NTbases_u, Nelems*NTbases_u);
% %     interface_u2v1 = interface_u2v1 + sparse(uv_Row(:), uv_Col(:), u2v1_temp(:), Nelems*NTbases_u, Nelems*NTbases_u);
% %     interface_u1v2 = interface_u1v2 + sparse(uv_Row(:), uv_Col(:), u1v2_temp(:), Nelems*NTbases_u, Nelems*NTbases_u);
% %     interface_u2v2 = interface_u2v2 + sparse(uv_Row(:), uv_Col(:), u2v2_temp(:), Nelems*NTbases_u, Nelems*NTbases_u);
% %     
% % end % for n


% output matrix setting
interface_u1v1_test = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
interface_u2v1_test = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
interface_u1v2_test = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
interface_u2v2_test = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
Ninterface_test = interfacemeshInfo.Ninterface;
for jj = 1:Ninterface_test
    ePoint1 = interfacemeshInfo.interfaceNode(interfacemeshInfo.interfaceEdge(jj,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = interfacemeshInfo.interfaceNode(interfacemeshInfo.interfaceEdge(jj,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = sqrt( sum((ePoint1-ePoint2).^2) );
    
    phyGpointsX = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = areaEdge * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
        
    elemS = interfacemeshInfo.interface2StokesElem(jj,1);
    eIndx_elemS = interfacemeshInfo.interface2StokesElem(jj,2);
        %> i.e., [n, e] = interface2StokesElem(ii,:), n standsfor the ii-the interface belongs to n-th elem, 
        %> and e standsfor the ii-th interface is the e-th edge of n-th elem. 

    elemS_nu = StokesmeshInfo.nuEdge0Elem{elemS}(:,eIndx_elemS); 
    elemS_tau = [-elemS_nu(2),elemS_nu(1)];
    
    elemS_xT = StokesmeshInfo.centroidElem(elemS,1);  
    elemS_yT = StokesmeshInfo.centroidElem(elemS,2);
    elemS_hT = StokesmeshInfo.hElem(elemS);
        %> the centroid(xing xin) and diameter of elemS.
    
    %-- get the value of Coeffs func at Gauss Points
    valueCoeff1 = Coeff1(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valueCoeff2 = Coeff2(phyGpointsX, phyGpointsY); % [Npoints x 1]
    
    %-- subPart -- get the bases on quad
    [u_eleP_Pb, ~, ~] = ...
        localBases2D(elemS_xT, elemS_yT, elemS_hT, phyGpointsX, phyGpointsY, degreek_u);
        %> u_eleP_Pb, u_eleP_Pbx, u_eleP_Pby, [Npoints x NTbases_u]
        
        
    % get the position in the big mat
    [uv_Row, uv_Col] = getROWandCOL(elemS, 0, NTbases_u, NTbases_u);
    
    u1v1_temp = elemS_tau(1) * elemS_tau(1) * u_eleP_Pb' * bsxfun(@times,valueCoeff1.*phyGweights,u_eleP_Pb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    u2v1_temp = elemS_tau(2) * elemS_tau(1) * u_eleP_Pb' * bsxfun(@times,valueCoeff2.*phyGweights,u_eleP_Pb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    u1v2_temp = elemS_tau(1) * elemS_tau(2) * u_eleP_Pb' * bsxfun(@times,valueCoeff1.*phyGweights,u_eleP_Pb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    u2v2_temp = elemS_tau(2) * elemS_tau(2) * u_eleP_Pb' * bsxfun(@times,valueCoeff2.*phyGweights,u_eleP_Pb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    interface_u1v1_test = interface_u1v1_test + sparse(uv_Row(:), uv_Col(:), u1v1_temp(:), Nelems*NTbases_u, Nelems*NTbases_u);
    interface_u2v1_test = interface_u2v1_test + sparse(uv_Row(:), uv_Col(:), u2v1_temp(:), Nelems*NTbases_u, Nelems*NTbases_u);
    interface_u1v2_test = interface_u1v2_test + sparse(uv_Row(:), uv_Col(:), u1v2_temp(:), Nelems*NTbases_u, Nelems*NTbases_u);
    interface_u2v2_test = interface_u2v2_test + sparse(uv_Row(:), uv_Col(:), u2v2_temp(:), Nelems*NTbases_u, Nelems*NTbases_u);
    

end % for jj

% aaa = interface_u1v1_test - interface_u1v1;

end % function



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
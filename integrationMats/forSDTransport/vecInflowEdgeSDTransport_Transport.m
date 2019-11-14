function C_K_inflow = ...
    vecInflowEdgeSDTransport_Transport(t, pde, ...
    DarcymeshInfo, StokesmeshInfo, ...
    lastStep_DarcyU, lastStep_StokesU, ...
    formulaGauss1D, degreek_Darcy, degreek_Stokes, degreek_testC)
%
%   %--------------------------------------------------
%       Just copy from vecInflowEdgeTransport.m
%   %--------------------------------------------------
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
%       C_vecInflowEdge, convection term inflow edges, [Nelems*NTbases_test x NTelems*NTbases_trial], 
%  
%
%   YcZhang 24/10/2017
%
%   Last modified 24/10/2017
%
%

if isempty(lastStep_DarcyU) || isempty(lastStep_StokesU)
    error('in vecInflowEdgeSDTransport_Transport.m, lastStep_ is []')
end

% mesh and dof information: outflow edges 
D_Nelems = DarcymeshInfo.Nelems;
S_Nelems = StokesmeshInfo.Nelems;
S_Nbases = (degreek_Stokes+1)*(degreek_Stokes+2)/2;
D_Nbases = (degreek_Darcy+1)*(degreek_Darcy+2)/2;
C_Nbases_test = (degreek_testC+1)*(degreek_testC+2)/2;

dof_S_u1 = S_Nelems*S_Nbases;

Nelems = D_Nelems +S_Nelems;

% mesh information: inflow edges 
if isfield(StokesmeshInfo,'inflowEdgeIndex')
    inflowEdgeIndex_S = StokesmeshInfo.inflowEdgeIndex;
else
    inflowEdgeIndex_S =[];
end 
if isfield(DarcymeshInfo,'inflowEdgeIndex')
    inflowEdgeIndex_D = DarcymeshInfo.inflowEdgeIndex;
else
    inflowEdgeIndex_D =[];
end 


C_K_inflow = zeros(Nelems*C_Nbases_test,1);

if ~isempty(inflowEdgeIndex_D)
    for CurrEdge = 1:length(inflowEdgeIndex_D)
        %% Part I, get the information about 
        % 1. physical GaussPoints, 
        % 2. different element bases on phy GaussPoints on ii-th edge.
        %
        %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
        inflowEindx = inflowEdgeIndex_D(CurrEdge); 
        ePoint1 = DarcymeshInfo.node(DarcymeshInfo.edge(inflowEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
        ePoint2 = DarcymeshInfo.node(DarcymeshInfo.edge(inflowEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
        areaEdge = DarcymeshInfo.areaEdge(inflowEindx);

        phyGpointsX = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
            %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
        phyGpointsY = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
            %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
        phyGweights = areaEdge * formulaGauss1D(:,2);
            %> [Npoints x 1], the physical Gauss weights on ii-th edge.

        elem1 = DarcymeshInfo.edge2elem(inflowEindx,1);
        elem2 = DarcymeshInfo.edge2elem(inflowEindx,2);
        local_e1 = DarcymeshInfo.edge2elem(inflowEindx,3);
        local_e2 = DarcymeshInfo.edge2elem(inflowEindx,4);
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

        %-- subPart -- get the bases on quad    
        [elemP_testPb, ~, ~] = ...
            localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_testC);
            %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]

        %<<-- End Part I ---------------------------------------------------------------------------------

    %% ------------------------------------------ Dividing line -------------------------------------------------- %%
    %-------------------------------------------------------------------------------------------------------------------%

        %% Part II, compute the inflow edges intgrations.
        %
        %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
        %
        %-- the inflow edges
        % get the position in the big mat    
        [~, elemPtrialPbx_D, elemPtrialPby_D] = localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_Darcy);
        elemPUh = lastStep_DarcyU((elemP-1)*D_Nbases+1:elemP*D_Nbases);
        value_vector_u_1 = elemPtrialPbx_D*elemPUh;
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
        value_vector_u_2 = elemPtrialPby_D*elemPUh;
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
        
        value_f_inflow = pde.f_inflow(t, phyGpointsX, phyGpointsY);
        
        Rhs_inflow_temp = elemP_nu(1) * elemP_testPb' * bsxfun(@times,value_vector_u_1.*phyGweights,value_f_inflow) ...
            + elemP_nu(2) * elemP_testPb' * bsxfun(@times,value_vector_u_2.*phyGweights,value_f_inflow);
        
        rhsInflowRow = ((elemP-1)*C_Nbases_test+1 : (elemP)*C_Nbases_test);
        
        C_K_inflow(rhsInflowRow,1) = C_K_inflow(rhsInflowRow,1) + Rhs_inflow_temp;

    end % for ii
end % if ~isempty(inflowEdgeIndex_D)


if ~isempty(inflowEdgeIndex_S)
    for CurrEdge = 1:length(inflowEdgeIndex_S)
        %% Part I, get the information about 
        % 1. physical GaussPoints, 
        % 2. different element bases on phy GaussPoints on ii-th edge.
        %
        %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
        inflowEindx = inflowEdgeIndex_S(CurrEdge); 
        ePoint1 = StokesmeshInfo.node(StokesmeshInfo.edge(inflowEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
        ePoint2 = StokesmeshInfo.node(StokesmeshInfo.edge(inflowEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
        areaEdge = StokesmeshInfo.areaEdge(inflowEindx);

        phyGpointsX = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
            %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
        phyGpointsY = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
            %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
        phyGweights = areaEdge * formulaGauss1D(:,2);
            %> [Npoints x 1], the physical Gauss weights on ii-th edge.

        elem1 = StokesmeshInfo.edge2elem(inflowEindx,1);
        elem2 = StokesmeshInfo.edge2elem(inflowEindx,2);
        local_e1 = StokesmeshInfo.edge2elem(inflowEindx,3);
        local_e2 = StokesmeshInfo.edge2elem(inflowEindx,4);
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

        %-- subPart -- get the bases on quad    
        [elemP_testPb, ~, ~] = ...
            localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_testC);
            %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]

        %<<-- End Part I ---------------------------------------------------------------------------------

    %% ------------------------------------------ Dividing line -------------------------------------------------- %%
    %-------------------------------------------------------------------------------------------------------------------%

        %% Part II, compute the inflow edges intgrations.
        %
        %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
        %
        %-- the inflow edges
        % get the position in the big mat    
        [elemPtrialPb_S, ~, ~] = localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_Darcy);
        elemPUh1 = lastStep_StokesU((elemP-1)*S_Nbases+1:elemP*S_Nbases);
        elemPUh2 = lastStep_StokesU(dof_S_u1+(elemP-1)*S_Nbases+1:dof_S_u1+elemP*S_Nbases);
        
        value_vector_u_1 = elemPtrialPb_S*elemPUh1;
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
        value_vector_u_2 = elemPtrialPb_S*elemPUh2;
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
        
        value_f_inflow = pde.f_inflow(t, phyGpointsX, phyGpointsY);
        
        Rhs_inflow_temp = elemP_nu(1) * elemP_testPb' * bsxfun(@times,value_vector_u_1.*phyGweights,value_f_inflow) ...
        + elemP_nu(2) * elemP_testPb' * bsxfun(@times,value_vector_u_2.*phyGweights,value_f_inflow);
        
        rhsInflowRow = ((D_Nelems+elemP-1)*C_Nbases_test+1 : (D_Nelems+elemP)*C_Nbases_test);
        
        C_K_inflow(rhsInflowRow,1) = C_K_inflow(rhsInflowRow,1) + Rhs_inflow_temp;

    end % for ii
end % if ~isempty(inflowEdgeIndex_S)

end % function
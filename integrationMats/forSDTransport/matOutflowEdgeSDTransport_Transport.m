function [C_K_outflow, DSC_K_outflow, C_K_outflow_RhsLastStep] = ...
    matOutflowEdgeSDTransport_Transport(... %
    DarcymeshInfo, StokesmeshInfo, ...
    lastStep_DarcyU, lastStep_StokesU, lastStep_C, ...
    formulaGauss1D, degreek_Darcy, degreek_Stokes, degreek_trialC, degreek_testC )
%
%   %-----------------------------------------------------------------
%       Just copy from matOutflowEdgeTransport.m
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
%                   Coeffs{3}, Coeffs{4} are the vector_u function vector_u_1, vector_u_2.
%                   let k12=0, k21=0, the Coeffs{1}=k11, Coeffs{2}= k22. 
%       meshInfo, the mesh information.
%       option, particularly, we need the \epsilon and the penalty parameter.
%       formulaGauss1D, the 1d Gauss quadrature formula, size: a matrix, [Npoints x 2],
%               the first column is the all Gauss-Points on [0,1],
%               the second is the weights of all Gauss-Points.
%       degreek_trial, the trial polynomial degree k.
%       degreek_test, the test polynomial degree k.
%
%   output:
%       C_matInflowEdge, convection term outflow edges, [Nelems*NTbases_test x NTelems*NTbases_trial], 
%  
%
%   YcZhang  22/10/2017
%
%   Last modified 23/10/2017
%
%

% mesh and dof information: outflow edges 
D_Nelems = DarcymeshInfo.Nelems;
S_Nelems = StokesmeshInfo.Nelems;
S_Nbases = (degreek_Stokes+1)*(degreek_Stokes+2)/2;
D_Nbases = (degreek_Darcy+1)*(degreek_Darcy+2)/2;

dof_S_u1 = S_Nelems*S_Nbases;
dof_S_p = S_Nelems*((degreek_Stokes-1+1)*(degreek_Stokes-1+2)/2);
dof_D_u = D_Nelems*D_Nbases;

dof_SD = (dof_D_u+2*dof_S_u1+dof_S_p);

Nelems = D_Nelems +S_Nelems;

if isfield(StokesmeshInfo,'outflowEdgeIndex')
    outflowEdgeIndex_S = StokesmeshInfo.outflowEdgeIndex;
else
    outflowEdgeIndex_S =[];
end 
if isfield(DarcymeshInfo,'outflowEdgeIndex')
    outflowEdgeIndex_D = DarcymeshInfo.outflowEdgeIndex;
else
    outflowEdgeIndex_D =[];
end 


% setting
C_Nbases_trial = (degreek_trialC+1)*(degreek_trialC+2)/2;
C_Nbases_test = (degreek_testC+1)*(degreek_testC+2)/2;

C_K_outflow = sparse(Nelems*C_Nbases_test,Nelems*C_Nbases_trial);
DC_K_outflow = sparse(Nelems*C_Nbases_test,dof_SD);
SC_K_outflow = sparse(Nelems*C_Nbases_test,dof_SD);
% % DSC_K_outflow = sparse(Nelems*C_Nbases_test,dof_SD);
C_K_outflow_RhsLastStep = zeros(Nelems*C_Nbases_test,1);

if ~isempty(outflowEdgeIndex_D)
    for CurrEdge = 1:length(outflowEdgeIndex_D)
        %% Part I, get the information about 
        % 1. physical GaussPoints, 
        % 2. different element bases on phy GaussPoints on ii-th edge.
        %
        %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
        outflowEindx = outflowEdgeIndex_D(CurrEdge); 
        ePoint1 = DarcymeshInfo.node(DarcymeshInfo.edge(outflowEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
        ePoint2 = DarcymeshInfo.node(DarcymeshInfo.edge(outflowEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
        areaEdge = DarcymeshInfo.areaEdge(outflowEindx);

        phyGpointsX = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
            %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
        phyGpointsY = DarcymeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
            %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
        phyGweights = areaEdge * formulaGauss1D(:,2);
            %> [Npoints x 1], the physical Gauss weights on ii-th edge.

        elem1 = DarcymeshInfo.edge2elem(outflowEindx,1);
        elem2 = DarcymeshInfo.edge2elem(outflowEindx,2);
        local_e1 = DarcymeshInfo.edge2elem(outflowEindx,3);
        local_e2 = DarcymeshInfo.edge2elem(outflowEindx,4);
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
        [elemP_trialPb, ~, ~] = ...
            localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_trialC);
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
        %-- the outflow edges
        % get the position in the big mat    
        [C_DarcyRow, C_DarcyCol] = getROWandCOL(elemP, 0, C_Nbases_trial, C_Nbases_test);

        %--------------------- term ----------------------------------
        [~, elemPtrialPbx_D, elemPtrialPby_D] = localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_Darcy);
        if ~isempty(lastStep_DarcyU)
            % in the Darcy domain
            elemPUh = lastStep_DarcyU((elemP-1)*D_Nbases+1:elemP*D_Nbases);
            
            value_vector_u_1 = elemPtrialPbx_D*elemPUh;
                %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
            value_vector_u_2 = elemPtrialPby_D*elemPUh;
                %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
                
            C_K_outflow_temp = elemP_nu(1) * elemP_testPb' * bsxfun(@times,value_vector_u_1.*phyGweights,elemP_trialPb) ...
                + elemP_nu(2) * elemP_testPb' * bsxfun(@times,value_vector_u_2.*phyGweights,elemP_trialPb);
                %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
            C_K_outflow = C_K_outflow ...
                + sparse(C_DarcyRow(:), C_DarcyCol(:), ...
                C_K_outflow_temp(:), ...
                Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
        else
            error('in matOutflowEdgeSDTransport_Transport.m, lastStep_DarcyU is []')
        end
        
        %---------------------- term --------------------------------------
        if ~isempty(lastStep_C)
            elemP_Ch = lastStep_C((elemP-1)*C_Nbases_trial+1 : elemP*C_Nbases_trial);
            value_elemCh = elemP_trialPb*elemP_Ch; % [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
            DC_K_outflow_temp = elemP_nu(1) * elemP_testPb' * bsxfun(@times,value_elemCh.*phyGweights,elemPtrialPbx_D) ...
                + elemP_nu(2) * elemP_testPb' * bsxfun(@times,value_elemCh.*phyGweights,elemPtrialPby_D);
            
            [DC_Row, DC_COL] = getROWandCOL(elemP, 0, D_Nbases, C_Nbases_test);
            
            DC_K_outflow = DC_K_outflow ...
                + sparse(DC_Row(:), DC_COL(:), ...
                DC_K_outflow_temp(:), ...
                Nelems*C_Nbases_test,dof_SD);
        else
            error('in matOutflowEdgeSDTransport_Transport.m, lastStep_C is []')
        end
        
        
        %----------------------- RhsLastStep --------------------------------
        elemP_Row_vec = ((elemP-1)*C_Nbases_test+1 : (elemP)*C_Nbases_test);

        C_K_outflow_RhsLastStep(elemP_Row_vec,1) = C_K_outflow_RhsLastStep(elemP_Row_vec,1) ...
            + elemP_nu(1)*elemP_testPb'*(value_elemCh.*value_vector_u_1.*phyGweights) ...
            + elemP_nu(2)*elemP_testPb'*(value_elemCh.*value_vector_u_2.*phyGweights);
        %--------------------------------------------------------------------------

    end % for ii
end % if ~isempty(outflowEdgeIndex_D)

if ~isempty(outflowEdgeIndex_S)
    for CurrEdge = 1:length(outflowEdgeIndex_S)
        %% Part I, get the information about 
        % 1. physical GaussPoints, 
        % 2. different element bases on phy GaussPoints on ii-th edge.
        %
        %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
        outflowEindx = outflowEdgeIndex_S(CurrEdge); 
        ePoint1 = StokesmeshInfo.node(StokesmeshInfo.edge(outflowEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
        ePoint2 = StokesmeshInfo.node(StokesmeshInfo.edge(outflowEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
        areaEdge = StokesmeshInfo.areaEdge(outflowEindx);

        phyGpointsX = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
            %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
        phyGpointsY = StokesmeshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
            %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
        phyGweights = areaEdge * formulaGauss1D(:,2);
            %> [Npoints x 1], the physical Gauss weights on ii-th edge.

        elem1 = StokesmeshInfo.edge2elem(outflowEindx,1);
        elem2 = StokesmeshInfo.edge2elem(outflowEindx,2);
        local_e1 = StokesmeshInfo.edge2elem(outflowEindx,3);
        local_e2 = StokesmeshInfo.edge2elem(outflowEindx,4);
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
        [elemP_trialPb, ~, ~] = ...
            localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_trialC);
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
        %-- the outflow edges
        % get the position in the big mat    
        [C_StokesRow, C_StokesCol] = getROWandCOL(D_Nelems+elemP, 0, C_Nbases_trial, C_Nbases_test);

        %--------------------- term ----------------------------------
        [elemPtrialPb_S, ~, ~] = localBases2D(elemP_xT, elemP_yT, elemP_hT, phyGpointsX, phyGpointsY, degreek_Stokes);
        if ~isempty(lastStep_StokesU)
            % in the Darcy domain
            elemPUh1 = lastStep_StokesU((elemP-1)*S_Nbases+1:elemP*S_Nbases);
            elemPUh2 = lastStep_StokesU(dof_S_u1+(elemP-1)*S_Nbases+1:dof_S_u1+elemP*S_Nbases);

            value_vector_u_1 = elemPtrialPb_S*elemPUh1;
                %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
            value_vector_u_2 = elemPtrialPb_S*elemPUh2;
                %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
                
            C_K_outflow_temp = elemP_nu(1) * elemP_testPb' * bsxfun(@times,value_vector_u_1.*phyGweights,elemP_trialPb) ...
                + elemP_nu(2) * elemP_testPb' * bsxfun(@times,value_vector_u_2.*phyGweights,elemP_trialPb);
                %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
            C_K_outflow = C_K_outflow ...
                + sparse(C_StokesRow(:), C_StokesCol(:), ...
                C_K_outflow_temp(:), ...
                Nelems*C_Nbases_test, Nelems*C_Nbases_trial);
        else
            error('in matOutflowEdgeSDTransport_Transport.m, lastStep_DarcyU is []')
        end
        
        %---------------------- term --------------------------------------
        if ~isempty(lastStep_C)
            elemP_Ch = lastStep_C((D_Nelems+elemP-1)*C_Nbases_trial+1 : (D_Nelems+elemP)*C_Nbases_trial);
            value_elemCh = elemP_trialPb*elemP_Ch; % [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
            
            SC_K_outflow_temp_u1 = elemP_nu(1) * elemP_testPb' * bsxfun(@times,value_elemCh.*phyGweights,elemPtrialPb_S);
            SC_K_outflow_temp_u2 = elemP_nu(2) * elemP_testPb' * bsxfun(@times,value_elemCh.*phyGweights,elemPtrialPb_S);
            
            
            SC_Row_u1 = ((D_Nelems+elemP-1)*C_Nbases_test+1 : (D_Nelems+elemP)*C_Nbases_test)' ... 
                *ones(1,C_Nbases_test); % [NTbases_test x NTbases_trial]
            SC_COL_u1 = ones(S_Nbases, 1) ...
                *(dof_D_u + (elemP-1)*S_Nbases+1 : dof_D_u + elemP*S_Nbases); % [NTbases_test x NTbases_trial]
            SC_Row_u2 = SC_Row_u1;
            SC_COL_u2 = dof_S_u1 + SC_COL_u1;
            
            SC_K_outflow = SC_K_outflow ...
                + sparse(SC_Row_u1(:), SC_COL_u1(:), ...
                SC_K_outflow_temp_u1(:), ...
                Nelems*C_Nbases_test,dof_SD);
            SC_K_outflow = SC_K_outflow ...
                + sparse(SC_Row_u2(:), SC_COL_u2(:), ...
                SC_K_outflow_temp_u2(:), ...
                Nelems*C_Nbases_test,dof_SD);
        else
            error('in matOutflowEdgeSDTransport_Transport.m, lastStep_C is []')
        end
        
        
        %----------------------- RhsLastStep --------------------------------
        elemP_Row_vec = ((D_Nelems+elemP-1)*C_Nbases_test+1 : (D_Nelems+elemP)*C_Nbases_test);

        C_K_outflow_RhsLastStep(elemP_Row_vec,1) = C_K_outflow_RhsLastStep(elemP_Row_vec,1) ...
            + elemP_nu(1)*elemP_testPb'*(value_elemCh.*value_vector_u_1.*phyGweights) ...
            + elemP_nu(2)*elemP_testPb'*(value_elemCh.*value_vector_u_2.*phyGweights);
        %--------------------------------------------------------------------------

    end % for ii
end % if ~isempty(outflowEdgeIndex_S)

DSC_K_outflow = DC_K_outflow + SC_K_outflow;

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

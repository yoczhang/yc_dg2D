%%>>-- Begin sub function 1 ---------------------------------------------------------
function [StokesmeshInfo, DarcymeshInfo] = ...
    getTransportBoundaryInfo2(pde, StokesmeshInfo, DarcymeshInfo, ...
    lastStep_DarcyU, lastStep_StokesU, degreek_Darcy, degreek_Stokes)
%
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       dgDecoupledSDTransport2, 
%       for the transport eqn, have inflow and outflow boundaryEdges.
%   %---------------------------------------------
%
%
%
%
%   YcZhang 28/10/2017
%
%   Last modified 28/10/2017
%

c_darcy = -pde.K;

S_Nelems = StokesmeshInfo.Nelems;
S_Nbases = (degreek_Stokes+1)*(degreek_Stokes+2)/2;
D_Nbases = (degreek_Darcy+1)*(degreek_Darcy+2)/2;

D_DirichletEdges = DarcymeshInfo.DirichletEdgeIndex; % here we need the DirichletEdges information.
S_DirichletEdges = StokesmeshInfo.DirichletEdgeIndex; % here we need the DirichletEdges information.

D_inflowEdgeIndex = [];
S_inflowEdgeIndex = [];
D_outflowEdgeIndex = [];
S_outflowEdgeIndex = [];

for CurrEdge = 1:length(D_DirichletEdges)
    DirichletEindx = D_DirichletEdges(CurrEdge); 
    ePoint1 = DarcymeshInfo.node(DarcymeshInfo.edge(DirichletEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = DarcymeshInfo.node(DarcymeshInfo.edge(DirichletEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    midPoint = (ePoint1 + ePoint2)/2;
    
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
        
        
    %-- subPart -- get the bases on quad
    [~, elemP_trialPbx, elemP_trialPby] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, midPoint(1), midPoint(2), degreek_Darcy);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [1 x NTg1bases_trial]
    
    elemPUh = lastStep_DarcyU((elemP-1)*D_Nbases+1:elemP*D_Nbases);
    value_vector_u_1 = c_darcy*elemP_trialPbx*elemPUh;
        %> [1 x Nbases]*[Nbases x 1] --> [1 x 1]
	value_vector_u_2 = c_darcy*elemP_trialPby*elemPUh;
        %> [1 x Nbases]*[Nbases x 1] --> [1 x 1]
    
    if(value_vector_u_1*elemP_nu(1) + value_vector_u_2*elemP_nu(2) < 0)
        D_inflowEdgeIdex_temp = [D_inflowEdgeIndex; DirichletEindx];
        D_inflowEdgeIndex = D_inflowEdgeIdex_temp;
    else
        D_outflowEdgeIdex_temp = [D_outflowEdgeIndex; DirichletEindx];
        D_outflowEdgeIndex = D_outflowEdgeIdex_temp;
    end 
    
end % for CurrE = 1:length(D_DirichletEdges)


for CurrEdge = 1:length(S_DirichletEdges)
    DirichletEindx = S_DirichletEdges(CurrEdge); 
    ePoint1 = StokesmeshInfo.node(StokesmeshInfo.edge(DirichletEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = StokesmeshInfo.node(StokesmeshInfo.edge(DirichletEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    midPoint = (ePoint1 + ePoint2)/2;
    
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
        
        
    %-- subPart -- get the bases on quad
    [elemP_trialPb, ~, ~] = ...
        localBases2D(elemP_xT, elemP_yT, elemP_hT, midPoint(1), midPoint(2), degreek_Stokes);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [1 x NTg1bases_trial]
    
    elemPUh1 = lastStep_StokesU((elemP-1)*S_Nbases+1:elemP*S_Nbases);
    elemPUh2 = lastStep_StokesU((S_Nelems+elemP-1)*S_Nbases+1:(S_Nelems+elemP)*S_Nbases);
    value_vector_u_1 = elemP_trialPb*elemPUh1;
        %> [1 x Nbases]*[Nbases x 1] --> [1 x 1]
	value_vector_u_2 = elemP_trialPb*elemPUh2;
        %> [1 x Nbases]*[Nbases x 1] --> [1 x 1]
    
    if(value_vector_u_1*elemP_nu(1) + value_vector_u_2*elemP_nu(2) < 0)
        S_inflowEdgeIdex_temp = [S_inflowEdgeIndex; DirichletEindx];
        S_inflowEdgeIndex = S_inflowEdgeIdex_temp;
    else
        S_outflowEdgeIdex_temp = [S_outflowEdgeIndex; DirichletEindx];
        S_outflowEdgeIndex = S_outflowEdgeIdex_temp;
    end 
    
end % for CurrE = 1:length(S_DirichletEdges)

DarcymeshInfo.inflowEdgeIndex = D_inflowEdgeIndex;
DarcymeshInfo.outflowEdgeIndex = D_outflowEdgeIndex;

StokesmeshInfo.inflowEdgeIndex = S_inflowEdgeIndex;
StokesmeshInfo.outflowEdgeIndex = S_outflowEdgeIndex;


end % function
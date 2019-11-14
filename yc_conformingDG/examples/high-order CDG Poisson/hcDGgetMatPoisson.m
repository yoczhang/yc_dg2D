function [sysStiffMat, vecRhsF] = ...
    hcDGgetMatPoisson(func_f, meshInfo, Gaussformulas, basis_k)
%
%
%   YcZhang  Jul.08.2019
%
%   Last modified  Jul.08.2019
%
%

Nelems = meshInfo.Nelems;

Coeff_one = @(x,y) 1+0.*x;
% Coeff_k11 = Coeff_one;
% Coeff_k22 = Coeff_one;

Gaussformula2D = Gaussformulas{1};
Gaussformula1D = Gaussformulas{2};


% original stiff-matrix
NTbases = nchoosek(basis_k+2,basis_k);
vecRhsF = zeros(Nelems*NTbases,1);

% reconstruction stiff-matrix
plusOrder = 1;
Rbasis_k = basis_k+plusOrder;
RNTbases = nchoosek(Rbasis_k+2,Rbasis_k);
RRstiffMat = sparse(Nelems*RNTbases,Nelems*RNTbases);
mdRRstiffMat = sparse(Nelems*RNTbases,Nelems*RNTbases); % modified RstiffMat
RrhsMat = sparse(Nelems*RNTbases, Nelems*NTbases);


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

for CurrElem = 1:Nelems
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    singleElem = meshInfo.elem{CurrElem,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges on singleElem.
    
    concavePointElem_ii = meshInfo.concavePointElem(CurrElem); % [1 x 1], the concavePointElem of ii-the elem.
    coordv = meshInfo.node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    centroidElem = meshInfo.centroidElem(CurrElem,:); % [1 x 2], the centroid(xing xin) (x-coord, y-coord) of element.
    
    coordTri0Elem = getcoordTri0Elem(singleNE, concavePointElem_ii, centroidElem, coordv); 
        %> if singleNE==3, coordTri0Elem, [3 x 2]. else coordTri0Elem, [3*singleNE x 2].
    
    elem_xT = centroidElem(1);
    elem_yT = centroidElem(2); 
    elem_hT = meshInfo.hElem(CurrElem); 
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    
    %- the original-reconstruction stiff-matrix
    %ORstiffMat_elem = zeros(RNTbases, NTbases);
    ORstiffMat_elem_dxdx = zeros(RNTbases, NTbases);
    ORstiffMat_elem_dydy = zeros(RNTbases, NTbases);
    
    %- the reconstruction-reconstruction stiff-matrix
    %RstiffMat_elem = zeros(RNTbases, RNTbases);
    %mdRstiffMat_elem = zeros(RNTbases, RNTbases); % modified RstiffMat_elem
    RRstiffMat_elem_dxdx = zeros(RNTbases, RNTbases);
    RRstiffMat_elem_dydy = zeros(RNTbases, RNTbases);
    
    %- the vec term
    V1 = zeros(RNTbases,1); % the reconstruction vector
    V2 = zeros(NTbases,1); % the original vector
    
    %- the rhs term
    vecRhsF_elem = zeros(NTbases,1);
    
    %-- the integration on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints2D, phyGweights2D] = getGaussLocalTri(coordTri_nt, Gaussformula2D);

        % get the bases values on quad points
        [Pb, Pbx, Pby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints2D(:,1), phyGpoints2D(:,2), basis_k);
            %> RPb, RPbx, RPby, [Npoints x RNTbases].

        [RPb, RPbx, RPby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints2D(:,1), phyGpoints2D(:,2), Rbasis_k);
            %> RPb, RPbx, RPby, [Npoints x RNTbases].

        %--- the funcValue may be chosen by case.
        valueCoeffOne = Coeff_one(phyGpoints2D(:,1), phyGpoints2D(:,2));
        %valueCoeff_k11 = Coeff_k11(phyGpoints2D(:,1), phyGpoints2D(:,2));
        %valueCoeff_k22 = Coeff_k22(phyGpoints2D(:,1), phyGpoints2D(:,2));
        valuefunc_f = func_f(phyGpoints2D(:,1), phyGpoints2D(:,2));
        
        %- the original-reconstruction stiff-matrix
        ORstiffMat_elem_dxdx = ORstiffMat_elem_dxdx + getMatOnElem(valueCoeffOne, Pbx, RPbx, phyGweights2D);
        ORstiffMat_elem_dydy = ORstiffMat_elem_dydy + getMatOnElem(valueCoeffOne, Pby, RPby, phyGweights2D);
        
        %- the reconstruction-reconstruction stiff-matrix
        RRstiffMat_elem_dxdx = RRstiffMat_elem_dxdx + getMatOnElem(valueCoeffOne, RPbx, RPbx, phyGweights2D);
        RRstiffMat_elem_dydy = RRstiffMat_elem_dydy + getMatOnElem(valueCoeffOne, RPby, RPby, phyGweights2D);
        
        %- 
        V1 = V1 + RPb' * (phyGweights2D.*valueCoeffOne); % [NRTbases x 1]
        V2 = V2 + Pb' * (phyGweights2D.*valueCoeffOne); % [NTbases x 1]
        
        %-- VecF
        vecRhsF_elem = vecRhsF_elem + Pb' * ( phyGweights2D.*valuefunc_f );
        
    end % for nt
    
    ORstiffMat_elem = ORstiffMat_elem_dxdx + ORstiffMat_elem_dydy;
    ORstiffMat_elem(1,:) = V2';
    
    RRstiffMat_elem = RRstiffMat_elem_dxdx + RRstiffMat_elem_dydy;
    mdRRstiffMat_elem = RRstiffMat_elem;
    mdRRstiffMat_elem(1,:) = V1';
    
    
    %- get row and column
    ORRow_elem = ((CurrElem-1)*RNTbases+1 : CurrElem*RNTbases)' ...
        * ones(1,NTbases); % [NTbases_test x NTbases_trial]
    ORCol_elem = ones(RNTbases, 1) ...
        *((CurrElem-1)*NTbases+1 : CurrElem*NTbases); % [NTbases_test x NTbases_trial]
    RRRow_elem = ((CurrElem-1)*RNTbases+1 : CurrElem*RNTbases)' ...
        * ones(1,RNTbases); % [NTbases_test x NTbases_trial]
    RRCol_elem = ones(RNTbases, 1) ...
        *((CurrElem-1)*RNTbases+1 : CurrElem*RNTbases); % [NTbases_test x NTbases_trial]
    
    
    
    %- original-reconstruction edge terms
%     ORMat_edge = zeros(RNTbases, NTbases);
%     Nei_ORMat_edge = zeros(RNTbases, NTbases); % Neighbour
    
    %-- the integration on edges
    for CurrEdge = 1:singleNE
        %--- edge setting
        edgeIndex = meshInfo.elem2edge{CurrElem}(CurrEdge);
        ePoint1 = meshInfo.node(meshInfo.edge(edgeIndex,1),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint1
        ePoint2 = meshInfo.node(meshInfo.edge(edgeIndex,2),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint2
        %edge_xE = (ePoint1(1) + ePoint2(1))/2;
        %edge_yE = (ePoint1(2) + ePoint2(2))/2;
        edge_hE = meshInfo.areaEdge(edgeIndex);
        edge_nu = meshInfo.nuEdge0Elem{CurrElem}(:,CurrEdge); % [1 x 2], the outward unit normal
        
        %- get neighbour element info
        edge2elem = meshInfo.edge2elem(edgeIndex,1:2);
        neiElemIndx = setdiff(edge2elem,CurrElem);
        
        %--- 1D Gauss Points
        phyGpoints1DX = (ePoint1(1)-ePoint2(1))/2*Gaussformula1D(:,1) + (ePoint1(1)+ePoint2(1))/2; 
        phyGpoints1DY = (ePoint1(2)-ePoint2(2))/2*Gaussformula1D(:,1) + (ePoint1(2)+ePoint2(2))/2; 
        phyGweights1D = edge_hE*Gaussformula1D(:,2)/2;
        
        %--- some values at 1D GaussPoints
        [~, RPbx, RPby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, Rbasis_k);
            %> RT_Pb, RT_Pbx, RT_Pby, [NFpoints x NRTbases].
        [Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, basis_k);
            %> T_Pb, T_Pbx, T_Pby, [NFpoints x NTbases].    
        valueCoeffOne_1D = Coeff_one(phyGpoints1DX, phyGpoints1DY); % [NFpoints x 1].
        
        %--- integration matrix
        ORMat_edge = edge_nu(1)*getMatOnElem(valueCoeffOne_1D, Pb, RPbx, phyGweights1D) ...
            + edge_nu(2)*getMatOnElem(valueCoeffOne_1D, Pb, RPby, phyGweights1D); % [NRTbases x NFbases]
        
        DirIndictor = -1;
        
        %- neighbour integration matrix
         if neiElemIndx > 0
            Nei_centroidElem = meshInfo.centroidElem(neiElemIndx,:); % [1 x 2], the centroid(xing xin) (x-coord, y-coord) of element.
            Nei_elem_xT = Nei_centroidElem(1);
            Nei_elem_yT = Nei_centroidElem(2); 
            Nei_elem_hT = meshInfo.hElem(neiElemIndx); 
            
            [Nei_Pb, ~, ~] = localBases2D(Nei_elem_xT, Nei_elem_yT, Nei_elem_hT, phyGpoints1DX, phyGpoints1DY, basis_k);
                %> T_Pb, T_Pbx, T_Pby, [NFpoints x NTbases].
                
            Nei_ORMat_edge = edge_nu(1)*getMatOnElem(valueCoeffOne_1D, Nei_Pb, RPbx, phyGweights1D) ...
                + edge_nu(2)*getMatOnElem(valueCoeffOne_1D, Nei_Pb, RPby, phyGweights1D); % [NRTbases x NTbases]
         
            Nei_ORCol_elem = ones(RNTbases, 1) ...
                *((neiElemIndx-1)*NTbases+1 : neiElemIndx*NTbases); % [NTg1bases_test x NTg1bases_trial]
    
            RrhsMat = RrhsMat ...                
                + sparse(ORRow_elem(:), Nei_ORCol_elem(:), ...
                (1/2)*Nei_ORMat_edge(:), ...
                Nelems*RNTbases, Nelems*NTbases);
            
            %- update DirIndictor
            DirIndictor = -1/2;
         end
         
         RrhsMat = RrhsMat ...                
             + sparse(ORRow_elem(:), ORCol_elem(:), ...
             DirIndictor*ORMat_edge(:), ...
             Nelems*RNTbases, Nelems*NTbases);
        
    end % CurrEdge
    
    %-- add MatElem to the position in the big mat
    mdRRstiffMat = mdRRstiffMat ...
        + sparse(RRRow_elem(:), RRCol_elem(:), ...
        mdRRstiffMat_elem(:), ...
        Nelems*RNTbases, Nelems*RNTbases);
    
    RRstiffMat = RRstiffMat ...
        + sparse(RRRow_elem(:), RRCol_elem(:), ...
        RRstiffMat_elem(:), ...
        Nelems*RNTbases, Nelems*RNTbases);
    
    RrhsMat = RrhsMat ...                
        + sparse(ORRow_elem(:), ORCol_elem(:), ...
        ORstiffMat_elem(:), ...
        Nelems*RNTbases, Nelems*NTbases);
    
    
    %-- value_rhs_f -- add to the position in the big vector
    Row_vec = ((CurrElem-1)*NTbases+1 : CurrElem*NTbases);
    vecRhsF(Row_vec,1) = vecRhsF(Row_vec,1) + vecRhsF_elem;
    
    %<<-- End Part II --------------------------------------------------------------------------------------
    
end % for ii


Rmat = mdRRstiffMat\RrhsMat;

sysStiffMat = Rmat' * RRstiffMat * Rmat;


end % function 



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
    
end % if n==1

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

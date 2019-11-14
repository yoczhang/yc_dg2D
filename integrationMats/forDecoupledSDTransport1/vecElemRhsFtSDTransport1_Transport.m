function vecRhsFt = vecElemRhsFtSDTransport1_Transport(t, ...
    pde, lastStep_DarcyUh, lastStep_StokesUh,...
    DarcymeshInfo, StokesmeshInfo, ...
    formulaGauss2D, degreek_Darcy, degreek_Stokes, degreek_C)
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       dgDecoupledSDTransport1, 
%       for the transport eqn, have NO inflow and outflow boundaryEdges,
%       i.e. for the transport eqn, all the boundaryEdges are Dirichlet Edges.
%   %---------------------------------------------
%
%   %--------------------------------------------------
%       Just copy from vecElemRhsFtTransport.m
%   %--------------------------------------------------
%
%   In this function we compute the integral f(t,x,y) on element T.
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases_trial denote the number of LOCAL trial bases on each K of Th.
%               NTbases_test denote the number of LOCAL test bases on each K of Th.
%
%   input:
%       rhs_f_func, the cell-type, here
%                   rhs_f_func{i} is the first rhs f func.
%       meshInfos, the mesh information.
%       formulaGauss2D, the 1d Gauss quadrature formula, size: a matrix, [Npoints x 3],
%               the first column is the x-coordinates of all Gauss-Points,
%               the second column is the y-coordinates of all Gauss-Points,
%               the third is the weights of all Gauss-Points.
%       degreek_test, the group1 test polynomial degree k.
%           %> here we consider that, in the Stokes equation, there will have VELOCITY u and PRESSURE p, 
%           %> so there have \int_e uv ds, then g1Trial_k=u_degree_k, g1Test_k=u_degree_k. 
%           %> And also have \int_e pv ds, then g2Trial_k=p_degree_k, g2Test_k=u_degree_k. 
%
%   output:
%       vecRhsFt, [Nelems*NTg1bases_test x 1].
%
%
%   YcZhang 27/10/2017
%
%   Last modified 27/10/2017
%
%

if isempty(lastStep_DarcyUh)
    error('in vecElemRhsFtSDTransport1_Transport.m, lastStep_DarcyUh is []')
end
if isempty(lastStep_StokesUh)
    error('in vecElemRhsFtSDTransport1_Transport.m, lastStep_StokesUh is []')
end

c_darcy = -pde.K;

D_Nelems = DarcymeshInfo.Nelems;
S_Nelems = StokesmeshInfo.Nelems;

Nelems = D_Nelems + S_Nelems;

% bases setting
D_Nbases = (degreek_Darcy+1)*(degreek_Darcy+2)/2;
S_Nbases = (degreek_Stokes+1)*(degreek_Stokes+2)/2;
C_Nbases = (degreek_C+1)*(degreek_C+2)/2;
vecRhsFt = zeros(Nelems*C_Nbases,1);


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
    
    elemUh = lastStep_DarcyUh((CurrElem-1)*D_Nbases+1:CurrElem*D_Nbases);
    %-- get the Diffusivity coefficients of elemP and elemM
    elem_DC = DarcymeshInfo.DiffusivityCoeffs(CurrElem);
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %-- the mat and vec setting
    rhs_f_temp = zeros(C_Nbases,1);
    
    %-- the integration on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D);

        % get the bases values on quad points
        [testPb_C, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_C);
            %> testPb, testPbx, testPby, [Npoints x NTbases_test].
        
        

        [~, trialPbx_D, trialPby_D] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_Darcy);
        value_vector_u_1 = c_darcy.*(trialPbx_D*elemUh);
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
        value_vector_u_2 = c_darcy.*(trialPby_D*elemUh);
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
            
        %--- the funcValue may be chosen by case.
%         valuefunc_f = pde.ct(t,phyGpoints(:,1),phyGpoints(:,2)) ...
%             + value_vector_u_1.*pde.cx(t,phyGpoints(:,1),phyGpoints(:,2)) ...
%             + value_vector_u_2.*pde.cy(t,phyGpoints(:,1),phyGpoints(:,2)) ...
%             - elem_DC*(pde.cxx(t,phyGpoints(:,1),phyGpoints(:,2))) ...
%             - elem_DC*(pde.cyy(t,phyGpoints(:,1),phyGpoints(:,2)));
        
        valuefunc_f = pde.fc_D(t,phyGpoints(:,1),phyGpoints(:,2));
        
        %-- VecF   
        rhs_f_temp = rhs_f_temp + testPb_C' * ( phyGweights.*valuefunc_f );
    end % for nt
    
    
    %-- value_rhs_f -- add to the position in the big vector
    Row_vec = ((CurrElem-1)*C_Nbases+1 : CurrElem*C_Nbases);
    vecRhsFt(Row_vec,1) = vecRhsFt(Row_vec,1) + rhs_f_temp;
    
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
    
    elemUh1 = lastStep_StokesUh((CurrElem-1)*S_Nbases+1:CurrElem*S_Nbases);
    elemUh2 = lastStep_StokesUh((S_Nelems+CurrElem-1)*S_Nbases+1:(S_Nelems+CurrElem)*S_Nbases);
    %-- get the Diffusivity coefficients of elemP and elemM
    elem_DC = StokesmeshInfo.DiffusivityCoeffs(CurrElem);
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %-- the mat and vec setting
    rhs_f_temp = zeros(C_Nbases,1);
    
    %-- the integration on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D);

        % get the bases values on quad points
        [testPb_C, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_C);
            %> testPb, testPbx, testPby, [Npoints x NTbases_test].

        [trialPb_S, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_Stokes);
        value_vector_u_1 = trialPb_S*elemUh1;
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
        value_vector_u_2 = trialPb_S*elemUh2;
            %> [Npoints x Nbases]*[Nbases x 1] --> [Npoints x 1]
            
        %--- the funcValue may be chosen by case.
%         valuefunc_f = pde.ct(t,phyGpoints(:,1),phyGpoints(:,2)) ...
%             + value_vector_u_1.*pde.cx(t,phyGpoints(:,1),phyGpoints(:,2)) ...
%             + value_vector_u_2.*pde.cy(t,phyGpoints(:,1),phyGpoints(:,2)) ...
%             - elem_DC*(pde.cxx(t,phyGpoints(:,1),phyGpoints(:,2))) ...
%             - elem_DC*(pde.cyy(t,phyGpoints(:,1),phyGpoints(:,2)));
        
        valuefunc_f = pde.fc_S(t,phyGpoints(:,1),phyGpoints(:,2));
        
        %-- VecF   
        rhs_f_temp = rhs_f_temp + testPb_C' * ( phyGweights.*valuefunc_f );
    end % for nt
    
    
    %-- value_rhs_f -- add to the position in the big vector
    Row_vec = ((D_Nelems+CurrElem-1)*C_Nbases+1 : (D_Nelems+CurrElem)*C_Nbases);
    vecRhsFt(Row_vec,1) = vecRhsFt(Row_vec,1) + rhs_f_temp;
    
    %<<-- End Part II --------------------------------------------------------------------------------------
    
end % for ii

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

function [matElemDuDv,matElemPVx,matElemPVy, matElempq,vecRhsF1,vecRhsF2,vecLambdaP] = ...
    matElemCoeffsDtrialDtestStokes(Coeffs_func,func_f1,func_f2, meshInfo, formulaGauss2D, degreek_u, degreek_q)
%
%   In this function we compute the elem integral.
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases_u denote the number of LOCAL velocity u bases on each K of Th.
%               NTbases_q denote the number of LOCAL pressure p bases on each K of Th.
%
%   input:
%       Coeffs, the cell-type, here 
%                   Coeffs{1} is the (\int_T Coeffs{1}uv dx),
%                   Coeffs{2}, and Coeffs{3} are the K=[k11, k12; k21, k22], 
%                   let k12=0, k21=0, the Coeffs{2}=k11, Coeffs{3}= k22. 
%                   Coeffs{2} is the (\int_T Coeffs{2} u_x v dx). 
%                   Coeffs{3} is the (\int_T Coeffs{3} u_y v dx).
%       rhs_f_func, the cell-type, here
%                   rhs_f_func{i} is the first rhs f func.
%       meshInfos, the mesh information.
%       formulaGauss2D, the 1d Gauss quadrature formula, size: a matrix, [Npoints x 3],
%               the first column is the x-coordinates of all Gauss-Points,
%               the second column is the y-coordinates of all Gauss-Points,
%               the third is the weights of all Gauss-Points.
%       g1Trial_k, the group1 trial polynomial degree k.
%       g1Test_k, the group1 test polynomial degree k.
%           %> here we consider that, in the Stokes equation, there will have VELOCITY u and PRESSURE p, 
%           %> so there have \int_e uv ds, then g1Trial_k=u_degree_k, g1Test_k=u_degree_k. 
%           %> And also have \int_e pv ds, then g2Trial_k=p_degree_k, g2Test_k=u_degree_k. 
%
%   output:
%       matElemDuDv, [Nelems*NTg1bases_test x NTelems*NTg1bases_trial], \sum_e \int_e {grad u}\cdot[v]. 
%       rhs_fh, [Nelems*NTg1bases_test x 1].
%
%
%   YcZhang 11/8/2017
%
%   Last modified 16/8/2017
%
%

Nelems = meshInfo.Nelems;

Coeff_k11 = Coeffs_func{1};
Coeff_k22 = Coeffs_func{2};
Coeff_one = @(x,y) 1+0.*x;


% group1 setting
NTbases_u = (degreek_u+1)*(degreek_u+2)/2;
NTbases_q = (degreek_q+1)*(degreek_q+2)/2;

matElemDuDv = sparse(Nelems*NTbases_u,Nelems*NTbases_u);
matElemPVx = sparse(Nelems*NTbases_u,Nelems*NTbases_q);
matElemPVy = sparse(Nelems*NTbases_u,Nelems*NTbases_q);
matElempq = sparse(Nelems*NTbases_q,Nelems*NTbases_q);
vecRhsF1 = zeros(Nelems*NTbases_u,1);
vecRhsF2 = zeros(Nelems*NTbases_u,1);
vecLambdaP = zeros(Nelems*NTbases_q,1);


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

for ii = 1:Nelems
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    singleElem = meshInfo.elem{ii,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges on singleElem.
    
    concavePointElem_ii = meshInfo.concavePointElem(ii); % [1 x 1], the concavePointElem of ii-the elem.
    coordv = meshInfo.node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    centroidElem = meshInfo.centroidElem(ii,:); % [1 x 2], the centroid(xing xin) (x-coord, y-coord) of element.
    
    coordTri0Elem = getcoordTri0Elem(singleNE, concavePointElem_ii, centroidElem, coordv); 
        %> if singleNE==3, coordTri0Elem, [3 x 2]. else coordTri0Elem, [3*singleNE x 2].
    
    elem_xT = centroidElem(1);
    elem_yT = centroidElem(2); 
    elem_hT = meshInfo.hElem(ii); 
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %-- the mat and vec setting
    DxTrialDxTest_temp = zeros(NTbases_u, NTbases_u);
    DyTrialDyTest_temp = zeros(NTbases_u, NTbases_u);
    DxuD0p_temp = zeros(NTbases_u, NTbases_q);
    DyuD0p_temp = zeros(NTbases_u, NTbases_q);
    pq_temp = zeros(NTbases_q, NTbases_q);
    rhs_f1_temp = zeros(NTbases_u,1);
    rhs_f2_temp = zeros(NTbases_u,1);
    lambdaP_temp = zeros(NTbases_q,1);
    
    %-- the integration on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D);

        % get the bases values on quad points
        [u_Pb, u_Pbx, u_Pby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_u);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases_trial].
        [q_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_q);
            %> testPb, testPbx, testPby, [Npoints x NTbases_test].
        
        %--- the funcValue may be chosen by case.
        valueCoeff_one = Coeff_one(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k11 = Coeff_k11(phyGpoints(:,1), phyGpoints(:,2));
        valueCoeff_k22 = Coeff_k22(phyGpoints(:,1), phyGpoints(:,2));
        valuefunc_f1 = func_f1(phyGpoints(:,1), phyGpoints(:,2));
        valuefunc_f2 = func_f2(phyGpoints(:,1), phyGpoints(:,2));
        
        
        %-- MatElem
        DxTrialDxTest_temp = DxTrialDxTest_temp + getMatOnElem(valueCoeff_k11, u_Pbx, u_Pbx, phyGweights);
        DyTrialDyTest_temp = DyTrialDyTest_temp + getMatOnElem(valueCoeff_k22, u_Pby, u_Pby, phyGweights);
        DxuD0p_temp = DxuD0p_temp + getMatOnElem(valueCoeff_one, q_Pb, u_Pbx, phyGweights);
        DyuD0p_temp = DyuD0p_temp + getMatOnElem(valueCoeff_one, q_Pb, u_Pby, phyGweights);
        pq_temp = pq_temp + getMatOnElem(valueCoeff_one, q_Pb, q_Pb, phyGweights);
        
        %-- VecF   
        rhs_f1_temp = rhs_f1_temp + u_Pb' * ( phyGweights.*valuefunc_f1 );
        rhs_f2_temp = rhs_f2_temp + u_Pb' * ( phyGweights.*valuefunc_f2 );
        
        %-- lambda term 
        lambdaP_temp = lambdaP_temp + q_Pb' * phyGweights;
    end % for nt
    
    %-- add MatElem to the position in the big mat
    uvRow_elem = ((ii-1)*NTbases_u+1 : ii*NTbases_u)' ...
        * ones(1,NTbases_u); % [NTbases_u x NTbases_u]
    uvCol_elem = ones(NTbases_u, 1) ...
        *((ii-1)*NTbases_u+1 : ii*NTbases_u); % [NTbases_u x NTbases_u]
    
    qvRow_elem = ((ii-1)*NTbases_u+1 : ii*NTbases_u)' ...
        * ones(1,NTbases_q); % [NTbases_u x NTbases_q]
    qvCol_elem = ones(NTbases_u, 1) ...
        *((ii-1)*NTbases_q+1 : ii*NTbases_q); % [NTbases_u x NTbases_q]
    
    pqRow_elem = ((ii-1)*NTbases_q+1 : ii*NTbases_q)' ...
        * ones(1,NTbases_q); % [NTbases_q x NTbases_q]
    pqCol_elem = ones(NTbases_q, 1) ...
        *((ii-1)*NTbases_q+1 : ii*NTbases_q); % [NTbases_q x NTbases_q]
    
    matElemDuDv = matElemDuDv ...
        + sparse(uvRow_elem(:), uvCol_elem(:), ...
        DxTrialDxTest_temp(:) + DyTrialDyTest_temp(:), ...
        Nelems*NTbases_u, Nelems*NTbases_u);
    matElemPVx = matElemPVx ...
        + sparse(qvRow_elem(:), qvCol_elem(:), ...
        DxuD0p_temp(:), ...
        Nelems*NTbases_u, Nelems*NTbases_q);
    matElemPVy = matElemPVy ...
        + sparse(qvRow_elem(:), qvCol_elem(:), ...
        DyuD0p_temp(:), ...
        Nelems*NTbases_u, Nelems*NTbases_q);
    matElempq = matElempq ...
        + sparse(pqRow_elem(:), pqCol_elem(:), ...
        pq_temp(:), ...
        Nelems*NTbases_q, Nelems*NTbases_q);
    
    %-- value_rhs_f -- add to the position in the big vector
    fRow_vec = ((ii-1)*NTbases_u+1 : ii*NTbases_u);
    vecRhsF1(fRow_vec,1) = vecRhsF1(fRow_vec,1) + rhs_f1_temp;
    vecRhsF2(fRow_vec,1) = vecRhsF2(fRow_vec,1) + rhs_f2_temp;
    
    %-- lambda term
    lamRow_vec = ((ii-1)*NTbases_q+1 : ii*NTbases_q);
    vecLambdaP(lamRow_vec,1) = vecLambdaP(lamRow_vec,1) + lambdaP_temp;
    
    %<<-- End Part II --------------------------------------------------------------------------------------
    
end % for ii

end % function matElemCoeffsTrialTest



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

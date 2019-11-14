function [g1_matEdgeAgradtrialJtest, g1_matEdgeAgradtestJtrial, g1_matEdgeJtrialJtest, varargout] = ...
    matInterEdgeAverJumpPoisson_backup(Coeffs_func, meshInfo, option,formulaGauss1D,g1Trial_k,g1Test_k,varargin)
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
%       g1Trial_k, the group1 trial polynomial degree k.
%       g1Test_k, the group1 test polynomial degree k.
%           %> here we consider that, in the Stokes equation, there will have VELOCITY u and PRESSURE p, 
%           %> so there have \int_e uv ds, then g1Trial_k=u_degree_k, g1Test_k=u_degree_k. 
%           %> And also have \int_e pv ds, then g2Trial_k=p_degree_k, g2Test_k=u_degree_k. 
%       varargin, if length(varargin)~=0, the varargin must include g2Trial_k, g2Test_k, (i.e. 
%                     group2 trial polynomian degree k, and group2 test polynominal degree k).
%
%   output:
%       g1_matEdgeAgradtrialJtest, [Nelems*NTg1bases_test x NTelems*NTg1bases_trial], \sum_e \int_e {grad u}\cdot[v]. 
%       g1_matEdgeAgradtestJtrial, [Nelems*NTg1bases_test x NTelems*NTg1bases_trial], \sum_e \int_e {grad v}\cdot[u]. 
%       g1_matEdgeJtrialJtest, [Nelems*NTg1bases_test x NTelems*NTg1bases_trial], \sum_e \int_e [u]\cdot[v]. 
%
%
%   YcZhang 5/8/2017
%
%   Last modified 10/8/2017
%
%

% % integrationAll = 0;
% % integrationAgradtrialJtest = 0; % \int_e {grad_trial_function}\cdot[test_function] ds
% % integrationAgradtestJtrial = 0; % \int_e {grad_test_function}\cdot[trial_function] ds
% % integrationJtrialJtest = 0; % \int_e [trial]\cdot[test] ds
% % 
% % for n = 1:length(varargin)
% %     if strcmp(varargin{n},'integrationAll')
% %         integrationAll = 1;
% %     end % if
% %     
% %     if strcmp(varargin{n},'integrationAgradtrialJtest') || integrationAll
% %         integrationAgradtrialJtest = 1;
% %     end
% %     
% %     if strcmp(varargin{n},'integrationAgradtestJtrial') || integrationAll
% %         integrationAgradtestJtrial = 1;
% %     end
% %     
% %     if strcmp(varargin{n},'integrationJtrialJtest') || integrationAll
% %         integrationJtrialJtest = 1;
% %     end
% % end % for n

% mesh information: interior edges 
Nelems = meshInfo.Nelems;
interEdges = meshInfo.interEdgeIndex; % here we need the interEdges information.

p_epsilon = option.p_epsilon;
p_sigma = option.p_sigma;
p_beta = option.p_beta;

Coeff_k11 = Coeffs_func{1};
Coeff_k22 = Coeffs_func{2};

% group1 setting
NTg1bases_trial = (g1Trial_k+1)*(g1Trial_k+2)/2;
NTg1bases_test = (g1Test_k+1)*(g1Test_k+2)/2;

g1_matEdgeAgradtrialJtest = sparse(Nelems*NTg1bases_test,Nelems*NTg1bases_trial);
g1_matEdgeAgradtestJtrial = sparse(Nelems*NTg1bases_test,Nelems*NTg1bases_trial);
g1_matEdgeJtrialJtest = sparse(Nelems*NTg1bases_test,Nelems*NTg1bases_trial);

% group2 setting
have_group2 = 0;
if length(varargin)==2
    have_group2 = 1;
    g2Trial_k = varargin{1};
    g2Test_k = varargin{2};   
    
    NTg2bases_trial = (g2Trial_k+1)*(g2Trial_k+2)/2;
    NTg2bases_test = (g2Test_k+1)*(g2Test_k+2)/2;
    
    g2_matEdgeAgradtrialJtest = saprse(Nelems*NTg2bases_test,Nelems*NTg2bases_trial);
    g2_matEdgeAgradtestJtrial = saprse(Nelems*NTg2bases_test,Nelems*NTg2bases_trial);
    g2_matEdgeJtrialJtest = saprse(Nelems*NTg2bases_test,Nelems*NTg2bases_trial);
end

for ii = 1:length(interEdges)
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th edge.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    interEindx = interEdges(ii); 
    ePoint1 = meshInfo.node(meshInfo.edge(interEindx,1),:); % [1 x 2], the (x-coord, y-coord) of the first point of ii-th edge.
    ePoint2 = meshInfo.node(meshInfo.edge(interEindx,2),:); % [1 x 2], the (x-coord, y-coord) of the second point of ii-th edge.
    areaEdge = meshInfo.areaEdge(interEindx);
    
    phyGpointsX = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(1),ePoint2(1));
        %> [Npoints x 1], the physical x-coordinates Gauss points on ii-th edge.
    phyGpointsY = meshInfo.mapRefE2PhyE(formulaGauss1D(:,1),ePoint1(2),ePoint2(2));
        %> [Npoints x 1], the physical y-coordinates Gauss points on ii-th edge.
    phyGweights = meshInfo.areaEdge(interEindx) * formulaGauss1D(:,2); 
        %> [Npoints x 1], the physical Gauss weights on ii-th edge.
    
%     [ele1, ele2, local_e1, local_e2] = meshInfo.edge2elem(interEindx,:); 
    ele1 = meshInfo.edge2elem(interEindx,1);
    ele2 = meshInfo.edge2elem(interEindx,2);
    local_e1 = meshInfo.edge2elem(interEindx,3);
    local_e2 = meshInfo.edge2elem(interEindx,4);
        %> [n1, n2]=edge2elem(k,1:2) stands for the n1-th elem and n2-th elem share the k-th edge.
        %> local_e1=edge2elem(k,3) stands for the local edge index in n1-th elem of the k-th edge .
        %> also, local_e2=edge2elem(k,4) stands for the local edge index in n2-th elem of the k-th edge .
    eleP = min(ele1, ele2); % let the smaller number of element as the element +.
    eleM = max(ele1, ele2); % let the larger number of element as the element -.
    
    if eleP == ele1
        eIndx_eleP = local_e1;
        eIndx_eleM = local_e2;
    else
        eIndx_eleP = local_e2;
        eIndx_eleM = local_e1;
    end 
    
    eleP_nu = meshInfo.nuEdge0Elem{eleP}(:,eIndx_eleP); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleP.
    eleM_nu = meshInfo.nuEdge0Elem{eleM}(:,eIndx_eleM); 
        %> [1 x 2], the outward unit normal vector of ii-th edge on eleM.
    
    eleP_xT = meshInfo.baryElem(eleP,1);  
    eleP_yT = meshInfo.baryElem(eleP,2);
    eleP_hT = meshInfo.diameters(eleP);
        %> the barycenter and diameter of eleP.
    eleM_xT = meshInfo.baryElem(eleM,1);  
    eleM_yT = meshInfo.baryElem(eleM,2);
    eleM_hT = meshInfo.diameters(eleM);
        %> the barycenter and diameter of eleM.
        
    %-- get the value of Coeffs func at Gauss Points
    valueCoeff_k11 = Coeff_k11(phyGpointsX, phyGpointsY); % [Npoints x 1]
    valueCoeff_k22 = Coeff_k22(phyGpointsX, phyGpointsY); % [Npoints x 1]
    
    %-- subPart -- get the GROUP1 bases on quad
    [eleP_g1trialPb, eleP_g1trialPbx, eleP_g1trialPby] = ...
        localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, g1Trial_k);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
    [eleM_g1trialPb, eleM_g1trialPbx, eleM_g1trialPby] = ...
        localBases2D(eleM_xT, eleM_yT, eleM_hT, phyGpointsX, phyGpointsY, g1Trial_k);
        %> eleMtrialPb, eleMtrialPbx, eleMtrialPby, [Npoints x NTg1bases_trial]
    
    [eleP_g1testPb, eleP_g1testPbx, eleP_g1testPby] = ...
        localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, g1Test_k);
        %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg1bases_test]
    [eleM_g1testPb, eleM_g1testPbx, eleM_g1testPby] = ...
        localBases2D(eleM_xT, eleM_yT, eleM_hT, phyGpointsX, phyGpointsY, g1Test_k);
        %> eleMtestPb, eleMtestPbx, eleMtestPby, [Npoints x NTg1bases_test]
    
    %-- sbuPart -- get the GROUP2 bases on quad
    if have_group2
        [eleP_g2trialPb, eleP_g2trialPbx, eleP_g2trialPby] = ...
            localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, g2Trial_k);
            %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg2bases_trial]
        [eleM_g2trialPb, eleM_g2trialPbx, eleM_g2trialPby] = ...
            localBases2D(eleM_xT, eleM_yT, eleM_hT, phyGpointsX, phyGpointsY, g2Trial_k);
            %> eleMtrialPb, eleMtrialPbx, eleMtrialPby, [Npoints x NTg2bases_trial]
    
        [eleP_g2testPb, eleP_g2testPbx, eleP_g2testPby] = ...
            localBases2D(eleP_xT, eleP_yT, eleP_hT, phyGpointsX, phyGpointsY, g2Test_k);
            %> elePtestPb, elePtestPbx, elePtestPby, [Npoints x NTg2bases_test]
        [eleM_g2testPb, eleM_g2testPbx, eleM_g2testPby] = ...
            localBases2D(eleM_xT, eleM_yT, eleM_hT, phyGpointsX, phyGpointsY, g2Test_k);
            %> eleMtestPb, eleMtestPbx, eleMtestPby, [Npoints x NTg2bases_test]
    end % if have_group2_k
        
    %<<-- End Part I ---------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%
    
    %% Part II, compute the (\int_e {.}\cdot[.] ds, \int_e [.]\cdot[.] ds)
    %> {u} = 0.5*(u^plus + u^minus), is a scaler.
    %> [u] = u^plus \nu^plus + u^minus \nu^minus, is a vector.
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    %
    %-- subPart I -- the GROUP1 
    % get the position in the big mat
    g1_pp_Row = ((eleP-1)*NTg1bases_test+1 : eleP*NTg1bases_test)' ...
        * ones(1,NTg1bases_trial); % [NTg1bases_test x NTg1bases_trial]
    g1_pp_Col = ones(NTg1bases_test, 1) ...
        *((eleP-1)*NTg1bases_trial+1 : eleP*NTg1bases_trial); 
        %> [NTg1bases_test x NTg1bases_trial]
    
    g1_pm_Row = g1_pp_Row; % [NTg1bases_test x NTg1bases_trial]
    g1_pm_Col = ones(NTg1bases_test, 1) ...
        *((eleM-1)*NTg1bases_trial+1 : eleM*NTg1bases_trial); 
        %> [NTg1bases_test x NTg1bases_trial]
    
    g1_mp_Row = ((eleM-1)*NTg1bases_test+1 : eleM*NTg1bases_test)' ...
        * ones(1,NTg1bases_trial); % [NTg1bases_test x NTg1bases_trial]
    g1_mp_Col = g1_pp_Col; % [NTg1bases_test x NTg1bases_trial]
    
    g1_mm_Row = g1_mp_Row; % [NTg1bases_test x NTg1bases_trial]
    g1_mm_Col = g1_pm_Col;
    
    g1_Row = [g1_pp_Row, g1_pm_Row;
        g1_mp_Row, g1_mm_Row];
    g1_Col = [g1_pp_Col, g1_pm_Col;
        g1_mp_Col, g1_mm_Col];
    
    % get the matrix
    [g1_A1, g1_A2, g1_A3] = getMatOnEdge(phyGweights, ...
        valueCoeff_k11, valueCoeff_k22, ...
        eleP_nu, ...
        eleP_g1trialPb, eleP_g1trialPbx, eleP_g1trialPby, ...
        eleP_g1testPb, eleP_g1testPbx, eleP_g1testPby, ...
        eleM_nu, ...
        eleM_g1trialPb, eleM_g1trialPbx, eleM_g1trialPby, ...
        eleM_g1testPb, eleM_g1testPbx, eleM_g1testPby);
        %> g1_A1, g1_A2, g1_A3, [2*NTg1bases_test  x 2*NTg1bases_trial]
        %> g1_A1=[A1_pp, A1_pm; 
        %>            A1_mp, A1_mm];
        %> the same to g1_A2, g1_A3.
    
    g1_matEdgeAgradtrialJtest = g1_matEdgeAgradtrialJtest ...
        + sparse(g1_Row(:), g1_Col(:), -g1_A1(:), Nelems*NTg1bases_test, Nelems*NTg1bases_trial);
    g1_matEdgeAgradtestJtrial = g1_matEdgeAgradtestJtrial ...
        + sparse(g1_Row(:), g1_Col(:), p_epsilon*g1_A2(:), Nelems*NTg1bases_test, Nelems*NTg1bases_trial);
    g1_matEdgeJtrialJtest = g1_matEdgeJtrialJtest ...
        + sparse(g1_Row(:), g1_Col(:), (p_sigma/areaEdge^p_beta)*g1_A3(:), Nelems*NTg1bases_test, Nelems*NTg1bases_trial);
    
    
    %-- subPart II -- the GROUP2
    if have_group2
        % get the position in the big mat
        g2_pp_Row = ((eleP-1)*NTg2bases_test+1 : eleP*NTg2bases_test)' ...
            * ones(1,NTg2bases_trial); % [NTg2bases_test x NTg2bases_trial]
        g2_pp_Col = ones(NTg2bases_test, 1) ...
            *((eleP-1)*NTg2bases_trial+1 : eleP*NTg2bases_trial); 
            %> [NTg2bases_test x NTg2bases_trial]

        g2_pm_Row = g2_pp_Row; % [NTg2bases_test x NTg2bases_trial]
        g2_pm_Col = ones(NTg2bases_test, 1) ...
            *((eleM-1)*NTg2bases_trial+1 : eleM*NTg2bases_trial); 
            %> [NTg2bases_test x NTg2bases_trial]

        g2_mp_Row = ((eleM-1)*NTg2bases_test+1 : eleM*NTg2bases_test)' ...
            * ones(1,NTg2bases_trial); % [NTg2bases_test x NTg2bases_trial]
        g2_mp_Col = g2_pp_Col; % [NTg2bases_test x NTg2bases_trial]

        g2_mm_Row = g2_mp_Row; % [NTg2bases_test x NTg2bases_trial]
        g2_mm_Col = g2_pm_Col;

        g2_Row = [g2_pp_Row, g2_pm_Row;
            g2_mp_Row, g2_mm_Row];
        g2_Col = [g2_pp_Col, g2_pm_Col;
            g2_mp_Col, g2_mm_Col];

        % get the matrix
        [g2_A1, g2_A2, g2_A3] = getMatOnEdge(phyGweights, ...
            valueCoeff_k11, valueCoeff_k22, ...
            eleP_nu, ...
            eleP_g2trialPb, eleP_g2trialPbx, eleP_g2trialPby, ...
            eleP_g2testPb, eleP_g2testPbx, eleP_g2testPby, ...
            eleM_nu, ...
            eleM_g2trialPb, eleM_g2trialPbx, eleM_g2trialPby, ...
            eleM_g2testPb, eleM_g2testPbx, eleM_g2testPby);
            %> g2_A1, g2_A2, g2_A3, [2*NTg2bases_test  x 2*NTg2bases_trial]
            %> g2_A1=[A1_pp, A1_pm; 
            %>            A1_mp, A1_mm];
            %> the same to g2_A2, g2_A3.

        g2_matEdgeAgradtrialJtest = g2_matEdgeAgradtrialJtest ...
            + sparse(g2_Row(:), g2_Col(:), -g2_A1(:), Nelems*NTg2bases_test, Nelems*NTg2bases_trial);
        g2_matEdgeAgradtestJtrial = g2_matEdgeAgradtestJtrial ...
            + sparse(g2_Row(:), g2_Col(:), p_epsilon*g2_A2(:), Nelems*NTg2bases_test, Nelems*NTg2bases_trial);
        g2_matEdgeJtrialJtest = g2_matEdgeJtrialJtest ...
            + sparse(g2_Row(:), g2_Col(:), (p_sigma/areaEdge^p_beta)*g2_A3(:), Nelems*NTg2bases_test, Nelems*NTg2bases_trial);
    end % if have_group2
        
    %<<-- End Part II --------------------------------------------------------------------------------------
end % for ii

if have_group2
    varargout{1} = g2_matEdgeAgradtrialJtest;
    varargout{2} = g2_matEdgeAgradtestJtrial;
    varargout{3} = g2_matEdgeJtrialJtest;
end % have_group2

end % function matInterEdgeAverJump







%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%% the sub function 
function [A1, A2, A3] = getMatOnEdge(phyGaussweights, ...
    valueCoeff_1, valueCoeff_2, ...
    eleP_nu, ...
    eleP_trialPb, eleP_trialPbx, eleP_trialPby, ...
    eleP_testPb, eleP_testPbx, eleP_testPby, ...
    eleM_nu, ...
    eleM_trialPb, eleM_trialPbx, eleM_trialPby, ...
    eleM_testPb, eleM_testPbx, eleM_testPby)

    %-- subPart I, \int_e {grad_trial}\cdot[test]
    % {grad_trial}\cdot[test] = 1/2(grad_trial^p + grad_trial^m)\cdot(test^p nu^p + test^m nu^m)
    %AgradtestJtrial_pp = zeros(NTbases_test, NTbases_trial);    
    A1_pp = 0.5*( eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A1_pm = 0.5*( eleP_nu(1) * eleP_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPbx) ...
        + eleP_nu(2) * eleP_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_mp = 0.5*( eleM_nu(1) * eleM_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPbx) ...
        + eleM_nu(2) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A1_mm = 0.5*( eleM_nu(1) * eleM_testPb' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPbx) ...
        + eleM_nu(2) * eleM_testPb' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPby) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A1 = [A1_pp, A1_pm;
        A1_mp, A1_mm];
        
    
    %-- subPart II, \int_e {grad_test}\cdot[trial]
    A2_pp = 0.5*( eleP_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + eleP_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A2_pm = 0.5*( eleM_nu(1) * eleP_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPb) ...
        + eleM_nu(2) * eleP_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_mp = 0.5*( eleP_nu(1) * eleM_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleP_trialPb) ...
        + eleP_nu(2) * eleM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleP_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A2_mm = 0.5*( eleM_nu(1) * eleM_testPbx' * bsxfun(@times,valueCoeff_1.*phyGaussweights,eleM_trialPb) ...
        + eleM_nu(2) * eleM_testPby' * bsxfun(@times,valueCoeff_2.*phyGaussweights,eleM_trialPb) );
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
       
    A2 = [A2_pp, A2_pm;
        A2_mp, A2_mm];    
    
        
    %-- subPart III, \int_e [trial]\cdot[test]    
    A3_pp = eleP_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
    
    A3_pm = -eleP_testPb' * bsxfun(@times,phyGaussweights,eleM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3_mp = -eleM_testPb' * bsxfun(@times,phyGaussweights,eleP_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3_mm = eleM_testPb' * bsxfun(@times,phyGaussweights,eleM_trialPb);
        %> [NTbases_test x Npoints] * [Npoints x NTbases_trial] = [NTbases_test x NTbases_trial].
        
    A3 = [A3_pp, A3_pm;
        A3_mp, A3_mm];
end % function getMatOnT

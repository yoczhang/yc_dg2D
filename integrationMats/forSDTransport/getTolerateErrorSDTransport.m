function L2err = getTolerateErrorSDTransport(...
        new_DarcyUh, new_StokesUh, new_Ch, lastStep_DarcyUh, lastStep_StokesUh, lastStep_Ch, ...
        DarcymeshInfo, StokesmeshInfo, formulaGauss2D, ...
        degreek_Darcy, degreek_Stokes, degreek_C)
%
%   This funciton is to get the tolerate error, in the iteration.
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases_u denote the number of LOCAL velocity u bases on each K of Th.
%
%   Input:
%       Uhnew, the new Uh.
%       Uhold, the old Uh.
%       meshInfo, the mesh information.
%       formulaGauss2D, the 2d Gauss quadrature formula, size: a matrix, [Npoints x 3],
%               the first column is the x-coordinates of all Gauss-Points,
%               the second column is the y-coordinates of all Gauss-Points,
%               the third is the weights of all Gauss-Points.
%       degreek_u, the polynomial degree k of Uh.
%
%
%   YcZhang 21/10/2017
%   Last modified 21/10/2017
%
%

S_Nelems = StokesmeshInfo.Nelems;
D_Nelems = DarcymeshInfo.Nelems;
S_Nbases = (degreek_Stokes+1)*(degreek_Stokes+2)/2;
D_Nbases = (degreek_Darcy+1)*(degreek_Darcy+2)/2;
C_Nbases = (degreek_C+1)*(degreek_C+2)/2;

% basic setting
dof_S_u1 = S_Nelems*S_Nbases;
% dof_S_p = S_Nelems*((degreek_Stokes-1+1)*(degreek_Stokes-1+2)/2);
% dof_D_u = D_Nelems*D_Nbases;

% dof_SD = (dof_D_u+2*dof_S_u1+dof_S_p);
% Nelems = S_Nelems + D_Nelems;

L2err = 0;
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
    
    D_UhnewElem = new_DarcyUh((CurrElem-1)*D_Nbases+1:CurrElem*D_Nbases);
    D_UholdElem = lastStep_DarcyUh((CurrElem-1)*D_Nbases+1:CurrElem*D_Nbases);
    C_UhnewElem = new_Ch((CurrElem-1)*C_Nbases+1:CurrElem*C_Nbases);
    C_UholdElem = lastStep_Ch((CurrElem-1)*C_Nbases+1:CurrElem*C_Nbases);
    
    %<<-- End Part I ---------------------------------------------------------------------------------------

%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II, compute the integration
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    
    %-- the integration on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D);

        % get the bases values on quad points
        [D_u_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_Darcy);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases_trial].
        [C_u_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_C);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases_trial].
            
        %--- the funcValue may be chosen by case.
        L2err = L2err + phyGweights' * (D_u_Pb*D_UhnewElem - D_u_Pb*D_UholdElem).^2 ...
            + phyGweights' * (C_u_Pb*C_UhnewElem - C_u_Pb*C_UholdElem).^2;
    end % for nt
    
    %<<-- End Part II ---------------------------------------------------------------------------------------

end % for CurrElem

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
    
    S_Uh1newElem = new_StokesUh((CurrElem-1)*S_Nbases+1:CurrElem*S_Nbases);
    S_Uh1oldElem = lastStep_StokesUh((CurrElem-1)*S_Nbases+1:CurrElem*S_Nbases);
    S_Uh2newElem = new_StokesUh(dof_S_u1+(CurrElem-1)*S_Nbases+1:dof_S_u1+CurrElem*S_Nbases);
    S_Uh2oldElem = lastStep_StokesUh(dof_S_u1+(CurrElem-1)*S_Nbases+1:dof_S_u1+CurrElem*S_Nbases);
    
    C_UhnewElem = new_Ch((D_Nelems+CurrElem-1)*C_Nbases+1:(D_Nelems+CurrElem)*C_Nbases);
    C_UholdElem = lastStep_Ch((D_Nelems+CurrElem-1)*C_Nbases+1:(D_Nelems+CurrElem)*C_Nbases);
    
    %<<-- End Part I ---------------------------------------------------------------------------------------

%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II, compute the integration
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    
    %-- the integration on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D);

        % get the bases values on quad points
        [S_u_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_Darcy);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases_trial].
        [C_u_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek_C);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases_trial].
            
        %--- the funcValue may be chosen by case.
        L2err = L2err + phyGweights' * (S_u_Pb*S_Uh1newElem - S_u_Pb*S_Uh1oldElem).^2 ...
            + phyGweights' * (S_u_Pb*S_Uh2newElem - S_u_Pb*S_Uh2oldElem).^2 ...
            + phyGweights' * (C_u_Pb*C_UhnewElem - C_u_Pb*C_UholdElem).^2;
    end % for nt
    
    %<<-- End Part II ---------------------------------------------------------------------------------------

end % for CurrElem

L2err = sqrt(L2err);

end % function
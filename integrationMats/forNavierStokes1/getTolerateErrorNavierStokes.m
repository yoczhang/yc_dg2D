function L2err = getTolerateErrorNavierStokes(Uhnew, Uhold, meshInfo, formulaGauss2D, degreek)
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

Nelems = meshInfo.Nelems;

% basic setting
NTbases = (degreek+1)*(degreek+2)/2;
dof_u = Nelems*NTbases;

L2err = 0;
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
    
    Uh1newElem = Uhnew((CurrElem-1)*NTbases+1:CurrElem*NTbases);
    Uh2newElem = Uhnew(dof_u+(CurrElem-1)*NTbases+1:dof_u+CurrElem*NTbases);
    Uh1oldElem = Uhold((CurrElem-1)*NTbases+1:CurrElem*NTbases);
    Uh2oldElem = Uhold(dof_u+(CurrElem-1)*NTbases+1:dof_u+CurrElem*NTbases);
    
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
        [u_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), degreek);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases_trial].
        
        %--- the funcValue may be chosen by case.
        L2err = L2err + phyGweights' * (u_Pb*Uh1newElem - u_Pb*Uh1oldElem).^2 ...
            + phyGweights' * (u_Pb*Uh2newElem - u_Pb*Uh2oldElem).^2;
    end % for nt
    
    %<<-- End Part II ---------------------------------------------------------------------------------------

end % for CurrElem

L2err = sqrt(L2err);

end % function
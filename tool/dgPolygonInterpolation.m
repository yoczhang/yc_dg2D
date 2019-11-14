function [DGM, DGT, Interp_uh, Interp_uhx, Interp_uhy] = dgPolygonInterpolation(exact_u, exact_ux, exact_uy, meshInfo)
%
%   NOTE that as for my present knowledge, I only can plot the P1 element
%   figure. So we use the trick that P2 elem can be treated as the 'refined-mesh'
%   P1 elem.
%
%   We let  Nvs is the number of veritces of the mesh Th.
%               Nelems denote the number of the elements of Th,
%               Nbases denote the number of LOCAL bases on each element.
%                       
%   input:
%       node, [Nvs x 2], the coordinates of vertices of the mesh Th.
%       elem, [Nelems x 3], golbal vertex indices of triangles accounting
%           for the counter-clockwise ordering.
%       bdFlag, only use to refine the mesh.
%       Uh, [Nbases*Nelems x 1], the HDG solution on Elems.
%       degreek, the polynomial degree.
%
%
%   YcZhang 3/8/2018
%
%   Last modified 4/8/2018
%

Nelems = meshInfo.Nelems;

% Lx = [0,1,0]; Ly=[0,0,1]; 
% refTri = [Lx', Ly', zeros(3,1)];
triCount = [1, 2, 3];

DGM = [];
DGT = [];
Interp_uh = [];
Interp_uhx = [];
Interp_uhy = [];
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
        
%     elem_xT = centroidElem(1);
%     elem_yT = centroidElem(2); 
%     elem_hT = meshInfo.hElem(ii);

    %<<-- End Part I ---------------------------------------------------------------------------------------
    
%% ------------------------------------------ Dividing line -------------------------------------------------- %%
%-------------------------------------------------------------------------------------------------------------------%

    %% Part II
    %
    %>>-- Begin Part II ------------------- THIS PART CAN BE MODIFIED --------------------------
    
    %-- on little triangles
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        
        % get the value of exact solution u, ux, uy at vertices of the triangle
        
        value_uh_onTri = exact_u(coordTri_nt(:,1),coordTri_nt(:,2));
        value_uhx_onTri = exact_ux(coordTri_nt(:,1),coordTri_nt(:,2));
        value_uhy_onTri = exact_uy(coordTri_nt(:,1),coordTri_nt(:,2));
        
        DGM_temp = [DGM; coordTri_nt];
        DGT_temp = [DGT; triCount];
        DGuh_temp = [Interp_uh; value_uh_onTri];
        DGuhx_temp = [Interp_uhx; value_uhx_onTri];
        DGuhy_temp = [Interp_uhy; value_uhy_onTri];
        
        DGM = DGM_temp;
        DGT = DGT_temp;
        Interp_uh = DGuh_temp;
        Interp_uhx = DGuhx_temp;
        Interp_uhy = DGuhy_temp;
        
        triCount = triCount + 3;
    end % for nt
   
    %<<-- End Part II --------------------------------------------------------------------------------------
    
end % for ii

%figure
%trisurf(DGT,DGM(:,1),DGM(:,2),DGuh);
%shading interp, xlabel('x'), ylabel('y'), colorbar , axis equal, axis off%, view(2)


% % if degreek == 1
% %     [LagrX, LagrY, dataLagr] = hdgProjDataDisc2DataLagr(Uh, elem, degreek);
% %         %> LagrX, [NLpoints x Nelems], map the Lx on REF elem TO PHY elem.
% %         %> LagrY, [NLpoints x Nelems], map the Ly on REF elem TO PHY elem.
% %         %> dataLagr, [NLpoints x Nelems], converts the representation matrix in the HDG/modal 
% %             %> basis to the respective representation matrix in a Lagrange/nodal basis.
% %     plotNbases = size(dataLagr,1);
% %     tri = bsxfun(@plus, ((1:g_mesh.numT)'-1)*plotNbases, 1:plotNbases);
% %     figure(1);
% %     trisurf(tri,LagrX(:),LagrY(:),dataLagr(:));
% %     shading interp, xlabel('x'), ylabel('y'), colorbar , view(2)
% % else
% %     [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
% %     
% % end % if

end % function
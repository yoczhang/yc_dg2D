function getStokesDarcyMesh()
clearvars;
close all
clc

addpath(genpath([ROOTDIR,'modules/upr/'])); 
addpath(genpath([ROOTDIR,'modules/book/'])); 

%%
% %------ mesh1 
% %-------------------- non-convex mesh ---------------------
% for nn = 1:4
%     close all;
%     meshN = num2str(2^(nn+1));
%     %--- Stokes Domain [0,1]x[1,2]
%     SmeshName = ['Smesh_nonconvex_[0,1]x[1,2]_',meshN];
%     [nodeOld,elemOld] = squarequadmesh([0,1,1,2],1/2^(nn+1)); 
%     [Snode,Selem] = non_convex_octagona_mesh(nodeOld,elemOld);
%     node = Snode; elem = Selem;
%     save(SmeshName, 'node', 'elem')
%  
%     %--- Darcy Domain [0,1]x[0,1]
%     DmeshName = ['Dmesh_nonconvex_[0,1]x[0,1]_',meshN];
%     [nodeOld,elemOld] = squarequadmesh([0,1,0,1],1/2^(nn+1)); 
%     [Dnode,Delem] = non_convex_octagona_mesh(nodeOld,elemOld);
%     node = Dnode; elem = Delem;
%     save(DmeshName, 'node', 'elem')
%    
%     PolyMshr_PlotMsh_twoDomain1(Snode,Selem,length(Selem),Dnode,Delem,length(Delem))
%     figure
%     PolyMshr_PlotMsh_twoDomain2(Snode,Selem,length(Selem),Dnode,Delem,length(Delem))
% end % nn


%%
% %------ mesh2
% %-------------------- Prrturbing Quadrilateral Mesh ---------------------
% for nn = 1:4
%     close all;
%     meshN = num2str(2^(nn+1));
%     %--- Stokes Domain [0,1]x[1,2]
%     SmeshName = ['Smesh_prrturbingQuad_[0,1]x[1,2]_',meshN];
%     [nodeOld,elemOld] = squarequadmesh([0,1,1,2],1/2^(nn+1)); 
%     meshInfo1 = polyMeshAuxStructure(nodeOld, elemOld);
%     [Snode, Selem] = generatePrrturbingQuadrilateralMesh(nodeOld, elemOld, meshInfo1.areaElem, 0.26);
%     node = Snode; elem = Selem;
%     col_elem = elem(:,4); elem(:,4) = elem(:,2); elem(:,2) = col_elem;
%     SmeshInfo = polyMeshAuxStructure(node, elem);
%     save(SmeshName, 'node', 'elem')
%  
%     %--- Darcy Domain [0,1]x[0,1]
%     DmeshName = ['Dmesh_prrturbingQuad_[0,1]x[0,1]_',meshN];
%     [nodeOld,elemOld] = squarequadmesh([0,1,0,1],1/2^(nn+1)); 
%     meshInfo1 = polyMeshAuxStructure(nodeOld, elemOld);
%     [Dnode, Delem] = generatePrrturbingQuadrilateralMesh(nodeOld, elemOld, meshInfo1.areaElem, 0.26);
%     node = Dnode; elem = Delem;
%     col_elem = elem(:,4); elem(:,4) = elem(:,2); elem(:,2) = col_elem;
%     DmeshInfo = polyMeshAuxStructure(node, elem);
%     save(DmeshName, 'node', 'elem')
%    
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
%     figure
%     PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
% end % nn


%%
% %------ mesh3
% %-------------------- Dual-Tri Mesh ---------------------
% for nn = 1:4
%     close all;
%     meshN = num2str(2^(nn+1));
%     %--- Stokes Domain [0,1]x[1,2]
%     SmeshName = ['Smesh_dualTri_[0,1]x[1,2]_',meshN];
%     %[nodeOld,elemOld] = squarequadmesh([0,1,1,2],1/2^(nn+1)); 
%     [nodeOld, elemOld] = generate_Tri_P_T(0,1,1,2,[1/2,1/2].^(nn+1));
%     meshInfo1 = polyMeshAuxStructure(nodeOld, elemOld);
%     [Snode, Selem] = dualmesh(nodeOld,elemOld);
% %     %--- contorted the mesh
% %     xx = Snode(:,1);
% %     yy = Snode(:,2);
% %     Snode(:,2) = yy.*((-4*abs(xx-0.5)+2).*(yy-2).^2.*(yy-1).^2+1);
%     
%     SmeshInfo = polyMeshAuxStructure(Snode, Selem);
%     node = SmeshInfo.node; elem = SmeshInfo.elem;
%     save(SmeshName, 'node', 'elem')
%  
%     %--- Darcy Domain [0,1]x[0,1]
%     DmeshName = ['Dmesh_dualTri_[0,1]x[0,1]_',meshN];
%     %[nodeOld,elemOld] = squarequadmesh([0,1,0,1],1/2^(nn+1)); 
%     [nodeOld, elemOld] = generate_Tri_P_T(0,1,0,1,[1/2,1/2].^(nn+1));
%     meshInfo1 = polyMeshAuxStructure(nodeOld, elemOld);
%     [Dnode, Delem] = dualmesh(nodeOld,elemOld);
%     DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
%     node = DmeshInfo.node; elem = DmeshInfo.elem;
%     save(DmeshName, 'node', 'elem')
%    
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
%     figure
%     PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
% end % nn


%%
% %------ mesh4
% %-------------------- Contorted Dual-Tri Mesh ---------------------
% for nn = 1:4
%     close all;
%     meshN = num2str(2^(nn+1));
%     %--- Stokes Domain [0,1]x[1,2]
%     SmeshName = ['Smesh_contortedDualTri_[0,1]x[1,2]_',meshN];
%     %[nodeOld,elemOld] = squarequadmesh([0,1,1,2],1/2^(nn+1)); 
%     [nodeOld, elemOld] = generate_Tri_P_T(0,1,1,2,[1/2,1/2].^(nn+1));
%     [Snode, Selem] = dualmesh(nodeOld,elemOld);
%     
%     %--- contorted the mesh
%     c_theta = 0.075;
%     xx = Snode(:,1);
%     yy = Snode(:,2);
%     Snode(:,1) = xx + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
%     Snode(:,2) = yy + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
%     
%     SmeshInfo = polyMeshAuxStructure(Snode, Selem);
%     node = SmeshInfo.node; elem = SmeshInfo.elem;
%     save(SmeshName, 'node', 'elem')
%  
%     %--- Darcy Domain [0,1]x[0,1]
%     DmeshName = ['Dmesh_contortedDualTri_[0,1]x[0,1]_',meshN];
%     %[nodeOld,elemOld] = squarequadmesh([0,1,0,1],1/2^(nn+1)); 
%     [nodeOld, elemOld] = generate_Tri_P_T(0,1,0,1,[1/2,1/2].^(nn+1));
%     [Dnode, Delem] = dualmesh(nodeOld,elemOld);
%     %--- contorted the mesh
%     c_theta = 0.075;
%     xx = Dnode(:,1);
%     yy = Dnode(:,2);
%     Dnode(:,1) = xx + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
%     Dnode(:,2) = yy + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
%     
%     DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
%     node = DmeshInfo.node; elem = DmeshInfo.elem;
%     save(DmeshName, 'node', 'elem')
%    
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
%     figure
%     PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
% end % nn


%%
% %------ mesh5
% %-------------------- Tri Mesh ---------------------
% for nn = 1:4
%     close all;
%     meshN = num2str(2^(nn+1));
%     %--- Stokes Domain [0,1]x[1,2]
%     SmeshName = ['Smesh_tri_[0,1]x[1,2]_',meshN];
%     %[nodeOld,elemOld] = squarequadmesh([0,1,1,2],1/2^(nn+1)); 
%     [Snode, Selem] = generate_Tri_P_T(0,1,1,2,[1/2,1/2].^(nn+1));
%     
% %     %--- contorted the mesh
% %     xx = Snode(:,1);
% %     yy = Snode(:,2);
% %     Snode(:,2) = yy.*((-4*abs(xx-0.5)+2).*(yy-2).^2.*(yy-1).^2+1);
%     
%     SmeshInfo = polyMeshAuxStructure(Snode, Selem);
%     node = SmeshInfo.node; elem = SmeshInfo.elem;
%     save(SmeshName, 'node', 'elem')
%  
%     %--- Darcy Domain [0,1]x[0,1]
%     DmeshName = ['Dmesh_tri_[0,1]x[0,1]_',meshN];
%     %[nodeOld,elemOld] = squarequadmesh([0,1,0,1],1/2^(nn+1)); 
%     [Dnode, Delem] = generate_Tri_P_T(0,1,0,1,[1/2,1/2].^(nn+1));
%     DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
%     node = DmeshInfo.node; elem = DmeshInfo.elem;
%     save(DmeshName, 'node', 'elem')
%    
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
%     figure
%     PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
% end % nn


%%
% %------ mesh6
% %-------------------- Quad Mesh ---------------------
% for nn = 1:4
%     close all;
%     meshN = num2str(2^(nn+1));
%     %--- Stokes Domain [0,1]x[1,2]
%     SmeshName = ['Smesh_quad_[0,1]x[1,2]_',meshN];
%     [Snode,Selem] = squarequadmesh([0,1,1,2],1/2^(nn+1)); 
%     %[Snode, Selem] = generate_quad_P_T(0,1,1,2,[1/2,1/2].^(nn+1));
%     
% %     %--- contorted the mesh
% %     xx = Snode(:,1);
% %     yy = Snode(:,2);
% %     Snode(:,2) = yy.*((-4*abs(xx-0.5)+2).*(yy-2).^2.*(yy-1).^2+1);
%     col_elem = Selem(:,4); Selem(:,4) = Selem(:,2); Selem(:,2) = col_elem;
%     SmeshInfo = polyMeshAuxStructure(Snode, Selem);
%     node = SmeshInfo.node; elem = SmeshInfo.elem;
%     
%     save(SmeshName, 'node', 'elem')
%  
%     %--- Darcy Domain [0,1]x[0,1]
%     DmeshName = ['Dmesh_quad_[0,1]x[0,1]_',meshN];
%     [Dnode,Delem] = squarequadmesh([0,1,0,1],1/2^(nn+1)); 
%     %[Dnode, Delem] = generate_Tri_P_T(0,1,0,1,[1/2,1/2].^(nn+1));
%     col_elem = Delem(:,4); Delem(:,4) = Delem(:,2); Delem(:,2) = col_elem;
%     DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
%     node = DmeshInfo.node; elem = DmeshInfo.elem;
%     save(DmeshName, 'node', 'elem')
%    
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
%     figure
%     PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
% end % nn


%%
% % %------ mesh7
% % %-------------------- Contorted Quad Mesh ---------------------
% for nn = 0:4
%     close all;
%     meshN = num2str(2^(nn+1));
%     %--- Stokes Domain [0,1]x[1,2]
%     SmeshName = ['Smesh_contortedQuad_[0,1]x[1,2]_',meshN];
%     [Snode,Selem] = squarequadmesh([0,1,1,2],1/2^(nn+1)); 
%     col_elem = Selem(:,4); Selem(:,4) = Selem(:,2); Selem(:,2) = col_elem;
%     %[Snode, Selem] = generate_Tri_P_T(0,1,1,2,[1/2,1/2].^(nn+1));
%     %--- contorted the mesh
%     xx = Snode(:,1);
%     yy = Snode(:,2);
%     mdxx = contortedX(xx, 1/2^(nn+1));
%     Snode(:,2) = yy + yy.*((1*mdxx).*(yy-2).^2.*(yy-1).^2);
%   
%     SmeshInfo = polyMeshAuxStructure(Snode, Selem);
%     node = SmeshInfo.node; elem = SmeshInfo.elem;
%     save(SmeshName, 'node', 'elem')
%  
%     %--- Darcy Domain [0,1]x[0,1]
%     DmeshName = ['Dmesh_contortedQuad_[0,1]x[0,1]_',meshN];
%     [Dnode,Delem] = squarequadmesh([0,1,0,1],1/2^(nn+1)); 
%     col_elem = Delem(:,4); Delem(:,4) = Delem(:,2); Delem(:,2) = col_elem;
%     %[Dnode, Delem] = generate_Tri_P_T(0,1,0,1,[1/2,1/2].^(nn+1));
%     %--- contorted the mesh
%     xx = Dnode(:,1);
%     yy = Dnode(:,2);
%     mdxx = contortedX(xx, 1/2^(nn+1));
%     Dnode(:,2) = yy + (yy+1).*((1*mdxx).*(yy-1).^2.*(yy-0).^2);
%     
%     DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
%     node = DmeshInfo.node; elem = DmeshInfo.elem;
%     save(DmeshName, 'node', 'elem')
%    
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
%     figure
%     PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
% end % nn


%%
% %------ mesh8
% %-------------------- PolyMesher-Polygon Mesh ---------------------
% for nn = 2:4
%     close all;
%     meshN = num2str(2^(nn+1));
%     SmeshName = ['Smesh_poly_[0,1]x[1,2]_',meshN];
%     DmeshName = ['Dmesh_poly_[0,1]x[0,1]_',meshN];
%     %--- Stokes Domain [0,1]x[1,2]
%     redo_StokesDomain = 1;
%     redo_DarcyDomain = 1;
%     
%     if redo_StokesDomain == 0
%         load(SmeshName);
%         SmeshInfo = polyMeshAuxStructure(node, elem);
%     end
%     if redo_DarcyDomain == 0
%         load(DmeshName);
%         DmeshInfo = polyMeshAuxStructure(node, elem);
%     end
%         
%     
%     if redo_StokesDomain
%         [Snode,Selem]=PolyMesher_toSDmesh(1/2^(nn+1),@poly_StokesDomain,4^(nn+1),1699);
%         Syy = Snode(:,2);
%         y_index = 1:length(Syy);
%         y_index = y_index(abs(Syy-1)<6e-8);
%         x_1 = Snode(y_index,1);
%         [~,b_2] = sort(x_1);
%         %Snode(y_index(b_2),1) = 0:1/2^(nn+1):1;
% 
%         SmeshInfo = polyMeshAuxStructure(Snode, Selem);
%         node = SmeshInfo.node; elem = SmeshInfo.elem;
%         save(SmeshName, 'node', 'elem')
%     end
% 
%     %--- Darcy Domain [0,1]x[0,1]
%     if redo_DarcyDomain
%         [Dnode,Delem]=PolyMesher_toSDmesh(1/2^(nn+1),@poly_DarcyDomain,4^(nn+1),1699);
%         Dyy = Dnode(:,2); 
%         y_index = 1:length(Dyy);
%         y_index = y_index(abs(Dyy-1)<6e-8);
%         x_1 = Dnode(y_index,1);
%         [~,b_2] = sort(x_1);
%         %Dnode(y_index(b_2),1) = 0:1/2^(nn+1):1;
% 
%         DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
%         node = DmeshInfo.node; elem = DmeshInfo.elem;
%         save(DmeshName, 'node', 'elem')
%     end
%    
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,SmeshInfo.Nelems,...
%         DmeshInfo.node,DmeshInfo.elem,DmeshInfo.Nelems)
%     figure
%     PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,SmeshInfo.Nelems,...
%         DmeshInfo.node,DmeshInfo.elem,DmeshInfo.Nelems)
% end % nn


%%
% %------ mesh 9
% %-------------------- PolyMesher-symmetryPolygon Mesh ---------------------
% for nn = 1:4
%     close all;
%     meshN = num2str(2^(nn+1));
%     SmeshName = ['Smesh_symmetryPoly_[0,1]x[1,2]_',meshN];
%     DmeshName = ['Dmesh_symmetryPoly_[0,1]x[0,1]_',meshN];
%     
%     [Snode,Selem]=PolyMesher_toSDmesh(1/2^(nn+1),@poly_StokesDomain,4^(nn+1),1699);
%     
%     left_bound = 0; right_bound = 1;
%     bottom_bound = 1; top_bound = 2;
%     Snode = modifiedNode(Snode,left_bound,right_bound,bottom_bound,top_bound);
%     SmeshInfo = polyMeshAuxStructure(Snode, Selem);
%     node = SmeshInfo.node; elem = SmeshInfo.elem;
%     save(SmeshName, 'node', 'elem')
%     
%     Dnode = zeros(size(Snode));
%     Dnode(:,1) = Snode(:,1); Dnode(:,2) = 2 - Snode(:,2); 
%     Delem = cell(length(Selem),1);
%     for nt = 1:length(Delem)
%         singleElem = Selem{nt};
%         singleNE = length(singleElem);
%         singleElem_1 = zeros(size(singleElem));
%         for ii = 1:length(singleElem)
%             singleElem_1(ii) = singleElem(singleNE-(ii-1));
%         end
%             Delem{nt} = singleElem_1;
%     end 
%     
%     DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
%     node = DmeshInfo.node; elem = DmeshInfo.elem;
%     save(DmeshName, 'node', 'elem')
%     
%     
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,SmeshInfo.Nelems,...
%         DmeshInfo.node,DmeshInfo.elem,DmeshInfo.Nelems)
%     figure
%     PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,SmeshInfo.Nelems,...
%         DmeshInfo.node,DmeshInfo.elem,DmeshInfo.Nelems)
% end


%%
% %------ mesh 10
% %-------------------------- Rectilinear Quad Mesh ---------------------------
% %--- (rectilinear, 直线运动的)
% for nn = 1:4
%     close all
%     meshN = num2str(2^(nn+1));
%     SmeshName = ['Smesh_rectilinearQuad_[0,1]x[1,2]_',meshN];
%     DmeshName = ['Dmesh_rectilinearQuad_[0,1]x[0,1]_',meshN];
%     
%     %--- set h
%     h = 1/(2^(nn+1));
%     
%     %--- Darcy mesh
%     % (1). --- mrst ---
%     dx =1 - (0.5) * cos((-1:(2*h):1) * pi);
%     dx = cumsum(dx);
%     maxdx = max(dx);
%     mindx = min(dx);
%     x = dx/(maxdx-mindx) + (-mindx)/(maxdx-mindx);
%     
%     %y = 0:h:1;
%     %G = tensorGrid(x, sqrt(y));
%     dy =1 - (0.5) * sin((-1:(2*h):1) * pi);
%     dy = cumsum(dy);
%     dy = sqrt(dy);
%     maxdy = max(dy);
%     mindy = min(dy);
%     y = dy/(maxdy-mindy) + (-mindy)/(maxdy-mindy);
%     G = tensorGrid(x, y);
%     %G = tensorGrid(x, sqrt(y));
%     G.nodes.coords = twister(G.nodes.coords);
%     plotGrid(G); axis([-0.05 1.05 -0.05 1.05]);
%     
%     % (2). --- get node and elem ---
%     Dnode = G.nodes.coords;
%     Delem = cell(G.cells.num,1);
%     for nElemD = 1:G.cells.num
%         faceIndex = G.cells.faces(G.cells.facePos(nElemD) : G.cells.facePos(nElemD+1)-1, 1);
%         
%         nodeIndex = [];
%         for ii = 1:length(faceIndex)
%             currNodeIndex = ...
%                 G.faces.nodes(G.faces.nodePos(faceIndex(ii)):G.faces.nodePos(faceIndex(ii)+1)-1);
%             nodeIndex = [currNodeIndex;
%                 nodeIndex];
%         end % for ii
%         currElemD = unique(nodeIndex','stable');
%         
%         coordXY = Dnode(currElemD,:);
%         K = convhull(coordXY(:,1),coordXY(:,2));
%         Delem{nElemD} = currElemD(K(1:end-1));
%         
%     end % nElemD
%     
%     DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
%     %plotPolyMsh(DmeshInfo);
%     node = DmeshInfo.node; elem = DmeshInfo.elem;
%     save(DmeshName, 'node', 'elem')
%     
%     %--- Stokes mesh
%     Snode = zeros(size(Dnode));
%     Snode(:,1) = Dnode(:,1); Snode(:,2) = 2 - Dnode(:,2); 
%     Selem = cell(length(Delem),1);
%     for nElemS = 1:length(Selem)
%         currElemS = Delem{nElemS};
%         coordXY = Snode(currElemS,:);
%         K = convhull(coordXY(:,1),coordXY(:,2));
%         Selem{nElemS} = currElemS(K(1:end-1));
%     end
%     SmeshInfo = polyMeshAuxStructure(Snode, Selem);
%     %plotPolyMsh(SmeshInfo);
%     node = SmeshInfo.node; elem = SmeshInfo.elem;
%     save(SmeshName, 'node', 'elem')
% 
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,SmeshInfo.Nelems,...
%         DmeshInfo.node,DmeshInfo.elem,DmeshInfo.Nelems)
%     figure
%     PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,SmeshInfo.Nelems,...
%         DmeshInfo.node,DmeshInfo.elem,DmeshInfo.Nelems)
%     
%     interfaceX = (0:0.1:1)';
%     interfaceY = ones(length(interfaceX),1);
%     plot(interfaceX,interfaceY,'-r','LineWidth',1.5);
% end % for nn



%%
%------- mesh 11
%---------------------- Delaunay Mesh ----------------------
fd=@(p) drectangle(p,0,1,0,1);
[p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);

%------------------------------------------------------------------




%-----------------------------------------------------------------------------

end


%%------- sub function
%
%--- sub function 1
function PolyMshr_PlotMsh_twoDomain1(Snode,Selem,SNelem,Dnode,Delem,DNelem)
Sxx = Snode(:,1);
Syy = Snode(:,2);

hold on;
Selem = Selem(1:SNelem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Selem));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Selem,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Snode,'FaceColor','w'); pause(1e-6)

hold on
Dxx = Dnode(:,1);
Dyy = Dnode(:,2);
Delem = Delem(1:DNelem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Delem));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Delem,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Dnode,'FaceColor','w'); pause(1e-6)

axis equal; axis([min([Sxx;Dxx])-0.1 max([Sxx;Dxx])+0.1 min([Syy;Dyy])-0.1 max([Syy;Dyy])+0.1])
end


%--- sub function 2
function PolyMshr_PlotMsh_twoDomain2(Snode,Selem,SNelem,Dnode,Delem,DNelem)
Sxx = Snode(:,1);
Syy = Snode(:,2);

hold on;
Selem = Selem(1:SNelem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Selem));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Selem,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Snode,'FaceColor','w'); pause(1e-6)

hold on
Dxx = Dnode(:,1);
Dyy = Dnode(:,2);
Delem = Delem(1:DNelem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Delem));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Delem,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Dnode,'FaceColor','w'); pause(1e-6)

axis equal; 
axis off;
%axis([min([Sxx;Dxx])-0.1 max([Sxx;Dxx])+0.1 min([Syy;Dyy])-0.1 max([Syy;Dyy])+0.1])
end


%--- sub function 3
function node = modifiedNode(node,left_bound,right_bound,bottom_bound,top_bound)

xx = node(:,1);
yy = node(:,2);

tol = 6e-9;

left_abs_xx = abs( xx - left_bound );
right_abs_xx = abs( xx - right_bound );
bottom_abs_yy = abs( yy - bottom_bound );
top_abs_yy = abs( yy - top_bound );

node(left_abs_xx<tol,1) = left_bound;
node(right_abs_xx<tol,1) = right_bound;
node(bottom_abs_yy<tol,2) = bottom_bound;
node(top_abs_yy<tol,2) = top_bound;

end % function modifiedNode
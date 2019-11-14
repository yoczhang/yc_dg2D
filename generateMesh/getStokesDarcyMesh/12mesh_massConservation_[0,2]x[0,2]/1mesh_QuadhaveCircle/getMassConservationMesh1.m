%%
%------ mesh 12
%---------------------- Darcy have Faults Mesh ----------------------
function getMassConservationMesh1()
close all
clearvars
clc

meshN = num2str(32);
%--- Stokes Domain [0,1]x[1,2]
SmeshName = ['Smesh_1massConservation_[0,2]x[1,2]_',meshN];
DmeshName = ['Dmesh_1massConservation_[0,2]x[0,1]_',meshN];
finalFigName = ['SD_1massConservation_(0,1)x(1,2)_',meshN,'_backup'];

if (true)
    if (true) % choose to get the Darcy mesh
        load(['meshMassConservation_',meshN,'.mat'])
        % Plot pebiGrid
        figure(); hold on
        %plotGrid(G,'facecolor','none')
        % plotGrid(G,G.cells.tag, 'facecolor','b')
        % centF = G.faces.centroids(G.faces.tag,:);
        % plot(centF(:,1), centF(:,2),'.','markersize',10)
        % axis equal tight
        % title('pebiGrid(...)')
        
        %-- remove short edges
        disp(['%#-- original G.faces.num = ',num2str(G.faces.num)]);
        G = removeShortEdges(G, 8e-4);
        disp(['%#-- after remove short edges: G.faces.num = ',num2str(G.faces.num)]);
        plotGrid(G,'facecolor','none')
        %--- next to get the Darcy meshInfo
        [Dnode,Delem] = mrstG_2_myMeshInfo(G);
        DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
        patchPlotMesh(DmeshInfo.node,DmeshInfo.elem);

        %-- retreat the mesh, aim to fixed the boudary points
        %load 'DmeshInfo'
        %patchPlotMesh(DmeshInfo.node,DmeshInfo.elem);

        bdNodeIndex = DmeshInfo.bdEdge(:);
        bdNodeIndex = unique(bdNodeIndex(:));
        bdNodeIndex_temp = (1:length(bdNodeIndex))';
        bdCoord = [DmeshInfo.node(bdNodeIndex,1),DmeshInfo.node(bdNodeIndex,2)];
        plot(bdCoord(:,1),bdCoord(:,2),'or','MarkerSize',4);
        bdEdgeIndex = DmeshInfo.bdEdgeIndex;
        barybdCoord = DmeshInfo.baryEdge(bdEdgeIndex,:);
        plot(barybdCoord(:,1),barybdCoord(:,2),'sb','MarkerSize',4);

        fixedNodeIndex_x0 = bdNodeIndex_temp(abs(bdCoord(:,1)-0)<1e-8);
        DmeshInfo.node(bdNodeIndex(fixedNodeIndex_x0),1) = 0;
        fixedNodeIndex_x1 = bdNodeIndex_temp(abs(bdCoord(:,1)-1)<1e-2);
        DmeshInfo.node(bdNodeIndex(fixedNodeIndex_x1),1) = 1;
        fixedNodeIndex_y0 = bdNodeIndex_temp(abs(bdCoord(:,2)-0)<1e-8);
        DmeshInfo.node(bdNodeIndex(fixedNodeIndex_y0),2) = 0;
        fixedNodeIndex_y1 = bdNodeIndex_temp(abs(bdCoord(:,2)-1)<1e-8);
        DmeshInfo.node(bdNodeIndex(fixedNodeIndex_y1),2) = 1;
        DmeshInfo = polyMeshAuxStructure(DmeshInfo.node,DmeshInfo.elem);
        patchPlotMesh(DmeshInfo.node,DmeshInfo.elem);

        bdNodeIndex = DmeshInfo.bdEdge(:);
        bdNodeIndex = unique(bdNodeIndex(:));
        bdCoord = [DmeshInfo.node(bdNodeIndex,1),DmeshInfo.node(bdNodeIndex,2)];
        plot(bdCoord(:,1),bdCoord(:,2),'or','MarkerSize',4);
        bdEdgeIndex = DmeshInfo.bdEdgeIndex;
        barybdCoord = DmeshInfo.baryEdge(bdEdgeIndex,:);
        plot(barybdCoord(:,1),barybdCoord(:,2),'sb','MarkerSize',4);

        % plotPolyMsh(DmeshInfo);
        % node = DmeshInfo.node; elem = DmeshInfo.elem;
        % save(DmeshName, 'node', 'elem')
        %save DmeshInfo DmeshInfo
    end 
    
    %--- then to get the Stokes meshInfo
    % we to get the (interface) matched mesh.
    %load 'DmeshInfo'
    interface_Y = 1;
    temp_1 = (1:DmeshInfo.Nnodes)';
    interface_nodeIndex = temp_1(abs(DmeshInfo.node(:,2)-interface_Y)<1e-8);
    interface_X = DmeshInfo.node(interface_nodeIndex,1);
    
    x_partition = sort(interface_X);
    y_partition = 1:2/(length(x_partition)+0):2;
    [x,y] = meshgrid(x_partition,y_partition);
    
    [Snode,Selem] = yc_squarequadmesh(x,y);
    % %--col_elem = Selem(:,4); Selem(:,4) = Selem(:,2); Selem(:,2) = col_elem;
    SmeshInfo = polyMeshAuxStructure(Snode, Selem);
    %patchPlotMesh(testmeshInfo.node,testmeshInfo.elem);
    
%     G_Smesh = triangleGrid([x(:) y(:)]);
%     [Snode,Selem] = mrstG_2_myMeshInfo(G_Smesh);
%     SmeshInfo = polyMeshAuxStructure(Snode, Selem);
%     %patchPlotMesh(SmeshInfo.node,SmeshInfo.elem);
    
    finalFig = figure;
    PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,SmeshInfo.Nelems,...
        DmeshInfo.node,DmeshInfo.elem,DmeshInfo.Nelems)
    
    interface_Y = ones(length(interface_X),1);
    %plot(interface_X,interface_Y,'or','MarkerSize',4);
    plot(interface_X,interface_Y,'-r','LineWidth',1);
    
    %saveas(finalFig,finalFigName,'fig')
    %saveas(finalFig,finalFigName,'eps')
    
    %-- plot the boundary
%     coordx1 = 1;
%     Nnode_temp = (1:DmeshInfo.Nnodes)';
%     fixed_x1 = Nnode_temp((abs(DmeshInfo.node(:,1)-coordx1)<1e-6));
%     fixed_coordy = DmeshInfo.node(fixed_x1,2);
%     xx = coordx1*ones(length(fixed_coordy),1);
%     plot(xx,fixed_coordy,'or','MarkerSize',4);
    
end % if false
node = DmeshInfo.node; elem = DmeshInfo.elem;
save(DmeshName, 'node', 'elem')
node = SmeshInfo.node; elem = SmeshInfo.elem;
save(SmeshName, 'node', 'elem')

%% other plot

% figure; hold on
% xx = [0,2,2,0,0];
% yy = [0,0,2,2,0];
% plot(xx,yy,'-k');
% plot([0,2],[1,1],'-r');
% r = exp(-4.1:0.35:-1.4);
% r0 = r(end-1);
% c1 = [0.3, 0.7]; c2 = [1.3, 0.4];
% plot(c1(1),c1(2),'.k','MarkerSize',4)
% plot(c2(1),c2(2),'.k','MarkerSize',4)
% rectangle('Position',[c1(1)-r0,c1(2)-r0,2*r0,2*r0] ,'Curvature',[1,1])
% rectangle('Position',[c2(1)-r0,c2(2)-r0,2*r0,2*r0] ,'Curvature',[1,1])
% axis off; axis equal
end % function


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

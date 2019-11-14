%%
%------ mesh 12
%---------------------- Darcy have Faults Mesh ----------------------
function getDarcyhaveOnlyFaultsMesh1()
close all
clearvars
clc

meshN = num2str(1);
%--- Stokes Domain [0,1]x[1,2]
SmeshName = ['Smesh_DarcyhaveFault_[0,1]x[1,2]_',meshN];
DmeshName = ['Dmesh_DarcyhaveFault_[0,1]x[0,1]_',meshN];
finalFigName = ['SD_meshHaveFault_(0,1)x(1,2)_',meshN,'_backup'];

if (true)
    if (true) % choose to get the Darcy mesh
        %load('G_fault_well_canbeused3')
        load('G_onefault')
        G = G_fault;
        
        faultFacesIndx = find(G.faces.tag);
        
        %- Plot pebiGrid
        figure(); hold on
        plotGrid(G,'facecolor','none')
        plotGrid(G,G.cells.tag, 'facecolor','b')
        centF = G.faces.centroids(G.faces.tag,:);
        plot(centF(:,1), centF(:,2),'.','markersize',10)
        axis equal tight
        title('pebiGrid(...)')
        
        %--- yc to get the cells
        faultNeighborCells = G.faces.neighbors(G.faces.tag,:);
%         [nrow,ncol] = size(faultNeighborCells);
%         faultNeighborCells = reshape(faultNeighborCells,nrow*ncol,1);
%         faultNeighborCells = faultNeighborCells(~(faultNeighborCells==0));
%         faultNeighborCells = unique(faultNeighborCells);
%         faultCells_tag = false(G.cells.num,1);
%         faultCells_tag(faultNeighborCells) = true;
        
%         figure
%         plotGrid(G,faultCells_tag, 'facecolor','b')
%         plotGrid(G,'FaceColor', 'none')
%         centF = G.faces.centroids(G.faces.tag,:);
%         plot(centF(:,1), centF(:,2),'.','markersize',5)
%         axis off
%         axis equal tight
        %title('pebiGrid(original)')
       
        G = removeShortEdges(G, 8e-8);
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
        
        %-- get the DmeshInfo faluts faces
        [nrow,~] = size(faultNeighborCells);
        if length(faultFacesIndx) == nrow
            D_faultFaces = zeros(nrow,1);
            for ii = 1:nrow
                NeiCell1 = faultNeighborCells(ii,1);
                NeiCell2 = faultNeighborCells(ii,2);
                faultsNeiCell1 = DmeshInfo.elem2edge{NeiCell1};
                faultsNeiCell2 = DmeshInfo.elem2edge{NeiCell2};
                D_faultFaces(ii) = intersect(faultsNeiCell1,faultsNeiCell2);
                
                baryFaultsCoord = DmeshInfo.baryEdge(D_faultFaces(ii),:);
                plot(baryFaultsCoord(:,1),baryFaultsCoord(:,2),'*b','MarkerSize',4);
            end
            DmeshInfo.faultFaces = D_faultFaces;
        else
            error('length(faultFacesIndx) != nrow.')
        end
        
        save DmeshInfo DmeshInfo
    end 
    
    %--- then to get the Stokes meshInfo
    % we to get the (interface) matched mesh.
    %load 'DmeshInfo'
    interface_Y = 1;
    temp_1 = (1:DmeshInfo.Nnodes)';
    interface_nodeIndex = temp_1(abs(DmeshInfo.node(:,2)-interface_Y)<1e-8);
    interface_X = DmeshInfo.node(interface_nodeIndex,1);
    
    x_partition = sort(interface_X);
    y_partition = 1:1/(length(x_partition)+1):2;
    [x,y] = meshgrid(x_partition,y_partition);
    
    [Snode,Selem] = yc_squarequadmesh(x,y);
    % need to check the Quad mesh node indices in one elem is anticlockwise
    % or clockwise, if is the clockwise, we need to rearrange the indices.
    % % col_elem = Selem(:,4); Selem(:,4) = Selem(:,2); Selem(:,2) = col_elem;
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
    
    saveas(finalFig,finalFigName,'fig')
    saveas(finalFig,finalFigName,'eps')
    
    %-- plot the boundary
%     coordx1 = 1;
%     Nnode_temp = (1:DmeshInfo.Nnodes)';
%     fixed_x1 = Nnode_temp((abs(DmeshInfo.node(:,1)-coordx1)<1e-6));
%     fixed_coordy = DmeshInfo.node(fixed_x1,2);
%     xx = coordx1*ones(length(fixed_coordy),1);
%     plot(xx,fixed_coordy,'or','MarkerSize',4);
    
end % if false
% node = DmeshInfo.node; elem = DmeshInfo.elem;
% save(DmeshName, 'node', 'elem', 'faultCells', 'wellCells')
% node = SmeshInfo.node; elem = SmeshInfo.elem;
% save(SmeshName, 'node', 'elem')

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

function main_getOnlyFaultsMeshInfo_case2()
%
%
%
%

close all
clearvars
clc

%%
[G_fault, originalFracture] = getFaultsMeshFromMRST_case2();
save G_mesh_backup2 G_fault originalFracture

NorigFracture = length(originalFracture);
G = G_fault;
splitedFraCells = splitedFractureFromMRST(G,originalFracture);

%- Plot pebiGrid
% figure(); hold on
% plotGrid(G,'facecolor','none')
% plotGrid(G,G.cells.tag, 'facecolor','b')
% centF = G.faces.centroids(G.faces.tag,:);
% plot(centF(:,1), centF(:,2),'.','markersize',10)
% axis equal tight
% title('pebiGrid(...)')

%--- next to get the Darcy meshInfo
[Dnode,Delem] = mrstG_2_myMeshInfo(G);
DmeshInfo = polyMeshAuxStructure(Dnode, Delem);

patchPlotMesh(DmeshInfo.node,DmeshInfo.elem);


%%
% %- plot the boundary nodes
bdNodeIndex = DmeshInfo.bdEdge(:);
bdNodeIndex = unique(bdNodeIndex(:));
bdNodeIndex_temp = (1:length(bdNodeIndex))';
bdCoord = [DmeshInfo.node(bdNodeIndex,1),DmeshInfo.node(bdNodeIndex,2)];
% bdEdgeIndex = DmeshInfo.bdEdgeIndex;
% barybdCoord = DmeshInfo.baryEdge(bdEdgeIndex,:);
% plot(bdCoord(:,1),bdCoord(:,2),'or','MarkerSize',4);
% plot(barybdCoord(:,1),barybdCoord(:,2),'sb','MarkerSize',4);

%- fixed the boundary nodes
fixedNodeIndex_x0 = bdNodeIndex_temp(abs(bdCoord(:,1)-0)<1e-8);
DmeshInfo.node(bdNodeIndex(fixedNodeIndex_x0),1) = 0;
fixedNodeIndex_x1 = bdNodeIndex_temp(abs(bdCoord(:,1)-1)<1e-2);
DmeshInfo.node(bdNodeIndex(fixedNodeIndex_x1),1) = 1;
fixedNodeIndex_y0 = bdNodeIndex_temp(abs(bdCoord(:,2)-0)<1e-8);
DmeshInfo.node(bdNodeIndex(fixedNodeIndex_y0),2) = 0;
fixedNodeIndex_y1 = bdNodeIndex_temp(abs(bdCoord(:,2)-1)<1e-8);
DmeshInfo.node(bdNodeIndex(fixedNodeIndex_y1),2) = 1;
DmeshInfo = polyMeshAuxStructure(DmeshInfo.node,DmeshInfo.elem);
% patchPlotMesh(DmeshInfo.node,DmeshInfo.elem);

% %- plot the boundary nodes
% bdNodeIndex = DmeshInfo.bdEdge(:);
% bdNodeIndex = unique(bdNodeIndex(:));
% bdCoord = [DmeshInfo.node(bdNodeIndex,1),DmeshInfo.node(bdNodeIndex,2)];
% bdEdgeIndex = DmeshInfo.bdEdgeIndex;
% barybdCoord = DmeshInfo.baryEdge(bdEdgeIndex,:);
% plot(bdCoord(:,1),bdCoord(:,2),'or','MarkerSize',4);
% plot(barybdCoord(:,1),barybdCoord(:,2),'sb','MarkerSize',4);

%%
%- get the DmeshInfo splited fractures faces index
splitedFractureFace = cell(NorigFracture,1);
fractureFace = [];
for ii = 1:NorigFracture
    currOrigFrac_cells = splitedFraCells{ii};
    [nrow,~] = size(currOrigFrac_cells);
    D_faultFaces = zeros(nrow,1);
    for kk = 1:nrow
        NeiCell1 = currOrigFrac_cells(kk,1);
        NeiCell2 = currOrigFrac_cells(kk,2);
        faultsNeiCell1 = DmeshInfo.elem2edge{NeiCell1};
        faultsNeiCell2 = DmeshInfo.elem2edge{NeiCell2};
        D_faultFaces(kk) = intersect(faultsNeiCell1,faultsNeiCell2);
        
        baryFaultsCoord = DmeshInfo.baryEdge(D_faultFaces(kk),:);
        plot(baryFaultsCoord(:,1),baryFaultsCoord(:,2),'*b','MarkerSize',4);
    end
    splitedFractureFace{ii} = D_faultFaces;
    temp = fractureFace;
    fractureFace = [temp; D_faultFaces]; 
end

%%
DmeshInfo.originalFracture = originalFracture;
DmeshInfo.splitedFractureFace = splitedFractureFace;
DmeshInfo.fractureFace = fractureFace;

%%
save DmeshInfo DmeshInfo


end % function getOnlyFaultsMeshInfo
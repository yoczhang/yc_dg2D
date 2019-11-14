function splitedFraCells = splitedFractureFromMRST(G,originalFracture)
% % close all
% % clearvars
% % clc

addpath(genpath([ROOTDIR,'modules/upr/'])); 

%-
faultNeighborCells = G.faces.neighbors(G.faces.tag,:);
Nfaults = size(faultNeighborCells,1);

% % figure(); hold on
% % plotGrid(G_fault,'facecolor','none')
% % centF = G_fault.faces.centroids(G_fault.faces.tag,:);
% % plot(centF(:,1), centF(:,2),'.','markersize',10)
% % axis equal tight

NoriginalFracture = length(originalFracture);
splitedFraCells = cell(NoriginalFracture,1);
for jj = 1:NoriginalFracture
    oF = originalFracture{jj};
    
    for ii = 1:Nfaults
        cell1 = faultNeighborCells(ii,1);
        cell2 = faultNeighborCells(ii,2);
        
        cell1_coord = G.cells.centroids(cell1,:);
        cell2_coord = G.cells.centroids(cell2,:);
        
        [~, fCut, ~] = splitAtInt({[cell1_coord;cell2_coord]},{oF});
        
        if length(fCut) > 1
            splitedFraCells{jj} = [splitedFraCells{jj};[cell1, cell2]];
        end
    end
    
    faultNeighborCells = setdiff(faultNeighborCells,splitedFraCells{jj},'rows');
    Nfaults = size(faultNeighborCells,1);
    
% %     sFcells = splitedF{jj};
% %     cellcoord1 = G.cells.centroids(sFcells(:,1),:);
% %     cellcoord2 = G.cells.centroids(sFcells(:,2),:);
% %     if jj ==3
% %         plot(cellcoord1(:,1),cellcoord1(:,2),'*r','MarkerSize',4);
% %         plot(cellcoord2(:,1),cellcoord2(:,2),'*g','MarkerSize',4);
% %     else
% %         plot(cellcoord1(:,1),cellcoord1(:,2),'.r','MarkerSize',4);
% %         plot(cellcoord2(:,1),cellcoord2(:,2),'.r','MarkerSize',4);
% %     end
    
end 


end % function
function sectorMeshTri2Poly
%
%
%
%

close all
clearvars; clc

initH = 0.001;


load('SDF_Dmesh_Sector_poly.mat','DmeshInfo')
load('SDF_Smesh_Sector_poly.mat','SmeshInfo')

%------- mrst --------
Dnode = DmeshInfo.node;
Delem = cell2mat(DmeshInfo.elem);
Dedge = DmeshInfo.edge;

Snode = SmeshInfo.node;
Selem = cell2mat(SmeshInfo.elem);

splitedFractureFace = DmeshInfo.splitedFractureFace;
fracture1 = splitedFractureFace{1};
fracture2 = splitedFractureFace{2};
fracture3 = splitedFractureFace{3};
fracture4 = splitedFractureFace{4};

% G_S = pebi(triangleGrid(Snode,Selem));
% G_D = pebi(triangleGrid(Dnode,Delem));
G_S = triangleGrid(Snode,Selem);
G_D = triangleGrid(Dnode,Delem);


plotGrid(G_D,'FaceColor','none');
hold on
plotGrid(G_S,'FaceColor','none');

%- plot interface
x1_5 = (0:initH:2*sqrt(2)/2)';
y1_5 = sqrt(4-x1_5.^2);
y1_6 = (0:initH:2*sqrt(2)/2)';
x1_6 = sqrt(4-y1_6.^2);
plot(x1_5, y1_5, '-r', 'LineWidth',1.2)
plot(x1_6, y1_6, '-r', 'LineWidth',1.2)



%- plot fractures
for k = 1:length(fracture1)
    nodeInd = Dedge(fracture1(k),:);
    coord = Dnode(nodeInd,:);
    plot(coord(:,1),coord(:,2),'-r', 'LineWidth',1.2)
end
for k = 1:length(fracture2)
    nodeInd = Dedge(fracture2(k),:);
    coord = Dnode(nodeInd,:);
    plot(coord(:,1),coord(:,2),'-r', 'LineWidth',1.2)
end
for k = 1:length(fracture3)
    nodeInd = Dedge(fracture3(k),:);
    coord = Dnode(nodeInd,:);
    plot(coord(:,1),coord(:,2),'-r', 'LineWidth',1.2)
end
for k = 1:length(fracture4)
    nodeInd = Dedge(fracture4(k),:);
    coord = Dnode(nodeInd,:);
    plot(coord(:,1),coord(:,2),'-r', 'LineWidth',1.2)


%---
axis equal;
view(2)
axis off;
print('SectorPolymesh.jpg','-djpeg', '-r750')
print('SectorPolymesh.eps','-depsc')




end 
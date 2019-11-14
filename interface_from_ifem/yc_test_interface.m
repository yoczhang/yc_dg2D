function yc_test_interface

close all
clearvars
clc

%- get the interfacemesh information
pde = interfacedata;
box = [ -1, 1, -1, 1];
h = 1/10;
% phi = @(p) sum(p.^2, 2) - 0.5.^2;
phi = pde.phi;
[node,elem,interfaceData] = interfacemesh(box,phi,h);

%---
%- square elems
sElem = interfaceData.sElem;
NTS = size(sElem,1);
isExteriorSElem = false(NTS,1);
c = (node(sElem(:,1),:) + node(sElem(:,2),:) + node(sElem(:,3),:)+node(sElem(:,4),:))/4;
isExteriorSElem(pde.phi(c)>0) = true;

%- triangle elems
tElem = interfaceData.tElem;
NTT = size(tElem,1);
isExteriorTElem = false(NTT,1);
c = (node(tElem(:,1),:) + node(tElem(:,2),:) + node(tElem(:,3),:))/3;
isExteriorTElem(pde.phi(c)>0) = true;

%- interface edges
interfaceEdge = interfaceData.interfaceEdge;

interfaceEdge2 =  findinterfaceedge_new(node, sElem, tElem, interfaceData.DomainEdge, pde.phi );

%%
totalNode = node;

[sRow, sCol] = size(sElem);
cell_sElem = mat2cell(sElem, ones(1,sRow), sCol);
[tRow, tCol] = size(tElem);
cell_tElem = mat2cell(tElem, ones(1,tRow), tCol);

totalElem = [cell_sElem; cell_tElem];
interfacemeshInfo = polyMeshAuxStructure(totalNode, totalElem);
% plotPolyMsh(interfacemeshInfo)
patchPlotMesh(interfacemeshInfo.node, interfacemeshInfo.elem)

%-
isExteriorElem = false(interfacemeshInfo.Nelems,1);
isExteriorElem(pde.phi(interfacemeshInfo.baryElem)>0) = true;
% % whichElem = (1:interfacemeshInfo.Nelems)';
% % whichElem = whichElem(isExteriorElem);
ExteriorElem = interfacemeshInfo.elem(isExteriorElem);
patchPlotMesh_forInterface(interfacemeshInfo.node, ExteriorElem)

%-
NinterfaceEdge = size(interfaceEdge,1);
for ii = 1:NinterfaceEdge
    hold on;
    vertex1 = interfaceEdge(ii,1);
    vertex2 = interfaceEdge(ii,2);
    xx1 = node(vertex1,1); yy1 = node(vertex1,2);
    xx2 = node(vertex2,1); yy2 = node(vertex2,2);
    
    plot([xx1, xx2], [yy1, yy2], '-b', 'LineWidth',2);
end


%-
% patchPlotMesh(interfacemeshInfo.node, interfacemeshInfo.elem)
% interfaceEdge_yc = interfaceData.interfaceEdge_yc;
% NinterfaceEdge_yc = size(interfaceEdge_yc,1);
% for ii = 1:NinterfaceEdge_yc
%     hold on;
%     vertex1 = interfaceEdge_yc(ii,1);
%     vertex2 = interfaceEdge_yc(ii,2);
%     xx1 = node(vertex1,1); yy1 = node(vertex1,2);
%     xx2 = node(vertex2,1); yy2 = node(vertex2,2);
%     
%     plot([xx1, xx2], [yy1, yy2], '-b', 'LineWidth',2);
% end

%-
ExteriorElem_bdEdge = findBoundaryFaces(ExteriorElem);
DomainbdEdge = interfaceData.DomainbdEdge;
interfaceEdge = setdiff(ExteriorElem_bdEdge,DomainbdEdge,'rows','stable');
patchPlotMesh(interfacemeshInfo.node, interfacemeshInfo.elem)
NinterfaceEdge = size(interfaceEdge,1);
for ii = 1:NinterfaceEdge
    hold on;
    vertex1 = interfaceEdge(ii,1);
    vertex2 = interfaceEdge(ii,2);
    xx1 = node(vertex1,1); yy1 = node(vertex1,2);
    xx2 = node(vertex2,1); yy2 = node(vertex2,2);
    
    plot([xx1, xx2], [yy1, yy2], '-b', 'LineWidth',2);
end

%- for example, find the location of first interfaceEdge in the totalEdge 
interfaceEdge1 = interfaceEdge(1,:);
[lia,locb] = ismember(interfaceEdge1,interfacemeshInfo.edge,'rows');
if lia
    disp(num2str(locb))
else
    interfaceEdge1_temp = interfaceEdge1;
    interfaceEdge1(1) = interfaceEdge1_temp(2);
    interfaceEdge1(2) = interfaceEdge1_temp(1);
    [~,locb] = ismember(interfaceEdge1,interfacemeshInfo.edge,'rows');
    disp(num2str(locb))
end

end 
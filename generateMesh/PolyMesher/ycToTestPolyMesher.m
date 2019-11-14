function ycToTestPolyMesher


%  [node,elem]=PolyMesher(@MichellDomain,20,10);
close all
clc
clearvars

date = datestr(now,31); 
setpath_pwd = pwd;
saveFilename=[setpath_pwd,'/temp_saveNodeElem/',date,'_','PolyNodeElem_','.mat'];

%% --- using distMesh to get the init mesh
% %--- part 1, ---
% xy_1 = [0,0.2];
% xy_2 = [0.8,0.2];
% xy_3 = [0.8,0.5];
% xy_4 = [1.3,0.5];
% xy_5 = [1.3,1];
% xy_6 = [2,1];
% xy_7 = [2,1.25];
% xy_8 = [0,1.25];
% 
% pv=[xy_1(1), xy_1(2); xy_2(1), xy_2(2); xy_3(1), xy_3(2); xy_4(1), xy_4(2); 
%     xy_5(1), xy_5(2); xy_6(1), xy_6(2); xy_7(1), xy_7(2); xy_8(1), xy_8(2); xy_1(1), xy_1(2); -0.4 -0.5];
% [p,t]=distmesh2d(@dpoly,@huniform,0.01,[-1,-1; 2,1],pv,pv);
% 
% meshInfo = polyMeshAuxStructure(p,t);
% P = meshInfo.baryElem;

%% -------------------------
% figure
% [node,elem]=PolyMesher_test(@SectorDomain,199,99);
% figure
[node,elem]=PolyMesher_test(@ycTestDomain2,699,999,P);
% figure
% [node,elem]=PolyMesher_test(@ycTestDomain3,190,99);
% figure
% [node,elem]=PolyMesher_test(@ycTestDomain4,190,199);
% figure
% [node,elem]=PolyMesher_test(@ycTestDomain5,190,99);
% figure
% [node,elem]=PolyMesher_test(@ycTestDomain6,190,99);
% figure
% [node,elem]=PolyMesher_test(@ycTestDomain7,190,699);
% figure
% [node,elem]=PolyMesher_test(@ycTestDomain8,190,699); %--- failed.

[node,elem]=PolyMesher_test(@MichellDomain,60,399);
% HornDomain
% SuspensionDomain
% MichellDomain
% WrenchDomain


%% --------------------------
NodeElem.node = node;
NodeElem.elem = elem;
save(saveFilename, 'NodeElem');

saveFigurename = [setpath_pwd,'/temp_saveNodeElem/',date,'_','figure_','.fig'];
saveas(gcf,saveFigurename);
end
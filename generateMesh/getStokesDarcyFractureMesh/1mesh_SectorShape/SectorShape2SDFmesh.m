function SectorShape2SDFmesh()

close all
clearvars
clc

%% example 1
doplot = 1;
initH = 0.09;
interfaceFix_x = (0:initH:2)';
interfaceFix_y = sqrt(4-interfaceFix_x.^2);
interfaceFix = [interfaceFix_x, interfaceFix_y];

%%
[~,~,interfaceFix] = StokesDomain(initH,interfaceFix,doplot);

%%
[p2,t2,interfaceFix] = DarcyDomain(initH,interfaceFix,doplot);

%%
[p1,t1,interfaceFix] = StokesDomain(initH,interfaceFix,doplot);

%% save to 
%- fractures
SmeshInfo = polyMeshAuxStructure(p1, t1);
DmeshInfo = polyMeshAuxStructure(p2, t2);

% %------ the fractures --------
% % fixF1x = (0.3:initH:1)';
% % fixF1y = 0.1*fixF1x + 1.25;
% % 
% % fixF2x = (0.2:initH:0.6)';
% % fixF2y = -0.3*fixF2x + 1.8;
% % 
% % fixF3x = (1.15:initH:1.5)';
% % fixF3y = -0.5*fixF3x + 1.;
% % 
% % fixF4y = (0.5:initH:0.9)';
% % fixF4x = 1.3*ones(size(fixF4y));
% %--------------------------

edge0 = DmeshInfo.baryEdge;
Nedge = size(edge0,1);
edgeIndx = 1:Nedge;
%- line 1
% % fixF1x = (0.3:initH:1)';
% % fixF1y = 0.1*fixF1x + 1.25;
line1x = [0.3, 1]';
line1y = 0.1*line1x + 1.25;
P1A = [line1x(1), line1y(1)];
P1B = [line1x(2), line1y(2)];

dis1 = DisPtToLine(edge0,P1A,P1B);
fracture1 = edgeIndx(abs(dis1)<1e-8);

%- line 2
% % fixF2x = (0.2:initH:0.6)';
% % fixF2y = -0.3*fixF2x + 1.8;
line2x = [0.2, 0.6]';
line2y = -0.3*line2x + 1.8;
P2A = [line2x(1), line2y(1)];
P2B = [line2x(2), line2y(2)];

dis2 = DisPtToLine(edge0,P2A,P2B);
fracture2 = edgeIndx(abs(dis2)<2e-8);

%- line 3
% % fixF3x = (1.15:initH:1.5)';
% % fixF3y = -0.5*fixF3x + 1.;
line3x = [1.15, 1.5]';
line3y = -0.5*line3x + 1.;
P3A = [line3x(1), line3y(1)];
P3B = [line3x(2), line3y(2)];

dis3 = DisPtToLine(edge0,P3A,P3B);
fracture3 = edgeIndx(abs(dis3)<2e-8);

%- line 4
% % fixF4y = (0.5:initH:0.9)';
% % fixF4x = 1.3*ones(size(fixF4y));
line4y = [0.5, 0.9]';
line4x = 1.3*ones(size(line4y));
P4A = [line4x(1), line4y(1)];
P4B = [line4x(2), line4y(2)];

dis4 = DisPtToLine(edge0,P4A,P4B);
fracture4 = edgeIndx(abs(dis4)<2e-8);

%- plot the fractures 
Dnode = DmeshInfo.node;
Delem = DmeshInfo.elem;
Dedge = DmeshInfo.edge;
patchPlotMesh(Dnode, Delem);
hold on
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
end

%----------
splitedFractureFace = cell(4,1);
splitedFractureFace{1} = fracture1;
splitedFractureFace{2} = fracture2;
splitedFractureFace{3} = fracture3;
splitedFractureFace{4} = fracture4;

DmeshInfo.splitedFractureFace = splitedFractureFace;

% save SDF_Smesh_Sector_[0,3]x[0,3]_1 SmeshInfo
% save SDF_Dmesh_Sector_[0,3]x[0,3]_1 DmeshInfo
save SDF_Smesh_Sector_poly SmeshInfo
save SDF_Dmesh_Sector_poly DmeshInfo


end % function



%--- sub function1
function [p1,t1,interfaceFix] = StokesDomain(initH,interfaceFix,doplot)
%--- paras
squarDomain1 = [0,0;3,3];
%---
fd1 = @(p) dintersect( drectangle(p,0,3,0,3), ddiff(dcircle(p,0,0,3),dcircle(p,0,0,2)) );
fh1 = @(p) 0.005*ones(size(p,1),1);
pfix1 = [0,0;3,0;3,3;0,3;interfaceFix];
pfix1 = [pfix1; 2,0; 3,0; 0,3; 0,2; ]; % add the four cornors to the fix points
[p1,t1] = distmesh( fd1, fh1, initH, squarDomain1, pfix1 );

%- repaire the boundary nodes
x0ind = (abs(p1(:,1))<1e-6);
p1(x0ind,1) = 0;
y0ind = (abs(p1(:,2))<1e-6);
p1(y0ind,2) = 0;
% [~,x2ind] = min(p1(y0ind,1));
% p1(y0ind(x2ind),1) = 2;



%- find the interface nodes
interface_indx = (abs(p1(:,1).^2 + p1(:,2).^2 - 4) < 1e-6);
interfaceFix = p1(interface_indx,:);

if doplot
    patch( 'vertices', p1, 'faces', t1, 'facecolor', [.9, .9, .9] )
    hold on
    plot(interfaceFix(:,1),interfaceFix(:,2),'or')
    plot(p1(y0ind,1),p1(y0ind,2),'.g')
    axis tight
    axis equal
end

end


%--- sub function21
function [p2,t2,interfaceFix] = DarcyDomain(initH,interfaceFix,doplot)
%--- paras
squarDomain2 = [0,0;2,2];
%---
fd2 = @(p) dintersect( drectangle(p,0,2,0,2), ddiff(dcircle(p,0,0,2),dcircle(p,0,0,1)) );
fh2 = @(p) 0.005*ones(size(p,1),1);
pfix2 = [0,0;2,0;2,2;0,2; interfaceFix];
pfix2 = [pfix2; 1,0; 2,0; 0,2; 0,1];

%- some fixed points for Fractures
fixF1x = (0.3:initH:1)';
% fixF1y = 1.1*ones(size(fixF1x))+0.02;
fixF1y = 0.1*fixF1x + 1.25;
fixF1 = [fixF1x,fixF1y];
pfix2 = [pfix2;fixF1];

fixF2x = (0.2:initH:0.6)';
fixF2y = -0.3*fixF2x + 1.8;
fixF2 = [fixF2x,fixF2y];
pfix2 = [pfix2;fixF2];

%- the intersect Fracture
fixF3x = (1.15:initH:1.5)';
fixF3y = -0.5*fixF3x + 1.;
fixF3 = [fixF3x,fixF3y];
pfix2 = [pfix2;fixF3];

fixF4y = (0.5:initH:0.9)';
fixF4x = 1.3*ones(size(fixF4y));
fixF4 = [fixF4x,fixF4y];
pfix2 = [pfix2;fixF4];

[p2,t2] = distmesh( fd2, fh2, initH, squarDomain2, pfix2 );

%- repaire the boundary nodes
x0ind = (abs(p2(:,1))<1e-6);
p2(x0ind,1) = 0;
y0ind = (abs(p2(:,2))<1e-6);
p2(y0ind,2) = 0;

if doplot
    patch( 'vertices', p2, 'faces', t2, 'facecolor', [.9, .9, .9] )
    hold on
    plot(interfaceFix(:,1),interfaceFix(:,2),'.r')
    plot(fixF1x,fixF1y,'-r')
    plot(fixF2x,fixF2y,'-r')
    plot(fixF3x,fixF3y,'-r')
    plot(fixF4x,fixF4y,'-r')
    axis tight
    axis equal
end

interface_indx = (abs(p2(:,1).^2 + p2(:,2).^2 - 4) < 1e-5);
interfaceFix = p2(interface_indx,:);

end



%%
% %------ mesh17
%-------------------- Combined StokesDarcy Nonconforming Quad Mesh ---------------------
function mesh18CombinedSDmesh_nonconformingQuad()
close all
clearvars
clc

for nn = 1:4
    close all;
    meshN = num2str(2^(nn+1));
    SDmeshName = ['SDmesh_nonconQuad_S[0,1][0,1]_D[0,1][-1,0]_',meshN];
    
    h = 1/(2^(nn+1));
    %-
    D_G_01 = cartGrid([1/(2*h), 1/(2*h)],[1/2,1/2]);
    D_G_02 = cartGrid([1/(h), 1/(h)],[1/2,1/2]);
    %-
    D_G1 = D_G_01;
    D_G2 = translateGrid(D_G_02,[1/2,0]);
    D_temp1 = glue2DGrid(D_G1, D_G2);
    %-
    D_G3 = translateGrid(D_G_01,[1/2,1/2]);
    D_G4 = translateGrid(D_G_02,[0,1/2]);
    D_temp2 = glue2DGrid(D_G3, D_G4);
    %-
    D_G = glue2DGrid(D_temp1, D_temp2);
    %plotGrid(D_G, 'FaceColor', 'none'); axis equal; axis off
    S_G = translateGrid(D_G,[0,1]);
    %-
    D_G.cells.tag = ones(D_G.cells.num,1);
    S_G.cells.tag = 2*ones(S_G.cells.num,1);
    SD_G = glue2DGrid(D_G, S_G);
    SD_G = translateGrid(SD_G,[0,-1]);
    theCells = SD_G.cells; 
    if isfield(theCells, 'indexMap'), theCells = rmfield(theCells, 'indexMap'); end
    SD_G.cells = theCells;
    %-
    allCells = (1:SD_G.cells.num)';
    subSD_D_cells = allCells(SD_G.cells.tag==1);
    subSD_S_cells = allCells(SD_G.cells.tag==2);
    %-
    subSD_D_G = extractSubgrid(SD_G,subSD_D_cells); %%--- extract subgrid
    subSD_S_G = extractSubgrid(SD_G,subSD_S_cells); %%--- extract subgrid
    %-
    %--- Stokes Domain [0,1]x[1,2]
    [Snode,Selem] = mrstG_2_myMeshInfo(subSD_S_G);
    node = Snode; elem = Selem;
    SmeshInfo = polyMeshAuxStructure(node, elem);
    
    %plotPolyMsh(SmeshInfo)
    %--- Darcy Domain [0,1]x[0,1]
    [Dnode,Delem] = mrstG_2_myMeshInfo(subSD_D_G);
    node = Dnode; elem = Delem;
    DmeshInfo = polyMeshAuxStructure(node, elem);
    
    %--- Combined StokesDarcy mesh
    [SDnode, SDelem] = mrstG_2_myMeshInfo(SD_G);
    node = SDnode; elem = SDelem;
    Selems = subSD_S_cells;
    Delems = subSD_D_cells;
    save(SDmeshName, 'node', 'elem', 'Selems', 'Delems')
 
    %--- plot
    PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
        DmeshInfo.node,DmeshInfo.elem,length(Delem))
    figure
    PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
        DmeshInfo.node,DmeshInfo.elem,length(Delem))
    
    %--- plot interface
    interfaceX = (0:0.1:1)';
    interfaceY = ones(length(interfaceX),1);
    plot(interfaceX,interfaceY,'-r','LineWidth',1.5);
    
end % nn
end % function getNonconformingQuadMesh()



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
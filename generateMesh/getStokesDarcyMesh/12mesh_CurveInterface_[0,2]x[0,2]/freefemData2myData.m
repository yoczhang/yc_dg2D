function freefemData2myData()
close all
clear
clc

addpath('../utils');
for n = 2:2
    addpath('../../')
    disp(['begin: n=',num2str(n)])
    close all
    folderpath = ['./SD_Mesh/0.7/',num2str(n)];
    addpath(folderpath);
    
    meshN = num2str(2^n);
    SmeshName = ['Smesh_curveInterface_[0,2]x[1,2]_',meshN];
    DmeshName = ['Dmesh_curveInterface_[0,2]x[0,1]_',meshN];
    finalFigName = ['SD_curveInterface_(0,1)x(1,2)_',meshN,'_backup'];
    
    load('my_elem_Darcy.dat')
    load('my_elem_Stokes.dat')
    load('my_nodes_Darcy.dat')
    load('my_nodes_Stokes.dat')
    rmpath(folderpath);
    
    node = my_nodes_Darcy(:,2:3);
    elem = my_elem_Darcy(:,end);
    SNelems = length(elem)/3;
    elem = reshape(elem,3,SNelems);
    elem = elem' + 1;
    Snode = node; Selem = elem;
    Selem_cell = cell(SNelems,1);
    for sn = 1:SNelems
        Selem_cell{sn} = elem(sn,:);
    end
    save(SmeshName,'node','elem');
    %patchPlotMesh(node,elem); hold on;
    
    node = my_nodes_Stokes(:,2:3); 
    elem = my_elem_Stokes(:,end);
    DNelems = length(elem)/3;
    elem = reshape(elem,3,DNelems);
    elem = elem' + 1;
    Dnode = node; Delem = elem;
    Delem_cell = cell(DNelems,1);
    for dn = 1:DNelems
        Delem_cell{dn} = elem(dn,:);
    end
% %     save(DmeshName,'node','elem');
% %     %patchPlotMesh(node,elem); hold off;
    
    
    %------- mrst --------
    G_S = pebi(triangleGrid(Snode,Selem));
    %G_S = triangleGrid(Snode,Selem);
    G_S = removeShortEdges(G_S,5e-4);
    G_S.nodes.coords(abs(G_S.nodes.coords(:,1)-0)<5e-7,1 ) = 0;
    G_S.nodes.coords(abs(G_S.nodes.coords(:,1)-2)<5e-7,1 ) = 2;
    G_S.nodes.coords(abs(G_S.nodes.coords(:,2)-2)<5e-7,2 ) = 2;
    plotGrid(G_S,'FaceColor','none'); hold on;
    G_D = pebi(triangleGrid(Dnode,Delem));
    %G_D = triangleGrid(Dnode,Delem);
    G_D = removeShortEdges(G_D,5e-4);
    G_D.nodes.coords(abs(G_D.nodes.coords(:,1)-0)<5e-7,1 ) = 0;
    G_D.nodes.coords(abs(G_D.nodes.coords(:,1)-2)<5e-7,1 ) = 2;
    G_D.nodes.coords(abs(G_D.nodes.coords(:,2)-0)<5e-7,2 ) = 0;
    
    xx0 = [0, 2, 2, 0, 0];
    yy0 = [0, 0, 2, 2,0];
    plot(xx0,yy0,'-k','LineWidth',1.5);
    
    finalFig = plotGrid(G_D,'FaceColor','none');
    xx1 = 0:0.05:2;
    yy1 = -sin(xx1*(pi/2)) .* cos(xx1*pi) * 0.15 + 1;
    plot(xx1,yy1,'-r','LineWidth',1.3);
    set(gca,'xlim',[0 2],'ylim',[0 2]);
    axis equal;
    axis off;
% %     saveas(finalFig,finalFigName,'fig')
    %-----------------------
    
    %------- glue2DGrid -------
    if n<=4
        Npartition = 10*2^(n-2)/2;
    else
        Npartition = 10*2^(n-3)/2;
    end
    
    G_D1 = G_D;
    G_D1.cells.tag = 1*ones(G_D1.cells.num,1);
    plotGrid(G_D1,'FaceColor','none');
    
    %--- G_D2 and G_D4
    G_D2 = cartGrid([ceil(0.8*Npartition*2),ceil(0.25*Npartition*2)],[0.8,0.25]);
    G_D2 = twister(G_D2);
    G_D2 = translateGrid(G_D2,[0,0.45]);
    G_D2.cells.tag = 2*ones(G_D2.cells.num,1);
    G_D4 = cartGrid([ceil(0.8*Npartition*2),ceil(0.25*Npartition*2)],[0.8,0.25]);
    G_D4 = twister(G_D4);
    G_D4 = translateGrid(G_D4,[1.2,0.2]);
    G_D4.cells.tag = 4*ones(G_D4.cells.num,1);
    plotGrid(G_D2,'FaceColor','none');
    plotGrid(G_D4,'FaceColor','none');
    
    %--- G_D3
    left = 0.8; right = 1.2; bottom = 0.2; top = 0.7;
    [p,t] = generate_Tri_P_T(left,right,bottom,top,[(right-left)/Npartition,(top-bottom)/Npartition]);
    G_D3_1 = triangleGrid(p,t);
    G_D3_1.cells.tag = 3*ones(G_D3_1.cells.num,1);
    plotGrid(G_D3_1,'FaceColor','none');
    
    left = 0; right = 0.8; bottom = 0.2; top = 0.45;
    [p,t] = generate_Tri_P_T(left,right,bottom,top,[(right-left)/Npartition/2,(top-bottom)/Npartition*2]);
    G_D3_2 = triangleGrid(p,t);
    G_D3_2.cells.tag = 3*ones(G_D3_2.cells.num,1);
    plotGrid(G_D3_2,'FaceColor','none');
    
    left = 1.2; right = 2; bottom = 0.45; top = 0.7;
    [p,t] = generate_Tri_P_T(left,right,bottom,top,[(right-left)/Npartition/2,(top-bottom)/Npartition*2]);
    G_D3_3 = triangleGrid(p,t);
    G_D3_3.cells.tag = 3*ones(G_D3_3.cells.num,1);
    plotGrid(G_D3_3,'FaceColor','none');
 
    %--- G_D5
    nx = ceil(2*Npartition*2); ny = ceil(0.2*Npartition*2); 
    ex = 2; ey = 0.2;
    G_D5 = cartGrid([nx,ny],[ex,ey]);
    c = G_D5.nodes.coords;
    I = or( or(any(c==0,2), any(c(:,1)==ex,2) ), any(c(:,2)==ey,2) );
    G_D5.nodes.coords(~I,:) = c(~I,:) + ey/(2^(n-1)*Npartition)*rand(sum(~I),2)-ex/(4^(n+1)*Npartition);
    G_D5.cells.tag = 5*ones(G_D5.cells.num,1);
    plotGrid(G_D5,'FaceColor','none');
    
    %--- glue
    G_D1 = addCellFacesTags(G_D1);
    G_D2 = addCellFacesTags(G_D2);
    G_D3_1 = addCellFacesTags(G_D3_1);
    G_D3_2 = addCellFacesTags(G_D3_2);
    G_D3_3 = addCellFacesTags(G_D3_3);
    G_D4 = addCellFacesTags(G_D4);
    G_D5 = addCellFacesTags(G_D5);

    GD_group1 = glue2DGrid(G_D2,G_D3_2);
    %figure; plotGrid(GD_group1,'FaceColor','none');
    
    GD_group1 = glue2DGrid(G_D3_1,GD_group1);
    %figure; plotGrid(GD_group1,'FaceColor','none');
    
    GD_group2 = glue2DGrid(G_D4,G_D3_3);
    %figure; plotGrid(GD_group2,'FaceColor','none');
    
    GD_group1 = glue2DGrid(GD_group1,GD_group2);
    %figure; plotGrid(GD_group1,'FaceColor','none');
    
    GD_group1 = glue2DGrid(GD_group1,G_D5);
    %figure; plotGrid(GD_group1,'FaceColor','none');
    
    G_D = glue2DGrid(GD_group1,G_D1);
    %G_D = removeShortEdges(G_D,5e-4);
    %figure; plotGrid(G_D,'FaceColor','none');
    %-------------------------------
    [node,elem] = mrstG_2_myMeshInfo(G_D);
	%DmeshInfo = polyMeshAuxStructure(node, elem);
    NcellsInd = (1:G_D.cells.num)';
    region1_cellsInd = NcellsInd(G_D.cells.tag==1);
    region2_cellsInd = NcellsInd(G_D.cells.tag==2);
    region3_cellsInd = NcellsInd(G_D.cells.tag==3);
    region4_cellsInd = NcellsInd(G_D.cells.tag==4);
    region5_cellsInd = NcellsInd(G_D.cells.tag==5);
% %     save(DmeshName,'node','elem','region1_cellsInd','region2_cellsInd',...
% %         'region3_cellsInd','region4_cellsInd','region5_cellsInd');
    
    %patchPlotMesh(DmeshInfo.node,DmeshInfo.elem);
    
    finalFig = figure; hold on;
    PolyMshr_PlotMsh_twoDomain2(Snode,Selem_cell,SNelems,...
        Dnode,Delem_cell,DNelems)
    xx1 = 0:0.05:2;
    yy1 = sin(xx1*(pi/2)) .* cos(xx1*pi) * 0.15 + 1;
    plot(xx1,yy1,'-r');
    axis equal;
    axis off;
% %     saveas(finalFig,finalFigName,'fig')
    

    
    disp(['complete: n=',num2str(n)])
end % for
end %function

%--- sub function1
function G = addCellFacesTags(G)
    G = computeGeometry(G);
    
    % Set orientation of faces
    hf = G.cells.faces(:,1);
    hf2cn = gridCellNo(G);
    sgn = 2*(hf2cn == G.faces.neighbors(hf, 1)) - 1;
    N   = bsxfun(@times, sgn, G.faces.normals(hf,:));
    N   = bsxfun(@rdivide, N, G.faces.areas(hf,:));
    n   = zeros(numel(hf),2); n(:,1)=1;

    % Add cell tags
    G.cells.faces(:,2) = zeros(size(hf));
    i = sum(N.*n,2)==-1; G.cells.faces(i,2) = 1;
    i = sum(N.*n,2)== 1; G.cells.faces(i,2) = 2;
    n = n(:,[2 1]);
    i = sum(N.*n,2)==-1; G.cells.faces(i,2) = 3;
    i = sum(N.*n,2)== 1; G.cells.faces(i,2) = 4;
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

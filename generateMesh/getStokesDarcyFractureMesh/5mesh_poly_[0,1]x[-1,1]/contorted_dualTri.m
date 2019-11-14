%------ mesh5
%-------------------- Poly Mesh ---------------------
function contorted_dualTri()
clearvars
clc
close all

for nn = 1:4
    close all;
    meshN = num2str(2^(nn+1));
    meshName = ['mesh_dualTri_[0,1]x[-1,1]_',meshN];
    
    initH = 1/(8*2^(nn-1));

    interfaceFix_x = (0:initH:1)';
    interfaceFix_y = 0*ones(size(interfaceFix_x));
    
    interfaceFix1 = [interfaceFix_x, interfaceFix_y]; % the Stokes-Darcy interface
    
    interfaceFix2 = [interfaceFix_x, -1/2+interfaceFix_y]; % the Darcy-Fracture interface
    interfaceFix2 = [interfaceFix1;interfaceFix2]; 
    %> there are two interface, 1). Stokes-Darcy interface. 2). the Darcy-Fracture interface.

    interfaceFix3 = [interfaceFix_x, -1/2+interfaceFix_y];
    
    %% Stokes Domain [0,1]x[0,1]
    squarDomain1 = [0,0;1,1];
    fd1 = @(p) drectangle(p,0,1,0,1);
    fh1 = @(p) 0.005*ones(size(p,1),1);
    pfix1 = [0,0;1,0;1,1;0,1;interfaceFix1];
    [p1,t1] = distmesh( fd1, fh1, initH, squarDomain1, pfix1 );
    %- repaire the boundary nodes
    x0ind = (abs(p1(:,1)-0)<1e-6);
    p1(x0ind,1) = 0;
    y0ind = (abs(p1(:,2)-0)<1e-6);
    p1(y0ind,2) = 0;
    x1ind = (abs(p1(:,1)-1)<1e-6);
    p1(x1ind,1) = 1;
    y1ind = (abs(p1(:,2)-1)<1e-6);
    p1(y1ind,2) = 1;
    
    [Snode, Selem] = dualmesh(p1,t1);
    %--- contorted the mesh
    c_theta = 0.075;
    xx = Snode(:,1);
    yy = Snode(:,2);
    Snode(:,1) = xx + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
    Snode(:,2) = yy + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
    
    SmeshInfo = polyMeshAuxStructure(Snode, Selem);
    
    %% Darcy Domain 
    % Domain-part 1: [0,1]x[-1/2,0]
    squarDomain2 = [0,-1/2;1,0];
    fd2 = @(p) drectangle(p,0,1,-1/2,0);
    fh2 = @(p) 0.002*ones(size(p,1),1);
    pfix2 = [0,0;0,-1/2;1,-1/2;1,0;interfaceFix2];
    [p2,t2] = distmesh( fd2, fh2, initH, squarDomain2, pfix2 );
    %- repaire the boundary nodes
    x0ind = (abs(p2(:,1)-0)<1e-6);
    p2(x0ind,1) = 0;
    y0ind = (abs(p2(:,2)-0)<1e-6);
    p2(y0ind,2) = 0;
    x1ind = (abs(p2(:,1)-1)<1e-6);
    p2(x1ind,1) = 1;
    y1ind = (abs(p2(:,2)+1/2)<1e-6);
    p2(y1ind,2) = -1/2;
    
    [Dnode, Delem] = dualmesh(p2,t2);
    
    %--- contorted the mesh
    c_theta = 0.075;
    xx = Dnode(:,1);
    yy = Dnode(:,2);
    Dnode(:,1) = xx + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
    Dnode(:,2) = yy + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
    
    DmeshInfo1 = polyMeshAuxStructure(Dnode, Delem);
    
    % Domain-part 2: [0,1]x[-1,-1/2]
    squarDomain3 = [0,-1;1,-1/2];
    fd3 = @(p) drectangle(p,0,1,-1,-1/2);
    fh3 = @(p) 0.002*ones(size(p,1),1);
    pfix3 = [0,-1/2;0,-1;1,-1;1,-1/2;interfaceFix3];
    [p3,t3] = distmesh( fd3, fh3, initH, squarDomain3, pfix3 );
    %- repaire the boundary nodes
    x0ind = (abs(p3(:,1)-0)<1e-6);
    p3(x0ind,1) = 0;
    y0ind = (abs(p3(:,2)+1/2)<1e-6);
    p3(y0ind,2) = -1/2;
    x1ind = (abs(p3(:,1)-1)<1e-6);
    p3(x1ind,1) = 1;
    y1ind = (abs(p3(:,2)+1)<1e-6);
    p3(y1ind,2) = -1;
    
    [Dnode, Delem] = dualmesh(p3,t3);
    %--- contorted the mesh
    c_theta = 0.075;
    xx = Dnode(:,1);
    yy = Dnode(:,2);
    Dnode(:,1) = xx + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
    Dnode(:,2) = yy + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
    
    DmeshInfo2 = polyMeshAuxStructure(Dnode, Delem);
    
    
    %% save the mesh-information
    %save(meshName, 'SmeshInfo', 'DmeshInfo')
   
    %% plot
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
%         DmeshInfo.node,DmeshInfo.elem,length(Delem))
%     figure
    PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
        DmeshInfo1.node,DmeshInfo1.elem,length(Delem))
    
    hold on
    DNelem = DmeshInfo2.Nelems;
    Delem = DmeshInfo2.elem;
    Delem = Delem(1:DNelem)';                 %Only plot the first block
    Dnode = DmeshInfo2.node;
    MaxNVer = max(cellfun(@numel,Delem));      %Max. num. of vertices in mesh
    PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
    ElemMat = cellfun(PadWNaN,Delem,'UniformOutput',false);
    ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
    patch('Faces',ElemMat,'Vertices',Dnode,'FaceColor','w'); pause(1e-6)
    
    plot([0,1],[0,0],'-r', 'LineWidth',1.2)
    plot([0,1],[-1/2,-1/2],'-r', 'LineWidth',1.2)
    axis tight
    axis off
    axis equal
    
    if nn == 2
        print('dualTri.jpg','-djpeg', '-r750')
        print('dualTri.eps','-depsc')
    end 
    
end % nn

end % function




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


%
%------ mesh4
%-------------------- Contorted Dual-Tri Mesh ---------------------
function getMassConservationMesh2()
for nn = 1:4
    close all;
    meshN = num2str(2^(nn+1));
    SmeshName = ['Smesh_2massConservation_[0,2]x[1,2]_',meshN];
    DmeshName = ['Dmesh_2massConservation_[0,2]x[0,1]_',meshN];
    finalFigName = ['SD_2massConservation_(0,1)x(1,2)_',meshN,'_backup'];
    %--- Stokes Domain [0,1]x[1,2]
    %[nodeOld,elemOld] = squarequadmesh([0,1,1,2],1/2^(nn+1)); 
    [nodeOld, elemOld] = generate_Tri_P_T(0,2,1,2,[1/2,1/2].^(nn+1));
    [Snode, Selem] = dualmesh(nodeOld,elemOld);
    
    %--- contorted the mesh
    c_theta = 0.075;
    xx = Snode(:,1);
    yy = Snode(:,2);
    Snode(:,1) = xx + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
    Snode(:,2) = yy + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
    
    SmeshInfo = polyMeshAuxStructure(Snode, Selem);
    node = SmeshInfo.node; elem = SmeshInfo.elem;
    save(SmeshName, 'node', 'elem')
 
    %--- Darcy Domain [0,1]x[0,1]
    %[nodeOld,elemOld] = squarequadmesh([0,1,0,1],1/2^(nn+1)); 
    [nodeOld, elemOld] = generate_Tri_P_T(0,2,0,1,[1/2,1/2].^(nn+1));
    [Dnode, Delem] = dualmesh(nodeOld,elemOld);
    %--- contorted the mesh
    c_theta = 0.075;
    xx = Dnode(:,1);
    yy = Dnode(:,2);
    Dnode(:,1) = xx + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
    Dnode(:,2) = yy + c_theta*sin(2*pi*xx).*sin(2*pi*yy);
    
    DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
    node = DmeshInfo.node; elem = DmeshInfo.elem;
    save(DmeshName, 'node', 'elem')
   
    PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
        DmeshInfo.node,DmeshInfo.elem,length(Delem))
    figure
    PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
        DmeshInfo.node,DmeshInfo.elem,length(Delem))
end % nn

end % function function getMassConservationMesh2()



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

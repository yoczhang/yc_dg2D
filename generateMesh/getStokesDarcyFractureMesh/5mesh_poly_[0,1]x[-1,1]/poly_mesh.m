%------ mesh 9
%-------------------- PolyMesher-symmetryPolygon Mesh ---------------------

function poly_mesh()

for nn = 2:4
    close all;
    meshN = num2str(2^(nn+1));
    meshName = ['poly_faults_[0,1]x[-1,1]_',meshN];
    
    %% Stokes mesh
    [Snode,Selem]=PolyMesher_toSDmesh(1/2^(nn+1),@poly_StokesDomain,4^(nn+1),1699);
    
    left_bound = 0; right_bound = 1;
    bottom_bound = 0; top_bound = 1;
    Snode = modifiedNode(Snode,left_bound,right_bound,bottom_bound,top_bound);
    SmeshInfo = polyMeshAuxStructure(Snode, Selem);
    
    %% Darcy mesh
    %--- Darcy domain1
    [Dnode,Delem]=PolyMesher_toSDmesh(1/2^(nn+1),@poly_DarcyDomain1,4^(nn+1)/2,800);
    
    left_bound = 0; right_bound = 1;
    bottom_bound = -1/2; top_bound = 0;
    Dnode = modifiedNode(Dnode,left_bound,right_bound,bottom_bound,top_bound);
    DmeshInfo1 = polyMeshAuxStructure(Dnode, Delem);
    
    %--- Darcy domain2
%     [Dnode,Delem]=PolyMesher_toSDmesh(1/2^(nn+1),@poly_DarcyDomain2,4^(nn+1)/2,800);
%     
%     left_bound = 0; right_bound = 1;
%     bottom_bound = -1; top_bound = -1/2;
%     Dnode = modifiedNode(Dnode,left_bound,right_bound,bottom_bound,top_bound);
%     DmeshInfo2 = polyMeshAuxStructure(Dnode, Delem);
    Dnode2 = Dnode;
    Dnode2(:,2) = -1 - Dnode2(:,2);
    Delem2 = Delem;
    DmeshInfo2 = polyMeshAuxStructure(Dnode2, Delem2);
    
    
    
    %% plot
%     PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,SmeshInfo.Nelems,...
%         DmeshInfo1.node,DmeshInfo1.elem,DmeshInfo1.Nelems)
%     figure
    PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,SmeshInfo.Nelems,...
        DmeshInfo1.node,DmeshInfo1.elem,DmeshInfo1.Nelems)
    
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
        print('polymesh.jpg','-djpeg', '-r750')
        print('polymesh.eps','-depsc')
    end 
end


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


function [x] = poly_StokesDomain(Demand,Arg)
    BdBox = [0 1 0 1];
    switch(Demand)
        case('Dist');  x = DistFnc(Arg,BdBox);
        case('BC');    x = BndryCnds(Arg{:},BdBox);
        case('BdBox'); x = BdBox;
        case('PFix');  x = FixedPoints(Arg,BdBox);
    end
end
function [x] = poly_DarcyDomain1(Demand,Arg)
    BdBox = [0 1 -1/2 0];
    switch(Demand)
        case('Dist');  x = DistFnc(Arg,BdBox);
        case('BC');    x = BndryCnds(Arg{:},BdBox);
        case('BdBox'); x = BdBox;
        case('PFix');  x = FixedPoints(Arg,BdBox);
    end
end
function [x] = poly_DarcyDomain2(Demand,Arg)
    BdBox = [0 1 -1 -1/2];
    switch(Demand)
        case('Dist');  x = DistFnc(Arg,BdBox);
        case('BC');    x = BndryCnds(Arg{:},BdBox);
        case('BdBox'); x = BdBox;
        case('PFix');  x = FixedPoints(Arg,BdBox);
    end
end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
end
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
    x = cell(2,1); % No boundary conditions specified for this problem
end
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(h,BdBox)
%   xx = (BdBox(1):h:BdBox(1))';
%   yy = BdBox(3) * ones(length(xx),1);
%   PFix = [xx, yy];
  PFix = [];
end
%-------------------------------------------------------------------------%


function node = modifiedNode(node,left_bound,right_bound,bottom_bound,top_bound)

    xx = node(:,1);
    yy = node(:,2);

    tol = 6e-9;

    left_abs_xx = abs( xx - left_bound );
    right_abs_xx = abs( xx - right_bound );
    bottom_abs_yy = abs( yy - bottom_bound );
    top_abs_yy = abs( yy - top_bound );

    node(left_abs_xx<tol,1) = left_bound;
    node(right_abs_xx<tol,1) = right_bound;
    node(bottom_abs_yy<tol,2) = bottom_bound;
    node(top_abs_yy<tol,2) = top_bound;

end % function modifiedNode
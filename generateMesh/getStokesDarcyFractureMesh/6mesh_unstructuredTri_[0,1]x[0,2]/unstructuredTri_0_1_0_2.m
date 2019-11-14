function unstructuredTri_0_1_0_2()
%
%
%

clearvars;
clc;
close all

for nn = 1:4
    meshN = num2str(2^(nn+1));
    meshName = ['unstructuredTri_[0,1]x[0,2]_',meshN];
    
    initH = 1/(5*2^(nn-1));

    interfaceFix_x = (0:initH:1)';
    interfaceFix_y = 0*ones(size(interfaceFix_x));
    interfaceFix1 = [interfaceFix_x, interfaceFix_y]; % the Stokes-Darcy interface
    
    interfaceFix2 = [interfaceFix_x, -1/2+interfaceFix_y]; % the Darcy-Fracture interface
    interfaceFix2 = [interfaceFix1;interfaceFix2]; 
    %> there are two interface, 1). Stokes-Darcy interface. 2). the Darcy-Fracture interface.


    %% Stokes domain
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
    
    p1(:,2) = p1(:,2) + 1;


    %% Darcy domain
    squarDomain2 = [0,-1;1,0];
    fd2 = @(p) drectangle(p,0,1,-1,0);
    fh2 = @(p) 0.002*ones(size(p,1),1);
    pfix2 = [0,0;0,-1;1,-1;1,0;interfaceFix2];
    [p2,t2] = distmesh( fd2, fh2, initH, squarDomain2, pfix2 );
    %- repaire the boundary nodes
    x0ind = (abs(p2(:,1)-0)<1e-6);
    p2(x0ind,1) = 0;
    y0ind = (abs(p2(:,2)-0)<1e-6);
    p2(y0ind,2) = 0;
    x1ind = (abs(p2(:,1)+1)<1e-6);
    p2(x1ind,1) = -1;
    y1ind = (abs(p2(:,2)-1)<1e-6);
    p2(y1ind,2) = 1;
    
     p2(:,2) = p2(:,2) + 1;

    SmeshInfo = polyMeshAuxStructure(p1, t1);
    DmeshInfo = polyMeshAuxStructure(p2, t2);


    %% find fracture
    edge0 = DmeshInfo.baryEdge;
    Nedge = size(edge0,1);
    edgeIndx = 1:Nedge;

    line1x = [0, 1]';
    line1y = 1/2*ones(size(line1x));
    P1A = [line1x(1), line1y(1)];
    P1B = [line1x(2), line1y(2)];

    dis1 = DisPtToLine(edge0,P1A,P1B);
    fracture1 = edgeIndx(abs(dis1)<1e-8);

    %----------
    splitedFractureFace = cell(1,1);
    splitedFractureFace{1} = fracture1;
    DmeshInfo.splitedFractureFace = splitedFractureFace;

    %% - save the the mat file
    save(meshName, 'SmeshInfo', 'DmeshInfo')

    %% plot
    doplot = true;
    if doplot
        Dnode = DmeshInfo.node;
        Dedge = DmeshInfo.edge;
        patch( 'vertices', p1, 'faces', t1, 'facecolor', [1, 1, 1] )
        hold on
        %plot(interfaceFix1(:,1),interfaceFix1(:,2),'or')

        patch( 'vertices', p2, 'faces', t2, 'facecolor', [1, 1, 1] )
        %plot(interfaceFix2(:,1),interfaceFix2(:,2),'or')

        plot([0,1],[1,1],'-r', 'LineWidth',1.2)

        for k = 1:length(fracture1)
            nodeInd = Dedge(fracture1(k),:);
            coord = Dnode(nodeInd,:);
            plot(coord(:,1),coord(:,2),'-r', 'LineWidth',1.2)
        end
        axis tight
        axis off
        axis equal
        if nn == 2
            print('unstructuredTri.jpg','-djpeg', '-r750')
            print('unstructuredTri.eps','-depsc')
        end
    end

end % for nn



end % function
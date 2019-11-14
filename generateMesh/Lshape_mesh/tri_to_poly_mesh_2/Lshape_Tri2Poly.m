function Lshape_Tri2Poly()
%
%
%	YcZhang 1/8/2018
%
%   Last modified 1/8/2018
%
clear
close all
clc

loadmesh_index = [4, 8, 16, 32, 64];
savemesh_index = [4, 8, 16, 32, 64];

for k = 1:length(loadmesh_index)
    close all
    loadfilename = ['./Lshape_tri_',num2str(loadmesh_index(k)),'.mat'];
    savefilename = ['./Lshape_poly_',num2str(savemesh_index(k)),'.mat'];
    load(loadfilename)
    %------- mrst --------
    G = pebi(triangleGrid(node,elem));
    G = removeShortEdges(G,5e-6);
%     G.nodes.coords(abs(G.nodes.coords(:,1)-0)<5e-7,1 ) = 0;
%     G.nodes.coords(abs(G.nodes.coords(:,1)-2)<5e-7,1 ) = 2;
%     G.nodes.coords(abs(G.nodes.coords(:,2)-2)<5e-7,2 ) = 2;
    
    plotGrid(G,'FaceColor','none'); hold on;
    
    [node,elem] = mrstG_2_myMeshInfo(G);
    save(savefilename,'node','elem')
    
end

end % function
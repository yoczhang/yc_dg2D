function test_plotLevelSetFunc
%
%
%   YcZhang 23/9/2017
%
%   Last modified 23/9/2017
%

clc
clearvars
close all

test_case = 3;

switch test_case
    case 3
        % get node and elem
        [node, elem] = generate_Tri_P_T(0,1,0,1,[1/4,1/4]);
        %patchPlotMesh(node, elem);
        meshInfo = polyMeshAuxStructure(node,elem);
        plotPolyMsh(meshInfo)
        
        % get the level-set func1
        aa = 1; bb = 0; cc = -0.5;
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];
        
        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        cutElem = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem)
        
        % get the level-set func2
        aa = 1; bb = 0; cc = -1.5;
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0,0.7];
        
        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        cutElem = getEnrichedElem(levelFunc, node, elem);
        patchPlotEnrichment(node, elem, cutElem)
        
        % get the level-set func3
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-1).^2 + (yy-0.6).^2 );
        ff = -(Dist - 0.1);
        x_domain=[0,1]; y_domain=[0,1];
        
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        cutElem = getEnrichedElem(levelFunc, node, elem);
        patchPlotEnrichment(node, elem, cutElem)
        
        % get the level-set func4
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-1).^2 + (yy-0.4).^2 );
        ff = -(Dist - 0.1);
        x_domain=[1,2]; y_domain=[0,1];
        
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        cutElem = getEnrichedElem(levelFunc, node, elem);
        patchPlotEnrichment(node, elem, cutElem)
        
        
    case 1
        % Get mesh.
        nElem = 4;
        [NodeNum, ElemNum, xx, yy, Mesh] = GetMesh(1, 1, nElem, nElem);
        xx = xx*2-1; yy = yy*2-1;


        yc_node = [xx,yy];
        yc_elem = Mesh;
        yc_meshInfo = polyMeshAuxStructure(yc_node, yc_elem);
        plotPolyMsh(yc_meshInfo)

        figure(2)
        PlotMesh_yctest(yc_node, yc_elem)
        cutElem = divideMesh2Submesh_test(yc_node, yc_elem);
        
        PlotEnrichment_yctest(yc_elem, yc_node, cutElem)
        
    case 2
        % Get mesh.
        nElem = 4;
        [NodeNum, ElemNum, xx, yy, Mesh] = GetMesh(1, 1, nElem, nElem);
        xx = xx*2-1; yy = yy*2-1;
        PlotMesh(Mesh, xx, yy, ElemNum)
        
        % Get level-set function.
        aa=0; bb=1; cc=-0.5; % Vertical disc. at x=0.5.
        ff = GetLevelSet(aa, bb, cc, xx, yy);
        PlotLevelSet(Mesh, xx, yy, ff, ElemNum)
        
        % Get enriched elements and nodes.
        [ElemsEnriched, NodesEnriched] = GetEnrichedNodesElems(Mesh, ff);
        PlotEnrichment(Mesh, xx, yy, ElemsEnriched, NodesEnriched)
        
end % switch


end % function
function enrichedElem = getEnrichedElemForExamples(node, elem, testcase)
%
%   Give some cases to get the enriched elems.
%
%   input: 
%       node, matrix.
%       elem, cell-type or matrix.
%       testcase, 1,2,3, or 4.
%
%       
%
%   YcZhang 2/10/2017
%
%   Last modified 2/10/2017
%
%

Nnodes = size(node,1);
patchPlotMesh(node, elem);
switch testcase
    case 1
        %% get sub-enriched elem and plot
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

        cutElem1 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1)

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

        cutElem2 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem2)

        % get the level-set func3
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-1).^2 + (yy-0.6).^2 );
        ff = -(Dist - 0.1);
        x_domain=[0,1]; y_domain=[0,1];

        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem3 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem3)

        % get the level-set func4
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-1).^2 + (yy-0.4).^2 );
        ff = -(Dist - 0.1);
        x_domain=[1,2]; y_domain=[0,1];

        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem4 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem4)
        
        %% get the all-enriched elem
        enrichedElem = unique([cutElem1; cutElem2; cutElem3; cutElem4]);
        
    case 2
        mesh_h = 1/40;
        %% get sub-enriched elem and plot
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

        
        [cutElem1,~] = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1)

        % get the level-set func2
        aa = 1; bb = 0; cc = -1.5;
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem2 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem2)
        
        %% get the all-enriched elem
        enrichedElem = unique([cutElem1; cutElem2]);
        
    case 3
        h = 1./(sqrt(Nnodes)-1);
        % get the level-set func1
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-1).^2 + (yy-0.5).^2 );
        ff = -(Dist - 0.3);

        levelFunc.ff = ff;

        cutElem1 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1)
        
        % get the level-set func2
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-1).^2 + (yy-0.5).^2 );
        ff = -(Dist - (0.3+h));

        levelFunc.ff = ff;

        cutElem2 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem2)
        
        % get the level-set func2
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-1).^2 + (yy-0.5).^2 );
        ff = -(Dist - (0.3-h));

        levelFunc.ff = ff;

        cutElem3 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem3)
        
        
        %% get the all-enriched elem
        enrichedElem = unique([cutElem1; cutElem2; cutElem3]);
        
    case 4
        mesh_h = 1./(sqrt(Nnodes)-1);
        %% get sub-enriched elem and plot
        % get the level-set func1
        aa = 1; bb = 0; cc = -(0.5-mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        [cutElem1,~] = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem1)
        
        % get the level-set func1_1
        aa = 1; bb = 0; cc = -(0.5-2*mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        [cutElem1_1,~] = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem1_1)
        
        
        % get the level-set func1_2
        aa = 1; bb = 0; cc = -(0.5-3*mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
%         [cutElem1_2,~] = getEnrichedElem(levelFunc, node, elem, 0);
%         patchPlotEnrichment(node, elem, cutElem1_2)
        cutElem1_2 = [];
        
        % get the level-set func2
        aa = 1; bb = 0; cc = -(0.5+mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        [cutElem2,~] = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem2)
        
        % get the level-set func2_1
        aa = 1; bb = 0; cc = -(0.5+2*mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        [cutElem2_1,~] = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem2_1)
        
        % get the level-set func2_2
        aa = 1; bb = 0; cc = -(0.5+3*mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        [cutElem2_2,~] = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem2_2)

        % get the level-set func3
        aa = 1; bb = 0; cc = -(1.5-mesh_h);
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem3 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem3)
        
        % get the level-set func4
        aa = 1; bb = 0; cc = -(1.5+mesh_h);
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem4 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem4)
        
        %% get the all-enriched elem
        enrichedElem = unique([cutElem1; cutElem1_1; cutElem1_2; cutElem2; cutElem2_1; cutElem2_2; cutElem3; cutElem4]);
        
    case 5
        h = 1./(sqrt(Nnodes)-1);
        
        %% get sub-enriched elem and plot
        % get the level-set func1
        aa = 1; bb = 0; cc = -1;
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem1 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1)
        
        % get the level-set func1_1
        aa = 1; bb = 0; cc = -(1+h);
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem1_1 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1_1)
        
         % get the level-set func1_2
        aa = 1; bb = 0; cc = -(1-h);
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem1_2 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1_2)

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

        cutElem2 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem2)
        
        % get the level-set func2_1
        aa = 1; bb = 0; cc = -(1.5+h);
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0,0.7];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem2_1 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem2_1)
        
        % get the level-set func2_2
        aa = 1; bb = 0; cc = -(1.5-h);
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0,0.7];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem2_2 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem2_2)

        % get the level-set func3
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-0.3).^2 + (yy-0.65).^2 );
        ff = -(Dist - 0.15);
        x_domain=[0,0.3]; y_domain=[0,1];

        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem3 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem3)
        
        % get the level-set func3_1
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-(0.3-h)).^2 + (yy-(0.65+h)).^2 );
        ff = -(Dist - 0.15);
        x_domain=[0,0.3+h]; y_domain=[0,1];

        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem3_1 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem3_1)

        % get the level-set func4
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-0.3).^2 + (yy-0.35).^2 );
        ff = -(Dist - 0.15);
        x_domain=[0.3,1]; y_domain=[0,1];

        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem4 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem4)
        
        % get the level-set func4_1
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-(0.3+h)).^2 + (yy-(0.35-h)).^2 );
        ff = -(Dist - 0.15);
        x_domain=[0.3,1]; y_domain=[0,1];

        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem4_1 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem4_1)
        
        %% get the all-enriched elem
        enrichedElem = unique([cutElem1; cutElem1_1; cutElem1_2; cutElem2; cutElem2_1; ...
            cutElem2_2; cutElem3; cutElem3_1; cutElem4; cutElem4_1]);
        
    case 6
        mesh_h = 1./(sqrt(Nnodes)-1);
        %% get sub-enriched elem and plot
        % get the level-set func1
        aa = 1; bb = 0; cc = -0.5;
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.203,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        
        [cutElem1,~] = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem1)

        % get the level-set func1
        aa = 1; bb = 0; cc = -(0.5-mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.203,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        
        [cutElem1_1,~] = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem1_1)
        
        % get the level-set func1
        aa = 1; bb = 0; cc = -(0.5+mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.203,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        
        [cutElem1_2,~] = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem1_2)
        
        % get the level-set func2
        aa = 1; bb = 0; cc = -1.5;
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.203,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem2 = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem2)
        
        % get the level-set func1
        aa = 1; bb = 0; cc = -(1.5-mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.203,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        
        [cutElem2_1,~] = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem2_1)
        
        % get the level-set func1
        aa = 1; bb = 0; cc = -(1.5+mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.203,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        
        [cutElem2_2,~] = getEnrichedElem(levelFunc, node, elem, 0);
        patchPlotEnrichment(node, elem, cutElem2_2)
        
        %% get the all-enriched elem
        enrichedElem = unique([cutElem1; cutElem1_1; cutElem1_2; cutElem2; cutElem2_1; cutElem2_2;]);
        
    case 7 
        % the random enrichment
        Nelems = size(elem,1);
        randNum = rand(Nelems,1); % generate the (0~1) random num.
        
        enrichedElem = (1:Nelems)';
        
        enrichedElem = enrichedElem(randNum>0.6);
        patchPlotEnrichment(node, elem, enrichedElem)
        
    case 8 % all domain
        Nelems = size(elem,1);
        enrichedElem = (1:Nelems)';
        patchPlotEnrichment(node, elem, enrichedElem)

end % switch

end % function

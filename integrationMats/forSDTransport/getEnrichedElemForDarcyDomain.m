function enrichedElem = getEnrichedElemForDarcyDomain(node, elem, testcase)
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
%   YcZhang 10/10/2017
%
%   Last modified 10/10/2017
%
%


%patchPlotMesh(node, elem);
switch testcase
    case 1 % darcy domain [0,1]x[0,1/2]
        
        Nnodes = size(node,1); h = 1./(sqrt(Nnodes)-1);
        
        %% get sub-enriched elem and plot
        
        % get the level-set func1
        % a = 0; b = 1; c =-0.5; => horizontal line at y = 0.5
        aa = 0; bb = 1; cc = -0.25;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,0.3]; y_domain=[0,1];
        
        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem1 = getEnrichedElem(levelFunc, node, elem, 1);
        
        % get the level-set func2
        % a = 0; b = 1; c =-0.5; => horizontal line at y = 0.5
        aa = 0; bb = 1; cc = -(0.25+h);
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,0.3]; y_domain=[0,1];
        
        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem2 = getEnrichedElem(levelFunc, node, elem, 1);
        
        % get the level-set func3
        % a = 0; b = 1; c =-0.5; => horizontal line at y = 0.5
        aa = 0; bb = 1; cc = -(0.25-h);
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,0.3]; y_domain=[0,1];
        
        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem3 = getEnrichedElem(levelFunc, node, elem, 1);
        
        %% get the all-enriched elem
        enrichedElem = unique([cutElem1; cutElem2; cutElem3]);
        patchPlotEnrichment(node, elem, enrichedElem)

    case 2 
        % the random enrichment
        Nelems = size(elem,1);
        randNum = rand(Nelems,1); % generate the (0~1) random num.
        
        enrichedElem = (1:Nelems)';
        
        enrichedElem = enrichedElem(randNum>0.6);
        patchPlotEnrichment(node, elem, enrichedElem)
        
    case 3 
        % all Darcy Domain is enriched elems
        Nelems = size(elem,1);
        enrichedElem = (1:Nelems)';

end % switch

end % function

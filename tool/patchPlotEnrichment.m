function patchPlotEnrichment(node, elem, ElemsEnriched)
%
%   %-----------------------------------------------------
%       Using the matlab func patch(..).
%       Just copy form PlotMesh.m (from SemXFEM_2d)
%   %-----------------------------------------------------
%
%
%
%   YcZhang 24/9/2017
%
%   Last modified 24/9/2017
%
%

% Plot the enriched elements and nodes.
xx = node(:,1);
yy = node(:,2);

for i = 1 : length(ElemsEnriched)
    CurrElem = ElemsEnriched(i);
    
    if iscell(elem)
        CurrNodes = elem{CurrElem};
    else
        CurrNodes = elem(CurrElem, :);
    end

	L = ones(1,length(CurrNodes));
    %patch(xx(CurrNodes), yy(CurrNodes), -0.004*L, 'r')
    patch(xx(CurrNodes), yy(CurrNodes), -0.004*L, [1, .5, .5])
        %> %> also may use: patch(x, y, [r g b]) to control the color.
end

end % function
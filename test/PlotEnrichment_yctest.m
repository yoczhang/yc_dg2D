% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function PlotEnrichment_yctest(elem, node, ElemsEnriched)

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
    patch(xx(CurrNodes), yy(CurrNodes), -0.004*L, 'r')
end

end % function
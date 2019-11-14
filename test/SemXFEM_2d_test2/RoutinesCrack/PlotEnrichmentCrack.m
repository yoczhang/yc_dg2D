% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function PlotEnrichmentCrack(Mesh, xx, yy, ElemsEnriched1, NodesEnriched1, ...
    ElemsEnriched2, NodesEnriched2)

% Plot the sign- and branch-enriched elements and nodes.

for i = 1 : size(ElemsEnriched1, 1)
    CurrElem = ElemsEnriched1(i);
    Nodes = Mesh(CurrElem, :);
    patch(xx(Nodes), yy(Nodes), -0.004*[1 1 1 1], 0.9*[1 1 1])
end

plot3(xx(NodesEnriched1), yy(NodesEnriched1), -0.001*ones(length(NodesEnriched1), 1), 'ko')

for i = 1 : size(ElemsEnriched2, 1)
    CurrElem = ElemsEnriched2(i);
    Nodes = Mesh(CurrElem, :);
    patch(xx(Nodes), yy(Nodes), -0.004*[1 1 1 1], 0.7*[1 1 1])
end

plot3(xx(NodesEnriched2), yy(NodesEnriched2), -0.001*ones(length(NodesEnriched2), 1), 'ks')

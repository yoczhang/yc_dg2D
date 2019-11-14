% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [ElemsEnriched, NodesEnriched] = GetEnrichedNodesElems(Mesh, ff)

% Get cut elements and correspoding nodes.

ElemsEnriched = find(min(sign(ff(Mesh))') ~= max(sign(ff(Mesh))'));
ElemsEnriched = ElemsEnriched';

NodesEnriched = Mesh(ElemsEnriched, :);
NodesEnriched = NodesEnriched(:);
NodesEnriched = unique(NodesEnriched);

function interfaceEdge = findinterfaceedge_new(node, sElem, tElem, DomainEdge, phi)
%
%
%
%
%
%	YcZhang 12/25/2018
%
%   Last modified 12/25/2018
%

NTS = size(sElem,1);
isExteriorSElem = false(NTS,1);
c = (node(sElem(:,1),:) + node(sElem(:,2),:) + node(sElem(:,3),:)+node(sElem(:,4),:))/4;
isExteriorSElem(phi(c)>0) = true;
exteriorSElem = sElem(isExteriorSElem,:);
exteriorSElemShift = circshift(exteriorSElem,[0,-1]);
exteriorSEdges = [exteriorSElem(:), exteriorSElemShift(:)];

NTT = size(tElem,1);
isExteriorTElem = false(NTT,1);
c = (node(tElem(:,1),:) + node(tElem(:,2),:) + node(tElem(:,3),:))/3;
isExteriorTElem(phi(c)>0) = true;
exteriorTElem = tElem(isExteriorTElem,:);
exteriorTElemShift = circshift(exteriorTElem,[0,-1]);
exteriorTEdges = [exteriorTElem(:), exteriorTElemShift(:)];

allexteriorEdges = [exteriorSEdges; exteriorTEdges];

allexteriorEdges = sort(allexteriorEdges,2);

[Ei, Ej, Es] = find( sparse(allexteriorEdges(:,2), allexteriorEdges(:,1), 1) );
exteriorbdFace = [Ej(Es==1), Ei(Es==1)];


DomainEdge = sort(DomainEdge,2);

interfaceEdge = setdiff(exteriorbdFace, DomainEdge,'rows','stable');



end % function 
%
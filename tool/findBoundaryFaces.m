function bdFace = findBoundaryFaces(elem)
%
%
%   %# ToDo: transform to C file.
%
%
%
%	YcZhang 10/6/2018
%
%   Last modified 10/6/2018
%

[r,c] = size(elem);
if r < c
    elem = elem'; % we need the [Nelem x Npoints_Of_One_element] form.
end % if r<c

if ~iscell(elem)
    if size(elem,2) < 3
        elem = elem';
    end
    N_elem = size(elem,1);
    elem_cell = cell(N_elem, 1);
    for i = 1:N_elem
        elem_cell{i} = elem(i,:);
    end % for i
    elem = elem_cell;
end % if

Nelems = size(elem, 1);

totalFaces = [];
for ii = 1:Nelems
    singleElem = elem{ii,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = length(singleElem); % get the number of edges on singleElem.
    
    %- get the all edges of each elem
    singleElem2E = zeros(singleNE,2);
    for jj = 1:singleNE-1
        % get the all edges of each elem
        singleElem2E(jj,:) = singleElem([jj,jj+1]);
    end % for jj
    singleElem2E(singleNE,:) = singleElem([singleNE,1]);
    singleElem2E = sort(singleElem2E,2);
    
    %- get the totalEdges 
    totalFaces_temp = totalFaces;
    totalFaces = [totalFaces_temp; singleElem2E];
end

%- bdEdge
[Ei, Ej, Es] = find( sparse(totalFaces(:,2), totalFaces(:,1), 1) );
bdFace = [Ej(Es==1), Ei(Es==1)];


end % function

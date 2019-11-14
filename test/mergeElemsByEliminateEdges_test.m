function mergeElemsByEliminateEdges_test
%
%
%   We let Nnodes denote the number of nodes of Th,
%               Nelems denote the number of the elements of Th.
%
%   input:
%       node, [Nnodes x 2], the (x-coord, y-coord) of nodes of Th.
%       interfaceFunc, the 
%
%
%   YcZhang 15/11/2017
%
%   Last modified 15/11/2017
%

clc;
clearvars;
close all

%--- mesh 1
% load('Delaunaymesh_80times40_[0_2]_[0_1]');

%--- mesh 2
%----------------- BR paper mesh ---------------------
n=3;
h_partition = 1/2^(n+1);
[node, elem] = get_BR_paper_mesh(0,1,0,1,h_partition);
%-----------------------------------------------------------

%--- mesh 3

%---------------------- polygon -------------------------
% load('Polygon_1')
% node = vertices'; elem = elements';
%-----------------------------------------------------------

Last_meshInfo = polyMeshAuxStructureForMergeElems(node, elem);
Nedges = Last_meshInfo.Nedges;

for nCiricle = 1:Nedges/6

    %%
    %--- Last meshInfo
    Last_meshInfo = polyMeshAuxStructureForMergeElems(node, elem);
    % plotPolyMsh(Last_meshInfo)
    patchPlotMeshForMergeElems(Last_meshInfo.node,Last_meshInfo.elem)


    %--- New meshInfo
    interEdgeIndex = Last_meshInfo.interEdgeIndex;
    randomEdgeIndex = randsrc(1,1,interEdgeIndex'); 
        %> get the random edge index form the interEdges.
    % randomEdgeIndex = 24;
    % randomEdgeIndex = 13;
    % randomEdgeIndex = 16;
    randomEdge = Last_meshInfo.edge(randomEdgeIndex,1:2);

    elem1_index = Last_meshInfo.edge2elem(randomEdgeIndex,1);
    elem2_index = Last_meshInfo.edge2elem(randomEdgeIndex,2);

    elemP_index = min(elem1_index,elem2_index); elemM_index = max(elem1_index,elem2_index);

    elemP = Last_meshInfo.elem{elemP_index};
    elemM = Last_meshInfo.elem{elemM_index};
    
    %--- control the max numEdges of elemP and elemM
    maxNumEdgeOfElem = 6;
    if length(elemP) > maxNumEdgeOfElem || length(elemM) > maxNumEdgeOfElem
        continue
        % if length(elemP) > maxNumEdgeOfElem or length(elemM) > maxNumEdgeOfElem,
        % we just skip the following programs
    end


    temp_P = 1:length(elemP);
    randomEdge_point1Locate_at_elemP = temp_P(elemP==randomEdge(1));
    randomEdge_point2Locate_at_elemP = temp_P(elemP==randomEdge(2));
    minLocate_at_elemP = min(randomEdge_point1Locate_at_elemP,randomEdge_point2Locate_at_elemP);
    maxLocate_at_elemP = max(randomEdge_point1Locate_at_elemP,randomEdge_point2Locate_at_elemP);

    temp_M = 1:length(elemM);
    randomEdge_point1Locate_at_elemM = temp_M(elemM==randomEdge(1));
    randomEdge_point2Locate_at_elemM = temp_M(elemM==randomEdge(2));
    minLocate_at_elemM = min(randomEdge_point1Locate_at_elemM,randomEdge_point2Locate_at_elemM);
    maxLocate_at_elemM = max(randomEdge_point1Locate_at_elemM,randomEdge_point2Locate_at_elemM);

    % Last_elemP = elemP;
    % length_P_1 = minLocate_at_elemP;
    % length_P_2 = length(elemP)-minLocate_at_elemP;
    % elemP(1:length_P_2+1) = Last_elemP(length_P_1:end);
    % elemP(length_P_2+2:length(elemP)) = Last_elemP(1:length_P_1-1);

    %--- rerange the elemM, to let the min nodeIndex at first.
    Last_elemM = elemM;
    length_M_1 = minLocate_at_elemM;
    length_M_2 = length(elemM)-minLocate_at_elemM;
    elemM(1:length_M_2+1) = Last_elemM(length_M_1:end);
    elemM(length_M_2+2:length(elemM)) = Last_elemM(1:length_M_1-1);

    elemP_remainNode_index = setdiff(elemP,randomEdge,'stable');
    elemM_remainNode_index = setdiff(elemM,randomEdge,'stable');

    if minLocate_at_elemP ==1 && maxLocate_at_elemP == length(elemP)
        Last_elemP_2 = elemP;
        elemP = [Last_elemP_2,elemM_remainNode_index];
    else
        Last_elemP_2 = elemP;
        elemP = [Last_elemP_2(1:minLocate_at_elemP),elemM_remainNode_index,Last_elemP_2(maxLocate_at_elemP:end)];
    end 
    elemP_unique = unique(elemP,'stable');

    if length(elemP_unique) ~= length(elemP)
        continue
        %elemP = elemP_unique;
    end 


    NewElem = cell(Last_meshInfo.Nelems-1,1);
    for CurrElem = 1:Last_meshInfo.Nelems

        if CurrElem < elemM_index && (CurrElem ~= elemP_index)
            NewElem{CurrElem} = Last_meshInfo.elem{CurrElem};
        end 
        if CurrElem == elemP_index
            NewElem{CurrElem} = elemP;
        end
        if CurrElem > elemM_index
            NewElem{CurrElem-1} = Last_meshInfo.elem{CurrElem};
        end 
    end

    %--- reset ---
    % node = node; 
    elem = NewElem;
    %--- plot the mesh
    patchPlotMeshForMergeElems(node,elem)
end % for nCiricle

%--------------------------- unused ----------------------------------
% elemP = zeros(1,length(elemP)+length(elemM_remainNode_index));
% 
% if minLocate_at_elemP == 1
%     elemP = [Last_elemP,elemM_remainNode_index];
% else
%     elemP(1:minLocate_at_elemP) = Last_elemP(1:minLocate_at_elemP);
%     elemP(minLocate_at_elemP+1:minLocate_at_elemP+length(elemM_remainNode_index)) = ...
%         elemM_remainNode_index;
%     elemP(minLocate_at_elemP+length(elemM_remainNode_index)+1:end) = elemP_remainNode_index;
% 
% end % if minLocate
%     
% if ismember(Last_elemP(1),randomEdge) || ismember(Last_elemP(2),randomEdge)
%     elemP = [Last_elemP,elemM_remainNode_index];
% 
% end % if ismember
%---------------------------------------------------------------------------

New_meshInfo = polyMeshAuxStructureForMergeElems(Last_meshInfo.node, NewElem);
patchPlotMeshForMergeElems(New_meshInfo.node,New_meshInfo.elem)
end

%% --- sub function ----
function g = polyMeshAuxStructureForMergeElems(node, elem)
%
%
%	YcZhang 15/11/2017
%
%   Last modified 15/11/2017
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

% we need the form of [Nnodes x 2]
[r, c] = size(node);
if r < c
    node = node';
end % if r<c

%%
Nnodes = size(node, 1);
Nelems = size(elem, 1);
elem2edge = cell(Nelems, 1);
baryElem = zeros(Nelems, 2);

%%
% we first need to get the number of allEdges TO  allocate memory of
% the following varibale 'allEdges'.
max_singleNE = 0; % to get the following: 15, edge2elem.
NtotalEdges = 0;
for ii = 1:Nelems
    singleElem = elem{ii,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges on singleElem.
    NtotalEdges = NtotalEdges + singleNE;
    
    if singleNE > max_singleNE
        max_singleNE = singleNE;
    end % if
end % for ii

%%
% to get the allEdges
totalEdges = zeros(NtotalEdges, 2);
NtotalEdges = 0; % here we need the NallEdges as the count index.
for ii = 1:Nelems
    singleElem = elem{ii,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges(nodes) on singleElem.
    
    %% get the barycenter of elem
    baryElem(ii,1) = sum( node(singleElem,1) )/singleNE;
    baryElem(ii,2) = sum( node(singleElem,2) )/singleNE;
    
    % get the all edges of each elem
    singleElem2E = zeros(singleNE,2);
    for jj = 1:singleNE-1
        % get the all edges of each elem
        singleElem2E(jj,:) = singleElem([jj,jj+1]);
    end % for jj
    singleElem2E(singleNE,:) = singleElem([singleNE,1]);
    %% rearrange the every row of singleElem2E in order from smallest to largest. 
    singleElem2E = sort(singleElem2E,2);
    
    %% get the totalEdges
    totalEdges(NtotalEdges+1:NtotalEdges+singleNE, :) = singleElem2E;
    NtotalEdges = NtotalEdges + singleNE;
end % for ii

%> 10, edge
[edge,~,indx2] = unique(totalEdges,'rows','legacy');

%> 18, interEdgeIndex, [NinterEdges x 1], the index of the interior edges
[~, ~, Es] = find( sparse(totalEdges(:,2), totalEdges(:,1), 1) );
interEdgeIndex = find(Es==2);

%%
% to get the elem2edge 
% and the "Matrix of edges" 
% and the diameter of each element, we choose the maximum distance between two vertices as the diameter.
NtotalEdges = 0; % here we need the NallEdges as the count index.
MofEdges = zeros(Nelems, max_singleNE); % to get the following edge2elem.
for ii = 1:Nelems
    singleElem = elem{ii,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges(nodes) on singleElem.
    
    %% get elem2edge
    elem2edge{ii} = indx2(NtotalEdges+1:NtotalEdges+singleNE)'; 
        %> fix the ii, then elem2edge{ii} is [1 x singleNE].
    NtotalEdges = NtotalEdges + singleNE;
    
    %% get the "matrix of edges"
    MofEdges(ii,1:singleNE) = elem2edge{ii}; % to get the following edge2elem.
    
end % for ii


%> 13, Nedges
Nedges = size(edge,1);


%% %> 15, edge2elem
edge2elem = zeros(Nedges, 4);
elem_nums = (1:Nelems)' * ones(1,max_singleNE); % [Nelems x max_singleNE], each column is 1:Nelems.
local_edge_nums = ones(Nelems, 1) * (1:max_singleNE); % [Nelems x max_singleNE], each row is 1:max_singleNE.
for ii = 1:Nedges
    elem_indx = elem_nums(MofEdges == ii);
    local_edge_indx = local_edge_nums(MofEdges == ii);
    if length(elem_indx)==2
        edge2elem(ii,1:2) = elem_indx;
        edge2elem(ii,3:4) = local_edge_indx;
    else
        edge2elem(ii,1) = elem_indx;
        edge2elem(ii,3) = local_edge_indx;
    end
end % for ii = 1: Nedges

g.node = node; %  [Nnodes x 2].
g.elem = elem; %  cell-type, elem{ii} is [1 x singleNE].
g.Nnodes = Nnodes; % 
g.Nelems = Nelems; % 
g.baryElem = baryElem; % 5, [Nelems x 2].
g.elem2edge = elem2edge; % elem2edge{ii} is [1 x singleNE].
g.edge = edge; % [Nedges x 2].
g.Nedges = Nedges; % 
g.edge2elem = edge2elem; %  [Nedges x 4].
g.interEdgeIndex = interEdgeIndex; %>  interEdgeIndex, [NinterEdges x 1], the index of the interior edges

end % function

function patchPlotMeshForMergeElems(node, elem)
%
%   YcZhang 24/9/2017
%   Last modified 24/9/2017
%
%

% Plot the mesh.
xx = node(:,1);
yy = node(:,2);

ElemNum = size(elem,1);

reset(cla), reset(clf), hold on

for CurrElem = 1 : ElemNum
    if iscell(elem)
        CurrNodes = elem{CurrElem};
    else
        CurrNodes = elem(CurrElem, :);
    end
    
    L = ones(1,length(CurrNodes));
    %patch(xx(CurrNodes), yy(CurrNodes), -0.005*L, 'w')
    patch(xx(CurrNodes), yy(CurrNodes), -0.005*L, [.95, .95, .95])
        %> also may use: patch(x, y, [r g b]) to control the color.
end

axis equal
axis([min(xx)-0.1 max(xx)+0.1 min(yy)-0.1 max(yy)+0.1])

end % function
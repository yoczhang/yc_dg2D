function g = polyMeshAuxStructure(node, elem)
% function g = polyMeshAuxStructure
%
%   [22/8/2017], add concavePointElem, see the details in the following %>21 .
%   [20/8/2017], add the centroid(xing xin) of each element.
%   [3/8/2017], add the diameter of each element. Here we choose the 
%   maximum distance between two vertices as the diameter.
%   [1/8/2017], add edge2elem.
%
%   In this function we get the auxiliary mesh structure based on the
%   polygon mesh structure Node and Elem.
%
%
%   We let Nnodes denote the number of nodes of the mesh Th,
%               Nedges denote the number of the EDGEs of Th,
%               Nelems denote the number of the elements of Th.
%
%
%   input:
%       node, [Nnodes x 2], the (x-coordinate, y-coordinate). or [2 x Nnodes], the (x-coordinate; y-coordinate). 
%       elem, cell-type data, [Nelems x 1], each row is the cell-type data, the points' order must be arranged in counter-clockwise.
%
%   output:
%       g.
%
%   Attention:
%       To save space all the data type in T is uint32. When use them as a input
%       of sparse(i,j,s,m,n), please change them into double type.
%
%
%	YcZhang 9/7/2017
%
%   Last modified 22/8/2017
%

% % clear
% % close all
% % [node,elem] = PolyMesher(@MbbDomain,30,6);

%%
% First we need to determine whether the input 'elem' is a cell-type.
% if 'elem' is not the cell-type, we need to tranfer it to cell-type.


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
%> 3, Nnodes
Nnodes = size(node, 1);

%> 4, Nelems
Nelems = size(elem, 1);

%> 5, to get barycenter of elem, baryElem
baryElem = zeros(Nelems, 2);

%> 6, elem2edge
elem2edge = cell(Nelems, 1);

%> 7, areaElem
areaElem = zeros(Nelems, 1);

%> 8, areaTri0Elem
areaTri0Elem = cell(Nelems, 1);

%>9, nuEdge0Elem
nuEdge0Elem = cell(Nelems, 1); % 9, nuEdge0Elem{ii} is [2 x singleNE], i.e., nuEdge0Elem{ii}(:,n), the norm-vector of n-th edge in ii-th elem.
 
%>21, concavePointElem
concavePointElem = zeros(Nelems,1); % In the program, we request that there at most has ONE concave point in each elem.
    %> 21, p=concavePointElem(n), stands for in the n-th elem, 
    %> if there has ONE concave point, then the concave point is the p-th point of the n-th elem,
    %> if there doesnot have concave point, then p is 0.

%%
% we first need to get the number of allEdges TO  allocate memory of
% the following varibale 'allEdges'.
max_singleNE = 0; % to get the following: 15, edge2elem.
NtotalEdges = 0;
for ii = 1:Nelems
    singleElem = elem{ii,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = length(singleElem); % get the number of edges on singleElem.
    NtotalEdges = NtotalEdges + singleNE;
    
    [nrow,ncol] = size(singleElem);
    if nrow > ncol
        elem{ii,1} = singleElem';
    end
    
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
    singleNE = length(singleElem); % get the number of edges(nodes) on singleElem.
    
    %% get the barycenter of elem
    baryElem(ii,1) = sum( node(singleElem,1) )/singleNE;
    baryElem(ii,2) = sum( node(singleElem,2) )/singleNE;
    
    %% 
    % get the ONE concave-point(凹点 ao dian)
    pointsElem = node(singleElem,:);
    p_forward = zeros(size(pointsElem));
    p_forward(1,:) = pointsElem(singleNE,:);
    p_forward(2:singleNE,:) = pointsElem(1:singleNE-1,:);
    p_backward = zeros(size(pointsElem));
    p_backward(1:singleNE-1,:) = pointsElem(2:singleNE,:);
    p_backward(singleNE,:) = pointsElem(1,:);
    vec1 = p_forward - pointsElem;
    vec2 = p_backward - pointsElem;
    
    % get the all edges of each elem
    singleElem2E = zeros(singleNE,2);
    for jj = 1:singleNE-1
        % get the all edges of each elem
        singleElem2E(jj,:) = singleElem([jj,jj+1]);
        
        % get the ONE concave-point
        r_temp = det([vec1(jj,:); vec2(jj,:)]);   
        if r_temp>0
            % request that the points' order must be arranged in counter-clockwise,
            % if in the clockwise, here should use r_temp<0.
            concavePointElem(ii) = jj; % concave-points (ao dian)         
        end
    end % for jj
    singleElem2E(singleNE,:) = singleElem([singleNE,1]);
    
    r_temp = det([vec1(singleNE,:); vec2(singleNE,:)]);
    if r_temp>0
        % request that the points' order must be arranged in counter-clockwise,
        % if in the clockwise, here should use r_temp<0.
        concavePointElem(ii) = singleNE; % concave-points (ao dian)         
    end
    
    
    %% get areaTri0Elem
    coordv1x = node(singleElem2E(:,1), 1);  coordv1y = node(singleElem2E(:,1), 2); % [singleNE x 1]
    coordv2x = node(singleElem2E(:,2), 1);  coordv2y = node(singleElem2E(:,2), 2); % [singleNE x 1]
        %> we let the two vertex of the 1-th and 2-th points of the small triangle,
        %> and let the barycenter of the polygon as the 3-th point.
% %     coordv3x = ones(singleNE,1)*baryElem(ii,1);  coordv3y = ones(singleNE,1)*baryElem(ii,2); 
% %         %> % [singleNE x 1], here we let the barycenter as the 3-th point of the triangle
% %     
% %     areaTri0Elem{ii} = 0.5*(coordv1x.*coordv2y + coordv2x.*coordv3y + coordv3x.*coordv1y ...
% %         - coordv1x.*coordv3y - coordv2x.*coordv1y - coordv3x.*coordv2y)'; % [1 x singleNE]
% %         %> If we know the coordiantes of vertex of triangle (v1x, v1y), (v2x, v2y), (v3x, v3y), then 
% %         %> the area of triangle = 0.5*(v1x*v2y + v2x*v3y + v3x*v1y - v1x*v3y - v2x*v1y - v3x*v2y).
% %     areaElem(ii) = sum(areaTri0Elem{ii});
    
    %% get nuEdge0Elem
    vectorEdge = [coordv2x-coordv1x, coordv2y-coordv1y]; % [singleNE x 2]
    nuEdge0Elem_temp = [vectorEdge(:,2), -vectorEdge(:,1)]; % [singleNE x 2]
    areaEdge0Elem = sqrt(vectorEdge(:,1).^2+vectorEdge(:,2).^2); % [singleNE x 1]
    arryAreaEdge0Elem = [areaEdge0Elem, areaEdge0Elem]; % [singleNE x 2]
    
    nuEdge0Elem{ii} = ( nuEdge0Elem_temp./arryAreaEdge0Elem )'; % [2 x singleNE]
    
    %% rearrange the every row of singleElem2E in order from smallest to largest. 
    singleElem2E = sort(singleElem2E,2);
    
    %% get the totalEdges
    totalEdges(NtotalEdges+1:NtotalEdges+singleNE, :) = singleElem2E;
    NtotalEdges = NtotalEdges + singleNE;
end % for ii

%> 10, edge
[edge,indx1,indx2] = unique(totalEdges,'rows','legacy');

%> 11, bdEdge
[Ei, Ej, Es] = find( sparse(totalEdges(:,2), totalEdges(:,1), 1) );
bdEdge = [Ej(Es==1), Ei(Es==1)];

%> 18, interEdgeIndex, [NinterEdges x 1], the index of the interior edges
interEdgeIndex = find(Es==2);

%> 19, bdEdgeIndex, [NinterEdges x 1], the index of the boundary edges
bdEdgeIndex = find(Es==1);

% % %> 9, bdEdge, because codegen does not support sparse function, so we chang
% % %> to the follwing the codes:
% % maxIndx = max(totalEdges(:,2));
% % matrix_temp = zeros(maxIndx,maxIndx);
% % loop_temp = size(totalEdges,1);
% % for ii = 1:loop_temp
% %     matrix_temp(totalEdges(ii,2), totalEdges(ii,1)) = ...
% %         matrix_temp(totalEdges(ii,2), totalEdges(ii,1)) + 1;
% % end % for ii
% % [Ei_1, Ej_1, Es_1] = find( matrix_temp );
% % bdEdge = [Ej_1(Es_1==1), Ei_1(Es_1==1)];


%%
% to get the elem2edge 
% and the "Matrix of edges" 
% and the diameter of each element, we choose the maximum distance between two vertices as the diameter.
NtotalEdges = 0; % here we need the NallEdges as the count index.
MofEdges = zeros(Nelems, max_singleNE); % to get the following edge2elem.
diameters = zeros(Nelems, 1); % 16, the diameters of each element.
centroidElem = zeros(Nelems,2); % 20, the centroid(xing xin) of each element.
hElem = zeros(Nelems, 1); 
    %> 22, the so-called 'h' of ecah element, in each element, h is defined by:
    %> the max distance between centroid(xing xin) and the element vertices.
    %> And hElem is used in the local bases, replacing the diameters used in the local bases.
for ii = 1:Nelems
    singleElem = elem{ii,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = length(singleElem); % get the number of edges(nodes) on singleElem.
    
    %% get elem2edge
    elem2edge{ii} = indx2(NtotalEdges+1:NtotalEdges+singleNE)'; 
        %> fix the ii, then elem2edge{ii} is [1 x singleNE].
    NtotalEdges = NtotalEdges + singleNE;
    
    %% get the "matrix of edges"
    MofEdges(ii,1:singleNE) = elem2edge{ii}; % to get the following edge2elem.
    
    %% 16, get the diamerters of each element.
    %% 20, get the centroid(xing xin) of each element.
    coordv = node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    baryElem_ii = baryElem(ii,:); % [1 x 2], the barycenter (x-coord, y-coord) of element.
    concavePointElem_ii = concavePointElem(ii); % [1 x 1], the concavePointElem of ii-the elem.
    coordTri0Elem = getcoordTri0Elem(singleNE, concavePointElem_ii, baryElem_ii, coordv);
    centroid_x = 0; centroid_y = 0;
    areaTri0Elem{ii} = zeros(1,size(coordTri0Elem,1)/3); % areaTri0Elem{ii} is [1 x size(coordTri0Elem,1)/3].
    for nt = 1:(size(coordTri0Elem,1)/3)
        % get the diameters
        vertices0singleElem = node(singleElem,:); % [singleNE x 2]
        firstVertices = ones(singleNE,1)*vertices0singleElem(nt,:);
        distances = (vertices0singleElem(:,1) - firstVertices(:,1)).^2 + (vertices0singleElem(:,2) - firstVertices(:,2)).^2;
        diameters(ii) = max([sqrt(distances);diameters(ii)]);
        
        % get the centroid
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt);
        centroid_x = centroid_x+phyGpoints(:,1)'*phyGweights;
        centroid_y = centroid_y+phyGpoints(:,2)'*phyGweights;
        
        % get areaTri0Elem
        areaTri0Elem{ii}(1,nt) = sum(phyGweights);
    end % for nt
    areaElem(ii) = sum(areaTri0Elem{ii});
    centroidElem(ii,1) = centroid_x/areaElem(ii);
    centroidElem(ii,2) = centroid_y/areaElem(ii);
    
    hElem(ii) = sqrt( max( sum((ones(singleNE,1)*centroidElem(ii,:)-coordv).^2,2) ) );
    
end % for ii


%% 
%> 12, areaEdge
areaEdge = sqrt( (node(edge(:,1),1) - node(edge(:,2),1)).^2 ...
    + (node(edge(:,1),2) - node(edge(:,2),2)).^2 ); % [Nedges x 1]

%> 13, Nedges
Nedges = size(edge,1);

%> 14, baryEdge
baryEdge = zeros(Nedges, 2);
baryEdge(:,1) = ( node(edge(:,1),1) + node(edge(:,2),1) )/2;
baryEdge(:,2) = ( node(edge(:,1),2) + node(edge(:,2),2) )/2;

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
%> edge2elem, [Nedges x 4]. 
%> in the k-th row, if there has 0, then stands for BoundaryEdge.
%> such as the k-th row: 
%> [n1, n2]=edge2elem(k,1:2) stands for the n1-th elem and n2-th elem share the k-th edge.
%> local_e1=edge2elem(k,3) stands for the local edge index in n1-th elem of the k-th edge .
%> also, local_e2=edge2elem(k,4) stands for the local edge index in n2-th elem of the k-th edge .

%% %> 17, the map, from ref [0,1] to phy edges
mapRefE2PhyE = @(s, P1, P2) s.*P1+(1-s).*P2;
    %> input: s, the ref Gauss points on [0,1].
    %>          P1, P2, [Npoints x 1], the two x(or y)-coordiants of Physical edge.
    %> output: [Npoints x 1], the x(or y) Gauss points on Physical e.



%% Method I: Construct the structure-type mesh g.
% % g = struct( 'node', node, ... % [Nnodes x 2].
% %     'elem', elem, ... % cell-type, elem{ii} is [1 x singleNE].
% %     'Nnodes', Nnodes, ... % 1
% %     'Nelems', Nelems, ...
% %     'baryElem', baryElem, ... % 3, [Nelems x 1].
% %     'elem2edge', elem2edge, ... % 4, elem2edge{ii} is [1 x singleNE].
% %     'areaElem', areaElem, ...
% %     'areaTri0Elem', areaTri0Elem, ... % 6, areaTri0Elem{ii} is [1 x singleNE].
% %     'nuEdge0Elem', nuEdge0Elem, ... % 7, nuEdge0Elem{ii} is [2 x singleNE].
% %     'edge', edge, ... % 8, [Nedges x 2]
% %     'bdEdge', bdEdge, ... % 9
% %     'areaEdge', areaEdge, ... % 10
% %     'Nedges', Nedges, ... % 11
% %     'baryEdge', baryEdge ... % 12, [Nedges x 2]
% %     );

g.node = node; % 1, [Nnodes x 2].
g.elem = elem; % 2, cell-type, elem{ii} is [1 x singleNE].
g.Nnodes = Nnodes; % 3
g.Nelems = Nelems; % 4
g.baryElem = baryElem; % 5, [Nelems x 2].
g.elem2edge = elem2edge; % 6, elem2edge{ii} is [1 x singleNE].
g.areaElem = areaElem; % 7 
g.areaTri0Elem = areaTri0Elem; % 8, areaTri0Elem{ii} is [1 x singleNE].
g.nuEdge0Elem = nuEdge0Elem; % 9, nuEdge0Elem{ii} is [2 x singleNE], i.e., nuEdge0Elem{ii}(:,n), the norm-vector of n-th edge in ii-th elem.
g.edge = edge; % 10, [Nedges x 2].
g.bdEdge = bdEdge; % 11, [NbdEdge x 2], the two-points-global-index of the edge.
g.areaEdge = areaEdge; % 12
g.Nedges = Nedges; % 13
g.baryEdge = baryEdge; % 14, [Nedges x 2].
g.edge2elem = edge2elem; % 15, [Nedges x 4].
    %> edge2elem, [Nedges x 4]. 
    %> in the k-th row, if there has 0, then stands for BoundaryEdge.
    %> such as the k-th row: 
    %> [n1, n2]=edge2elem(k,1:2) stands for the n1-th elem and n2-th elem share the k-th edge.
    %> local_e1=edge2elem(k,3) stands for the local edge index in n1-th elem of the k-th edge .
    %> also, local_e2=edge2elem(k,4) stands for the local edge index in n2-th elem of the k-th edge .
g.diameters = diameters; % 16, [Nelems x 1], the diameter of each element.
    %> we choose the maximum distance between two vertices as the diameter.
g.mapRefE2PhyE = mapRefE2PhyE; % 17, the map, from ref [0,1] to phy edges.
    %> input: s, the ref Gauss points on [0,1].
    %>          P1, P2, [Npoints x 1], the two x(or y)-coordiants of Physical edge.
    %> output: [Npoints x 1], the x(or y) Gauss points on Physical e.
g.interEdgeIndex = interEdgeIndex; %> 18, interEdgeIndex, [NinterEdges x 1], the index of the interior edges
g.bdEdgeIndex = bdEdgeIndex; %> 19, bdEdgeIndex, [NinterEdges x 1], the index of the boundary edges
g.centroidElem = centroidElem; %>20, centroidElem, [Nelems x 2], the centroid(xing xin) of each elem.
g.concavePointElem = concavePointElem; 
    %>21, p=concavePointElem(n), stands for in the n-th elem, 
    %> if there has ONE concave point, then the concave point is the p-th point of the n-th elem,
    %> if there doesnot have concave point, then p is 0.
g.hElem = hElem;


%---------------------------


% % %% Method II: Construct the cell-type mesh g.
% % %> For example, if we want to get the var node, node = g{g{1}.node}.
% % g = cell(15,1);
% % varName = struct('node', 2, ...
% %     'elem', 3, ...
% %     'Nnodes', 4, ...
% %     'Nelems', 5, ...
% %     'baryElem', 6, ...
% %     'elem2edge', 7, ...
% %     'areaElem', 8, ...
% %     'areaTri0Elem', 9, ...
% %     'nuEdge0Elem', 10, ...
% %     'edge', 11, ...
% %     'bdEdge', 12, ...
% %     'areaEdge', 13, ...
% %     'Nedges', 14, ...
% %     'baryEdge', 15 ...
% %     );
% % 
% % g{1} = varName;
% % g{2} = node; % [Nnodes x 2].
% % g{3} = elem; % cell-type, elem{ii} is [1 x singleNE].
% % g{4} = Nnodes; % 1
% % g{5} = Nelems; % 2
% % g{6} = baryElem; % 3, [Nelems x 1].
% % g{7} = elem2edge; % 4, elem2edge{ii} is [1 x singleNE].
% % g{8} = areaElem;
% % g{9} = areaTri0Elem; % 6, areaTri0Elem{ii} is [1 x singleNE].
% % g{10} = nuEdge0Elem; % 7, nuEdge0Elem{ii} is [2 x singleNE].
% % g{11} = edge; % 8, [Nedges x 2].
% % g{12} = bdEdge; % 9
% % g{13} = areaEdge; % 10
% % g{14} = Nedges; % 11
% % g{15} = baryEdge; % 12, [Nedges x 2].


end % function


%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function coordTri0Elem = getcoordTri0Elem(singleNE, beginP_n, baryElem, coordv)
%
%
%   input:
%       beginP_n, the n-th point of singleElem, and from the n-th point to
%       construct the little triangles.
%

m = @(x) mod(x,singleNE)+(x==singleNE)*singleNE;
n = beginP_n;

if n == 0
    if singleNE == 3
        coordTri0Elem = coordv; % [3 x 2]
    elseif singleNE == 4
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 5
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(1,:)];
    elseif singleNE == 6
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 7
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 8
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(8,:);
            baryElem; coordv(8,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 9
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(8,:);
            baryElem; coordv(8,:); coordv(9,:);
            baryElem; coordv(9,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 10
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(8,:);
            baryElem; coordv(8,:); coordv(9,:);
            baryElem; coordv(9,:); coordv(10,:);
            baryElem; coordv(10,:); coordv(1,:)]; % [3*singleNE x 2]
    end % if singleNE == 3
    
else
    if singleNE == 3
        coordTri0Elem = coordv; % [3 x 2]
    elseif singleNE == 4
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:)];
    elseif singleNE == 5
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:)];
    elseif singleNE == 6
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:)];
    elseif singleNE == 7
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:)];  
    elseif singleNE == 8
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:); ...
            coordv(m(n),:); coordv(m(n+6),:); coordv(m(n+7),:)];  
    elseif singleNE == 9
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:); ...
            coordv(m(n),:); coordv(m(n+6),:); coordv(m(n+7),:); ...
            coordv(m(n),:); coordv(m(n+7),:); coordv(m(n+8),:)];
    elseif singleNE == 10
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:); ...
            coordv(m(n),:); coordv(m(n+6),:); coordv(m(n+7),:); ...
            coordv(m(n),:); coordv(m(n+7),:); coordv(m(n+8),:); ...
            coordv(m(n),:); coordv(m(n+8),:); coordv(m(n+9),:)];
    end % if singleNE == 3
    
end % if n==0

end % function getcoordTri0Elem
%%<<-- End sub function 1 ---------------------------------------------------------------


%%>> -- Begin sub function 2 -------------------------------------------------------------
function [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt)
%
%   output:
%       phyGpoints, [Npoints x 2]
%       phyGweights, [Npoints x 1]
%
%	YcZhang 20/8/2017
%   Last modified 20/8/2017
%
x1=coordTri_nt(1,1);
y1=coordTri_nt(1,2);
x2=coordTri_nt(2,1);
y2=coordTri_nt(2,2);
x3=coordTri_nt(3,1);
y3=coordTri_nt(3,2);
JacobiTri=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));

refX = [0.063089014491502, 0.873821971016996, 0.063089014491502, ...
    0.249286745170910, 0.501426509658179, 0.249286745170910, ...
    0.310352451033785, 0.053145049844816, 0.636502499121399, ...
    0.053145049844816, 0.636502499121399, 0.310352451033785];
refY = [0.063089014491502, 0.063089014491502, 0.873821971016996, ...
    0.249286745170910, 0.249286745170910, 0.501426509658179, ...
    0.053145049844816, 0.310352451033785, 0.053145049844816, ...
    0.636502499121399, 0.310352451033785, 0.636502499121399];
refW  = [0.025422453185103, 0.025422453185103, 0.025422453185103, ...
    0.058393137863189, 0.058393137863189, 0.058393137863189, ...
    0.041425537809187, 0.041425537809187, 0.041425537809187, ...
    0.041425537809187, 0.041425537809187, 0.041425537809187];

phyGweights = JacobiTri * refW';
phyGpoints(:,1)=x1+(x2-x1)*refX' + (x3-x1)*refY';
phyGpoints(:,2)=y1+(y2-y1)*refX' + (y3-y1)*refY';
end % function getGaussLocalTri

%%>> -- Begin sub function 2 -------------------------------------------------------------

%%<<-- End sub function 2 ---------------------------------------------------------------

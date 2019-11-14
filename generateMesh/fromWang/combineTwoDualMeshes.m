function [node, elem] = combineTwoDualMeshes(refine, elementType)

%% Two meshes
node1 = [0,1; 1,1; 1,2; 0,2]; elem1 = [2,3,1; 4,1,3];
for i = 1:refine, [node1, elem1] = uniformrefine(node1, elem1); end
[node1,elem1] = dualmesh(node1,elem1);

node2 = [0,0; 1,0; 1,1; 0,1]; elem2 = [2,3,1; 4,1,3];
for i = 1:refine, [node2, elem2] = uniformrefine(node2, elem2); end
[node2,elem2] = dualmesh(node2,elem2);

%% Find the intersect points
eps = 1e-10;
[Inte1, unused] = find(abs(node1(:,2) - 1)<=eps);
[Inte2, unused] = find(abs(node2(:,2) - 1)<=eps);

[unused, sortInte1] = sort(node1(Inte1,1)); Inte1 = Inte1(sortInte1);
[unused, sortInte2] = sort(node2(Inte2,1)); Inte2 = Inte2(sortInte2);


%% Modify by combining the overlapping points
NV1 = size(node1,1); NT1 = size(elem1,1); NV2 = size(node2,1); NT2 = size(elem2,1);
NV = NV1 + NV2 - length(Inte1); 

otherNode2 = setdiff((1:NV2)', Inte2);

node = zeros(NV,2); node(1:NV1,:) = node1; node(NV1+1:end,:) = node2(otherNode2,:);
elem = cell(NT1+NT2,1); elem(1:NT1) = elem1; elem(NT1+1:end) = elem2;

for n = 1:length(elem2)
    verts = elem2{n};
    for i = 1:length(verts)
        j = verts(i);
        if ismember(j,Inte2)
            [j2, unused] = find(Inte2==j);
            elem{n+NT1}(i) = Inte1(j2);
        else
            [j2, unused] = find(otherNode2==j);
            new = j2 + NV1;
            node(new,:) = node2(j,:);
            elem{n+NT1}(i) = new;
        end
    end
end

%% Plot
if strcmp(elementType, 'distortedPolygon')
    x = node(:,1); y = node(:,2);
    node(:,1) = x + 1/10 * sin(2*pi*x).*sin(2*pi*y);
    node(:,2) = y + 1/10 * sin(2*pi*x).*sin(2*pi*y);
end
% PolyMshr_PlotMsh(node, elem, length(elem));

    function PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load)
        axis equal; axis off; hold on;
        % Element = Element(1:NElem)';                 %Only plot the first block
        MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
        PadWNaN = @(E) [E' NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
        ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
        ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
    
        patch('Faces',ElemMat,'Vertices',Node,'FaceColor','w'); pause(1e-6)
        if exist('Supp','var')&&~isempty(Supp) %Plot Supp BC if specified
           plot(Node(Supp(:,1),1),Node(Supp(:,1),2),'b>','MarkerSize',8);
        end
        if exist('Load','var')&&~isempty(Load) %Plot Load BC if specified
           plot(Node(Load(:,1),1),Node(Load(:,1),2),'m^','MarkerSize',8);
        end
    end

end
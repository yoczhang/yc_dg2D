function [dNode,dElem] = dualmesh(node,elem)
% the input mesh should be convex triangle mesh or convex quad mesh
[NT,NV4e] = size(elem);
N = size(node,1); dNT = N; dN = NT;

dNode = zeros(NT,2);
for nd = 1:NV4e
    dNode = dNode + node(elem(:,nd),:);
end
dNode = dNode/NV4e; % the center of each element
[bdNode,bdEdge,isBdNode] = findboundary(elem);
NbEd = size(bdEdge,1);
dN = dN + NbEd;
bdEdge = myunique(bdEdge);
dNode(NT+1:dN,:) = (node(bdEdge(:,1),:) + node(bdEdge(:,2),:))/2;


E2n = sparse(repmat(1:NT,1,NV4e), elem, 1, NT, N);
valence = accumarray(elem(:),ones(NV4e*NT,1),[N 1]);
tNode = 1:N; tbdEdge = 1:NbEd;
intrNode = tNode(~isBdNode);
dElem = cell(dNT,1); 
for i = 1:numel(bdNode)
    cardi = valence(bdNode(i));
    nodeStar = find(E2n(:,bdNode(i)));
    dN1D = tbdEdge((bdEdge(:,1)==bdNode(i))|(bdEdge(:,2)==bdNode(i)));
    dNode(dN+1,:) = node(bdNode(i),:);
    dN = dN + 1;
    dN4E = [nodeStar;dN1D'+NT;dN];
    order = ccworder(dNode(dN4E,:));
    dElem{bdNode(i)} = dN4E(order);
end

for i = 1:numel(intrNode)
    dN4E = find(E2n(:,intrNode(i)));
    order = ccworder(dNode(dN4E,:));
    dN4E = dN4E(order);
    dElem{intrNode(i)} = dN4E;
end

% x = dNode(:,1); y = dNode(:,2);
% dNode(:,1) = x + 1/10 * sin(2*pi*x).*sin(2*pi*y);
% dNode(:,2) = y + 1/10 * sin(2*pi*x).*sin(2*pi*y);
% PolyMshr_PlotMsh(dNode,dElem,dNT);
end

function order = ccworder(node)
    center = mean(node,1);
    a = atan2(node(:,2) - center(2), node(:,1) - center(1));
    [~, order] = sort(a);
end
function PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load)
clf; axis equal; axis off; hold on;
%Element = Element(1:NElem)';                 %Only plot the first block
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
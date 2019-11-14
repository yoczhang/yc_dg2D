function [node,elem] = non_convex_octagona_mesh(node,elem)
    % 5*5 quad mesh
    
    T = auxstructurequad(elem);
    N = size(node,1);
    edge = T.edge;
    NEd = size(edge,1); NT = size(elem,1);
    m4Ed = (node(edge(:,1),:) + node(edge(:,2),:))/2;
    elem2edge = T.elem2edge;
    elem = [elem elem2edge+N];
   
    %tmp = (1:NEd)';
    %IEd = tmp(edge2elem(:,1)~=edge2elem(:,2));
    nx = ceil(sqrt(NT)); ny = nx;
  %  idxmat = reshape(uint32(1:nx*ny), nx, ny);
 %  idxmat = flipud(transpose(reshape(uint32(1:nx*ny), nx, ny)));
    idxmat = transpose(reshape(uint32(1:nx*ny), nx, ny));
    hx = 1/nx/2; hy = 1/ny/2;
      % modified all the iterior midpoints
    ide1 = idxmat(1:nx-1,1:ny-1);
    ide1 = ide1(:); n1 = length(ide1);
    ed1 = elem2edge(ide1,3); ed2 = elem2edge(ide1,2);
    m4Ed(ed1,:) = node(elem(ide1,3),:) + 0.5*hx*[ones(n1,1),-ones(n1,1)];
    m4Ed(ed2,:) = node(elem(ide1,3),:) - 0.5*hx*[ones(n1,1),-ones(n1,1)];
    
    ide2 = idxmat(nx,1:ny-1);
    n2 = length(ide2);
    ed2 = elem2edge(ide2,3);
    m4Ed(ed2,:) = node(elem(ide2,3),:) + 0.5*hx*[ones(n2,1),-ones(n2,1)];
   
    ide3 = idxmat(1:nx-1,ny);
    n3 =  length(ide3);
    ed3 = elem2edge(ide3,2);
    m4Ed(ed3,:) = node(elem(ide3,3),:) - 0.5*hx*[ones(n3,1),-ones(n3,1)];
    
    node = [node;m4Ed];
    elem = elem(:,[3 6 2 5 1 8 4 7]);
    
    elem = reshape(elem',NT*8,1);idx = 8*ones(NT,1);
    elem = mat2cell(elem,idx);

   figure(2)
   PolyMshr_PlotMsh(node,elem,NT);
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
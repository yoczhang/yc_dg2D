function [nodeNew, elemNew] = cutTriangularMesh(refine)


%% node, elem
domain = [0 0; 1 0; 1 2; 0 2]; eps = 1e-10; Y = 1; % the interface
gd = [2; size(domain,1); domain(:,1); domain(:,2)]; gl = decsg(gd);
[p, e, t] = initmesh(gl); node = p'; elem = t(1:3,:)'; 
for i = 1:refine
    [node, elem] = uniformrefine(node, elem);
end
clear domain gd gl p e t

% [node, elem] = squaremesh([0,1,0,2], 1/9); eps = 1e-10; Y = 1; % the interface

%% edge, left, righ, righPlusLeft, righMinuLeft, t
[~, edge] = dofedge(elem);
NE = size(edge,1); N = size(node,1); NT = size(elem,1);

left = node(edge(:,1),:); righ = node(edge(:,2),:); 
righPlusLeft = (righ + left)/2; righMinuLeft = (righ - left)/2;


%% inte, objecCoor
allEdge = false(NE,1); allEdge((left(:,2)==righ(:,2))) = true;
horiEdge = find(allEdge); otherEdge = find(~allEdge); % find the non-vertical edges
t = (Y - righPlusLeft(otherEdge,2)) ./ (righMinuLeft(otherEdge,2)); % for otheredge
clear allEdge

[used, unused] = find(abs(t)<1-eps); objecEdge = otherEdge(used); 
% objec is the edge cut by the interface y = 0
clear unused

objecEdgeCoor = [righMinuLeft(objecEdge,1) .* t(used) + righPlusLeft(objecEdge,1) Y*ones(length(objecEdge),1)];

[objecNode, unused] = find(abs(node(:,2) - Y) <= eps); 

objecCoorOld = [node(objecNode,:); objecEdgeCoor];
[i, j] = sort(objecCoorOld(:,1)); objecCoor = objecCoorOld(j,:); 
objecN = size(objecCoor,1); clear i j

inte = zeros(objecN,3); addN = N;
for i = 1:objecN
    coor = objecCoor(i,:);
    disCor = sqrt(sum((objecCoorOld - repmat(coor,size(objecCoorOld,1),1)).^2, 2));
    [used, unused] = find(disCor <= eps);
    if used<=length(objecNode)
        inte(i,1) = 1; inte(i,2) = objecNode(used);
    else
        addN = addN + 1; 
        inte(i,1) = 2; inte(i,2) = objecEdge(used-length(objecNode)); inte(i,3) = addN;
    end
    clear used unused disCor coor
end
clear objecCoorOld objecEdgeCoor objecNode objecEdge

%% Modify
objecCenter = [(objecCoor(1:end-1,1)+objecCoor(2:end,1))/2 (objecCoor(1:end-1,2)+objecCoor(2:end,2))/2];
elemOld = elem; addNT = size(objecCenter,1); addElem = cell(addNT,1);
for i = 1:addNT
    inin = zeros(NT,1);
    for n = 1:NT
        inin(n) = inpolygon(objecCenter(i,1), objecCenter(i,2), node(elemOld(n,:),1)', node(elemOld(n,:),2)');
%         if in==1
%             break;
%         end
    end
    [n, unused] = find(inin);
    % the current point is in the n-th element
    
    j1 = inte(i,1); j2 = inte(i+1,1);
    if j1==1&&j2==2
        n1 = inte(i,2);
        [unused, l1] = find(elemOld(n,:)==n1);
        if l1==1
            elem(n,:) = [n1, inte(i+1,3), elemOld(n,3)]; addElem{i} = [n1, elemOld(n,2), inte(i+1,3)];
        elseif l1==2
            elem(n,:) = [n1, inte(i+1,3), elemOld(n,1)]; addElem{i} = [n1, elemOld(n,3), inte(i+1,3)];
        elseif l1==3
            elem(n,:) = [n1, inte(i+1,3), elemOld(n,2)]; addElem{i} = [n1, elemOld(n,1), inte(i+1,3)];
        end
    elseif j1==2&&j2==1
        n1 = inte(i+1,2);
        [unused, l1] = find(elemOld(n,:)==n1); 
        if l1==1
            elem(n,:) = [n1, elemOld(n,2), inte(i,3)]; addElem{i} = [n1, inte(i,3), elemOld(n,3)];
        elseif l1==2
            elem(n,:) = [n1, elemOld(n,3), inte(i,3)]; addElem{i} = [n1, inte(i,3), elemOld(n,1)];
        elseif l1==3
            elem(n,:) = [n1, elemOld(n,1), inte(i,3)]; addElem{i} = [n1, inte(i,3), elemOld(n,2)];
        end
    elseif j1==2&&j2==2
        n1 = inte(i,3); n2 = inte(i+1,3);
        flag = sign(node(elemOld(n,:),2) - Y) + 1;
        [flagi, unused] = find(abs(flag-2) <= eps); 
        if length(flagi)==1
            if flagi==1
                elem(n,:) = [elemOld(n,flagi), n1, n2]; addElem{i} = [n1, elemOld(n,2), elemOld(n,3), n2];
            elseif flagi==2
                elem(n,:) = [elemOld(n,flagi), n1, n2]; addElem{i} = [n1, elemOld(n,3), elemOld(n,1), n2];
            elseif flagi==3
                elem(n,:) = [elemOld(n,flagi), n1, n2]; addElem{i} = [n1, elemOld(n,1), elemOld(n,2), n2];
            end
        elseif length(flagi)==2
            flagj = setdiff(1:3, flagi);
            if flagj==1
                elem(n,:) = [elemOld(n,flagj), n2, n1]; addElem{i} = [n2, elemOld(n,2), elemOld(n,3), n1];
            elseif flagj==2
                elem(n,:) = [elemOld(n,flagj), n2, n1]; addElem{i} = [n2, elemOld(n,3), elemOld(n,1), n1];
            elseif flagj==3
                elem(n,:) = [elemOld(n,flagj), n2, n1]; addElem{i} = [n2, elemOld(n,1), elemOld(n,2), n1];
            end
        end
    end
end

elemNew = cell(NT+addNT,1);
for i = 1:NT
    elemNew{i} = elem(i,:)';
end
for i = 1:addNT
    elemNew{i+NT} = addElem{i}';
end

[i, j] = find(inte(:,1)==2);
nodeNew = [node; objecCoor(i, :)]; 

% figure; 
% subplot(1,3,1); PolyMshr_PlotMsh(nodeNew, elemNew, NT+addNT, [], [], 1, NT+addNT);
% subplot(1,3,2); PolyMshr_PlotMsh(nodeNew, elemNew, NT+addNT, [], [], 1, NT);
% subplot(1,3,3); PolyMshr_PlotMsh(nodeNew, elemNew, NT+addNT, [], [], NT+1, NT+addNT);

    function PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load,a,b)
        axis equal; axis off; hold on;
        Element = Element(a:b)';                 %Only plot the first block
        MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
        PadWNaN = @(E) [E' NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
        ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
        ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
        
        patch('Faces',ElemMat,'Vertices',Node,'FaceColor','w'); pause(1e-6)
%         if exist('Supp','var')&&~isempty(Supp) %Plot Supp BC if specified
%            plot(Node(Supp(:,1),1),Node(Supp(:,1),2),'b>','MarkerSize',8);
%         end
%         if exist('Load','var')&&~isempty(Load) %Plot Load BC if specified
%             plot(Node(Load(:,1),1),Node(Load(:,1),2),'m^','MarkerSize',8);
%         end
    end
end
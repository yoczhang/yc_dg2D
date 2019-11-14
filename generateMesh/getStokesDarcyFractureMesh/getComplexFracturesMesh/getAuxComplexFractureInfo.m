function meshInfo = getAuxComplexFractureInfo(meshInfo)
%
%
%	YcZhang 19/1/2019
%

    %-
    %- get the unit-normal-vector of fracture edges, other edges set 0;
    %- and we get the rearranged fractureFace.
    fractureFace = [];
    nuFracture = [];
    fractureBeginEndPoint = [];
    fractureFaceNode = [];
    splitedFractureFace = meshInfo.splitedFractureFace;
    NsplitedF = length(splitedFractureFace);
    for ii = 1:NsplitedF
        currFrac = splitedFractureFace{ii};
        [R,C] = size(currFrac);
        if R<C
            currFrac = currFrac'; % let the currFrac is the column vector
        end
        [currNuFra,currFraFace,currFraBeginEndPoint,currFraFaceNode] = ...
            getComplexFactureNormalVector(currFrac,meshInfo.edge,meshInfo.node);
        
        temp = [nuFracture; currNuFra];
        nuFracture = temp;
        
        temp = [fractureFace; currFraFace];
        fractureFace = temp;
        
        temp = [fractureBeginEndPoint; currFraBeginEndPoint];
        fractureBeginEndPoint = temp;
        
        temp = [fractureFaceNode; currFraFaceNode];
        fractureFaceNode = temp;
    end % for ii
    
    meshInfo.nuFracture = nuFracture;
    meshInfo.fractureFace = fractureFace;
    meshInfo.fractureBeginEndPoint = fractureBeginEndPoint;
    meshInfo.fractureFaceNode = fractureFaceNode;
    
    
    %-------------------------------
    fractureFace = meshInfo.fractureFace;
    
    %% following is the auxiliary fracture informations
    FractureElem = meshInfo.edge2elem(fractureFace,1:2);
    isFractureElem = false(meshInfo.Nelems,1);
    isFractureElem(FractureElem(:)) = true;
    meshInfo.isFractureElem = isFractureElem;
    
    %- whether add edge in element
    Nedges = meshInfo.Nedges;
    map_edge_2_addedge = zeros(Nedges,2);
    map_edge_2_addedge(:,1) = 1:Nedges;
    map_edge_2_addedge(:,2) = 1:Nedges;
    map_edge_2_addedge(fractureFace,2) = Nedges + (1:length(fractureFace));
    
    edge2elem = meshInfo.edge2elem;
    addedge2elem = edge2elem(fractureFace,1:2);
    addedge_in_elem = max(addedge2elem,[],2); % we let the bigger elem to add the edge.
    isAddedge_in_elem = false(meshInfo.Nelems,1);
    isAddedge_in_elem(addedge_in_elem) = true;
    
    meshInfo.map_edge_2_addedge = map_edge_2_addedge;
    meshInfo.isAddedge_in_elem = isAddedge_in_elem;
    
% %     %- get the unit-normal-vector of fracture edges, other edges set 0;
% %     %- and we get the rearranged fractureFace.
% %     [nuFracture,fractureFace,fractureBeginEndPoint,fractureFaceNode] = ...
% %         getComplexFactureNormalVector(fractureFace,meshInfo.edge,meshInfo.node);
% %     meshInfo.nuFracture = nuFracture;
% %     meshInfo.fractureFace = fractureFace;
% %     meshInfo.fractureBeginEndPoint = fractureBeginEndPoint;
% %     meshInfo.fractureFaceNode = fractureFaceNode;

end % function



%% sub functions
function [nuFracture, rearrange_fractureFace,fractureBeginEndPoint, fractureFaceNode] = ...
    getComplexFactureNormalVector(fractureFace, edge, node)
%
%   Why we need to re-arrange the fractureFace?
%   We consider the following case:
%   The actual fracture faces are in order: 3-2-5-9-8, 
%   but, we store in the array structure fractureFace is: 2-3-5-8-9, (from small to big values).
%   So, we need to fix the order in fractureFace accord with the actual case.
%

fractureTotal = edge(fractureFace,:);
% oneVal = unique(fractureTotal(:),'stable');
oneVal = unique(fractureTotal(:));
counts_temp = histc(fractureTotal(:),oneVal);

bdNode = oneVal((counts_temp==1));
fractureBeginEndPoint = bdNode; % the bdNode(1) is begin point, bdNode(2) is end point. 

fractureFace_temp = [fractureFace, fractureFace];

currentNode = bdNode(1);
currentEdge = fractureFace_temp(fractureTotal == currentNode);

rearrange_edge = zeros(size(edge));
nuFracture = zeros(size(edge));
rearrange_fractureFace = zeros(size(fractureFace));
fractureFaceNode = zeros(length(fractureFace),2);

for CurrEdge = 1:length(fractureFace)
    if isempty(currentEdge) || isempty(currentNode)
        error('MyDebug: "currentEdge" or "currentNode" is empty')
    end
    
    rearrange_edge(currentEdge,1) = currentNode;
    temp_edge = edge(currentEdge,:);
    nextNode = setdiff(temp_edge,currentNode);
    rearrange_edge(currentEdge,2) = nextNode;
    
    fractureFaceNode(CurrEdge,1) = currentNode;
    fractureFaceNode(CurrEdge,2) = nextNode;
    
    rearrange_fractureFace(CurrEdge,1) = currentEdge;
    
    firstnodeXY = node(currentNode,:);
    secondnodeXY = node(nextNode,:);
    areaedge = sqrt( sum( (firstnodeXY-secondnodeXY).^2 ) );
    
    nuFracture(currentEdge,1) = -(secondnodeXY(2)-firstnodeXY(2))/areaedge;
    nuFracture(currentEdge,2) = (secondnodeXY(1)-firstnodeXY(1))/areaedge;
    
    
    nextEdge = fractureFace_temp(fractureTotal == nextNode);
    nextEdge = setdiff(nextEdge,currentEdge);
    
    currentEdge = nextEdge;
    currentNode = nextNode;
end

nuFracture = sparse(nuFracture);

end % function
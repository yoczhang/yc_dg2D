function [node,elem] = mrstG_2_myMeshInfo(mestG)
%   
%
%
%

mestG = sortEdges(mestG);
node = mestG.nodes.coords;
elem = cell(mestG.cells.num,1);
for nElem = 1:mestG.cells.num
    faceIndex = mestG.cells.faces(mestG.cells.facePos(nElem) : mestG.cells.facePos(nElem+1)-1, 1);
    
    nodeIndex = [];
    upperNodeIndex = [];
    for ii = 1:length(faceIndex)
        currNodeIndex = ...
            mestG.faces.nodes(mestG.faces.nodePos(faceIndex(ii)):mestG.faces.nodePos(faceIndex(ii)+1)-1);
        if ii == 1
            upperNodeIndex = currNodeIndex;
            continue;
        end
        intersectNodeIndex = intersect(upperNodeIndex,currNodeIndex);
        temp_1 = 1:2;
        
        if ii == 2
            if temp_1(upperNodeIndex==intersectNodeIndex) == 1
                aa = upperNodeIndex(1);
                upperNodeIndex(1) = upperNodeIndex(2);
                upperNodeIndex(2) = aa;
            end
            nodeIndex = upperNodeIndex;
        end
        
        if temp_1(currNodeIndex==intersectNodeIndex) == 2
            aa = currNodeIndex(1);
            currNodeIndex(1) = currNodeIndex(2);
            currNodeIndex(2) = aa;
        end
        
        upperNodeIndex = currNodeIndex;
        nodeIndex = [nodeIndex;
            currNodeIndex];
    end % for ii
    currElem = unique(nodeIndex','stable');
    elem{nElem} = currElem;
end % nElemD

end % function
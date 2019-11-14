% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DCrackTipElem(xxElem, ...
    yyElem, ffElem, xxTip, yyTip, nQ1D)

% Divide the element containing the crack tip into six integration triangles.
% Each triangle contains the crack tip as a node.
% 
%  4 +------------------+ 3        4 +------------------+ 3
%    |                  |            |--_     [5]    _- |
%    |                  |   ###      |   --_       _-   |
%    |          7       |###         |      --_ 7_- [4] |
%    |        X       ### 6          | [6]   _-X-_------+ 6
%    |             ###  |            |     _- /   -_[3] |
%    |          ###     |            |   _-  /      -_  |
%    |       ###        |            | _-[1]/   [2]   -_|
%  1 +----###-----------+ 2        1 +-----+------------+ 2
%      ### 5                               5
% 
% Input data are the real element coordinates, the level set function and the 
% crack tip position. Note, that the level set is assumed a straight line in 
% the reference element. But the crack tip position is defined in the real 
% domain (not in the reference domain). 

% Example call:
% % [x, y, w] = IntPoints2DCrackTipElem([-2 2 2 -2], [-2 -2 2 2], [-1 -1 1 1]', 0, 0, 5)
% xxElem = [5 -3 -3  5];
% yyElem = [2  2 -2 -2];
% ffElem = [1/7 -1 3 5];
% xxTip = 4; yyTip = 1; nQ1D = 3;
% [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DCrackTipElem(xxElem, ...
%     yyElem, ffElem, xxTip, yyTip, nQ1D)

% Extract the cut element segments of this element.
SignVect = sign(ffElem);

Count = 0;
if SignVect(1) ~= SignVect(2)
    Count = Count + 1; 
    CutSegm(Count) = 1;
    xxS(Count) = Interpolate(-1,  1, ffElem(1), ffElem(2));
    yyS(Count) = Interpolate(-1, -1, ffElem(1), ffElem(2));
end 
if SignVect(2) ~= SignVect(3)
    Count = Count + 1; 
    CutSegm(Count) = 2;
    xxS(Count) = Interpolate( 1, 1, ffElem(2), ffElem(3));
    yyS(Count) = Interpolate(-1, 1, ffElem(2), ffElem(3));
end 
if SignVect(3) ~= SignVect(4)
    Count = Count + 1; 
    CutSegm(Count) = 3;
    xxS(Count) = Interpolate(1, -1, ffElem(3), ffElem(4));
    yyS(Count) = Interpolate(1,  1, ffElem(3), ffElem(4));
end 
if SignVect(4) ~= SignVect(1)
    Count = Count + 1; 
    CutSegm(Count) = 4;
    xxS(Count) = Interpolate(-1, -1, ffElem(4), ffElem(1));
    yyS(Count) = Interpolate( 1, -1, ffElem(4), ffElem(1));
end 

if Count ~= 2
    error('Internal error.')
end

% Construct the six integration triangles ...
% ... number the nodes.
Count = 1;
for i = 1 : 4
    HelpNodes(i+Count-1) = i;
    if Count <= 2
        if CutSegm(Count) == i
            HelpNodes(i+Count) = 4+Count;
            Count = Count + 1;
        end
    end
end

% ... set node coordinates.
[xxTipRef, yyTipRef, CaseIntPointInElem] = ...
    ProjectRealToRefElem(xxElem, yyElem, xxTip, yyTip);
if CaseIntPointInElem == 0
    error('Crack tip is not in supposed crack-tip element!')
end
xxTri = [[-1  1 1 -1], xxS(1), xxS(2), xxTipRef];
yyTri = [[-1 -1 1  1], yyS(1), yyS(2), yyTipRef];

% ... construct the triangular integration mesh.
for i = 1 : 6
    MeshTri(i, 1) = 7;
    MeshTri(i, 2) = HelpNodes(i);
    if i == 6
        MeshTri(i, 3) = HelpNodes(1);
    else
        MeshTri(i, 3) = HelpNodes(i+1);
    end
end

% Set integration points in each triangle according to almost polar integration.
xxIntRef = zeros(1, 6*nQ1D*nQ1D);
yyIntRef = zeros(1, 6*nQ1D*nQ1D);
wwIntRef = zeros(1, 6*nQ1D*nQ1D);
[xxIntRefTri, yyIntRefTri, wwIntRefTri] = IntPoints2DRefElemTri(nQ1D);
for i = 1 : 6
    NodesTri = MeshTri(i, :);
    xxElemTri = xxTri(NodesTri);
    yyElemTri = yyTri(NodesTri);
    [xxInt, yyInt, wwInt] = IntPoints2DRealElemTri(xxElemTri, yyElemTri, ...
        xxIntRefTri, yyIntRefTri, wwIntRefTri, nQ1D*nQ1D);
    nQ2D = nQ1D * nQ1D;
    xxIntRef((i-1)*nQ2D+1 : i*nQ2D) = xxInt;
    yyIntRef((i-1)*nQ2D+1 : i*nQ2D) = yyInt;
    wwIntRef((i-1)*nQ2D+1 : i*nQ2D) = wwInt;
end

% % Plot situation.
% reset(cla), reset(clf)
% subplot(1,2,1)
% hold on
% patch(xxElem, yyElem, 'y')
% plot(xxTip, yyTip, 'b*')
% NN = 0.25*[(1-xxIntRef).*(1-yyIntRef); (1+xxIntRef).*(1-yyIntRef);  ...
%     (1+xxIntRef).*(1+yyIntRef); (1-xxIntRef).*(1+yyIntRef)];
% plot(xxElem * NN, yyElem * NN, 'k*')
% axis equal
% axis([min(xxElem)-0.1 max(xxElem)+0.1 min(yyElem)-0.1 max(yyElem)+0.1])
% 
% subplot(1,2,2)
% hold on
% patch([-1 1 1 -1], [-1 -1 1 1], 'y')
% for i = 1 : 6
%     NodesTri = MeshTri(i, :);
%     patch(xxTri(NodesTri), yyTri(NodesTri), 'y')
% end
% plot(xxS(1), yyS(1), 'b*')
% plot(xxS(2), yyS(2), 'b*')
% plot(xxTipRef, yyTipRef, 'b*')
% plot(xxIntRef, yyIntRef, 'k*')
% axis equal
% axis([-1.1 1.1 -1.1 1.1])

function [xStar] = Interpolate(x1, x2, f1, f2)

    xStar = x1 + (x2-x1) * f1 / (f1-f2);
 
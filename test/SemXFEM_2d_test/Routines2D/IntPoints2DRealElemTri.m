% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [xxInt, yyInt, wwInt] = IntPoints2DRealElemTri(xxElem, yyElem, ...
    xxIntRef, yyIntRef, wwIntRef, nQ)

% Get nQ integration points and weights in a real tri-element with  
% element nodes at (xxElem, yyElem).
%
% Example call:
% nQ = 11;
% [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DRefElemTri(nQ);
% [xxInt, yyInt, wwInt] = IntPoints2DRealElemTri([4 -5 -3], ...
%     [3 4 -4], xxIntRef, yyIntRef, wwIntRef, nQ*nQ)

xxInt = (1-xxIntRef-yyIntRef)*xxElem(1) + ...
    (xxIntRef)*xxElem(2) + (yyIntRef)*xxElem(3);
yyInt = (1-xxIntRef-yyIntRef)*yyElem(1) + ...
    (xxIntRef)*yyElem(2) + (yyIntRef)*yyElem(3);
wwInt = wwIntRef * (...
    (xxElem(2)-xxElem(1)) * (yyElem(3)-yyElem(1)) - ...
    (yyElem(2)-yyElem(1)) * (xxElem(3)-xxElem(1)));

% % Plot situation.
% reset(cla), reset(clf), hold on
% patch(xxElem, yyElem, 'y')
% a = plot(xxInt, yyInt, 'ko');
% set(a, 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% axis equal

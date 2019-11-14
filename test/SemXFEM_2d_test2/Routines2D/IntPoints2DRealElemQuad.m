% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [xxInt, yyInt, wwInt] = IntPoints2DRealElemQuad(xxElem, yyElem, ...
    xxIntRef, yyIntRef, wwIntRef, nQ)

% Get nQ integration points and weights in a real quad-element with  
% element nodes at (xxElem, yyElem).
%
% Example call:
% nQ = 11;
% [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DRefElemQuad(nQ);
% [xxInt, yyInt, wwInt] = IntPoints2DRealElemQuad([4 -5 -3 2], ...
%     [3 4 -4 -3], xxIntRef, yyIntRef, wwIntRef, nQ*nQ)

[N, dNdx, dNdy, xxInt, yyInt, wwInt] = ShapeFctsStrdFEM(...
    xxElem, yyElem, xxIntRef, yyIntRef, wwIntRef, nQ);

% % Plot situation.
% reset(cla), reset(clf), hold on
% patch(xxElem, yyElem, 'y')
% a = plot(xxInt, yyInt, 'ko');
% set(a, 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% axis equal

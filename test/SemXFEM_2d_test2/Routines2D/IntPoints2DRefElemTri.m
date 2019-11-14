% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DRefElemTri(nQ)

% Get (nQ*nQ) integration points and weights in a triangular
% reference element.

[xxQuad, yyQuad, wwQuad] = IntPoints2DRefElemQuad(nQ);

xxIntRef = 0.25*(1+xxQuad).*(1+yyQuad);
yyIntRef = 0.25*(1-xxQuad).*(1+yyQuad);
wwIntRef = 0.125*(1+yyQuad).*wwQuad;

% % Plot situation.
% reset(cla), reset(clf), hold on
% patch([0 1 0], [0 0 1], 'y')
% a = plot(xxIntRef, yyIntRef, 'ko');
% set(a, 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% axis equal
% axis([0 1 0 1])

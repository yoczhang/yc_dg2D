% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DRefElemQuad(nQ)

% Get (nQ*nQ) integration points and weights in a quadrilateral
% reference element.

[xx1d, ww1d] = IntPoints1DGauss(nQ);

% Project 1D integration points to 2D.
[XX, YY] = meshgrid(xx1d, xx1d);
[W1, W2] = meshgrid(ww1d, ww1d);
WW = W1 .* W2;

xxIntRef = XX(:)';
yyIntRef = YY(:)';
wwIntRef = WW(:)';

% % Plot situation.
% reset(cla), reset(clf), hold on
% patch([-1 1 1 -1], [-1 -1 1 1], 'y')
% a = plot(xxIntRef, yyIntRef, 'ko');
% set(a, 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% axis equal
% axis([-1 1 -1 1])

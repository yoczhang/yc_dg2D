% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [ff] = GetLevelSet(a, b, c, xx, yy)

% Get the level-set function such that the discontinuity is 
% defined by the linear function yy = -(a/b)*xx - c/b.
%
% Some exmples:
% a = 1; b =-1; c = 0.5; => y = x + 0.5
% a = 0; b = 1; c =-0.5; => horizontal line at y = 0.5
% a = 1; b = 0; c =-0.5; => vertical line at x = 0.5

ff = (a*xx + b*yy + c) / (sqrt(a^2 + b^2));

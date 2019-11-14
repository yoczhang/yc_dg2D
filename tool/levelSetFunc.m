function ff = levelSetFunc(aa,bb,cc,node)
%
%   %-----------------------------------------------------
%       Just copy form GetLevelSet.m (from SemXFEM_2d)
%   %-----------------------------------------------------
%
%
%   Get the level-set function such that the discontinuity is 
%   defined by the linear function yy = -(a/b)*xx - c/b.
%
%   Some exmples:
%   a = 1; b =-1; c = 0.5; => y = x + 0.5
%   a = 0; b = 1; c =-0.5; => horizontal line at y = 0.5
%   a = 1; b = 0; c =-0.5; => vertical line at x = 0.5
%
%
%   YcZhang 24/9/2017
%
%   Last modified 24/9/2017
%

xx = node(:,1);
yy = node(:,2);

ff = (aa*xx + bb*yy + cc) / (sqrt(aa^2 + bb^2));


end 
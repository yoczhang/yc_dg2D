% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [u1, du1dx, du1dy, u2, du2dx, du2dy, u3, du3dx, du3dy, u4, du4dx, du4dy] = ...
    EvalCrackTipEnr(CrackTip, xxInt, yyInt)

% Evaluate the branch enrichment functions for the crack-tip enrichment.

xx = xxInt - CrackTip.xx;
yy = yyInt - CrackTip.yy;

% Get polar coordinates (rr, th) out of (xx, yy).
rr = sqrt(xx*xx+yy*yy);

drrdx = xx/rr;
drrdy = yy/rr;

dthdx = -yy/(xx*xx+yy*yy);
if xx == 0
    dthdy = 0;
else
    dthdy = 1/(xx+yy*yy/xx);
end

if xx == 0 
    if yy == 0
        [xxInt yyInt]
        error('Integration point directly at crack tip!')
        %th = 0;
    elseif yy > 0
        th = 0.5*pi;
    else
        th = 1.5*pi;
    end
elseif yy == 0
    if xx > 0
        th = 0;
    else
        th = pi;
    end
elseif xx > 0 & yy > 0
    th = atan(yy/xx);
elseif xx < 0 & yy > 0
    th = pi + atan(yy/xx);
elseif xx < 0 & yy < 0
    th = pi + atan(yy/xx);
elseif xx > 0 & yy < 0
    th = 2*pi + atan(yy/xx);
else
    error('Internal error!')
end

thStar = (th - CrackTip.th);
% thStar must be in the range from -pi to +pi.
% The following mod. is important due to the multiplication with 1/2 below!!!
if abs(thStar) > pi 
    if thStar > 0
        thStar = thStar - 2*pi;
    else
        thStar = thStar + 2*pi;
    end
end

if abs(thStar) > pi 
    error('Internal error.')
end

% Define the four branch functions and derivatives.
u1 = sqrt(rr)*sin(0.5*thStar);
du1drr = 1/(2*sqrt(rr))*sin(0.5*thStar);
du1dth = 0.5*sqrt(rr)*cos(0.5*thStar);
du1dx = du1drr * drrdx + du1dth * dthdx;
du1dy = du1drr * drrdy + du1dth * dthdy;

u2 = sqrt(rr)*sin(0.5*thStar)*sin(thStar);
du2drr = 1/(2*sqrt(rr))*sin(0.5*thStar)*sin(thStar);
du2dth = sqrt(rr)*(0.5*cos(0.5*thStar)*sin(thStar)+sin(0.5*thStar)*cos(thStar));
du2dx = du2drr * drrdx + du2dth * dthdx;
du2dy = du2drr * drrdy + du2dth * dthdy;

u3 = sqrt(rr)*cos(0.5*thStar);
du3drr = 1/(2*sqrt(rr))*cos(0.5*thStar);
du3dth = -0.5*sqrt(rr)*sin(0.5*thStar);
du3dx = du3drr * drrdx + du3dth * dthdx;
du3dy = du3drr * drrdy + du3dth * dthdy;

u4 = sqrt(rr)*cos(0.5*thStar)*sin(thStar);
du4drr = 1/(2*sqrt(rr))*cos(0.5*thStar)*sin(thStar);
du4dth = sqrt(rr)*(-0.5*sin(0.5*thStar)*sin(thStar)+cos(0.5*thStar)*cos(thStar));
du4dx = du4drr * drrdx + du4dth * dthdx;
du4dy = du4drr * drrdy + du4dth * dthdy;

function d = DisPtToLine( P, PA, PB) 
%
%   Give the distence form P to line [PA,PB].
%   - PA: (x,y)-coordinate at point A,
%   - PB: (x,y)-coordinate at point B.
%   
%   The examples: compute the distance from (0,0) to line1,
%   line1x = [0.3, 1]';
%   line1y = 0.1*line1x + 1.25;
%   P1A = [line1x(1), line1y(1)];
%   P1B = [line1x(2), line1y(2)];
%   dis1 = DisPtToLine([0,0],P1A,P1B);
%

Px = P(:,1);
Py = P(:,2);
PAx = PA(:,1);
PAy = PA(:,2);
PBx = PB(:,1);
PBy = PB(:,2);

d=abs((PAx-Px).*(PBy-Py)-(PBx-Px).*(PAy-Py))/sqrt((PAx-PBx)^2+(PAy-PBy)^2);

end
% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [xxIntRef, yyIntRef, wwInt, nxVect, nyVect, xx1, yy1, xx2, yy2, ...
    nxRef, nyRef] = IntPoints2DAlongLevelSet(xxElem, yyElem, ff, nQ)

% Set integration points in the 2D REFERENCE element along the 1D 
% discontinuity. Compute the integration weights for the REAL element.
% Consequently, this routine works for cut elements only.
% At each integration point, the normal vector in the REAL element
% is computed (nxVect, nyVect). The normal vector in the REFERENCE element
% (nxRef, nyRef) is assumed to be constant for all integration points.
% The coordinates (xx1, yy1) and (xx2, yy2) are the intersection points
% of the discontinuity with the REAL element edges.
% The input are the element nodes at (xxElem, yyElem) and the level-set 
% function at the element nodes ff. Note that nQ refers to the number of 
% integration points in the 1D-reference element.

% Example call:
% xxElem = [-1 1 1 -1];
% yyElem = [-1 -1 1 1];
% ff = [-1 1 -1 -1]/2;
% [xxIntRef, yyIntRef, wwInt, nxVect, nyVect, xx1, yy1, xx2, yy2, ...
%     nxRef, nyRef] = IntPoints2DAlongLevelSet(xxElem, yyElem, ff, 21);

SignVect = sign(ff);
if isempty(find(SignVect==0)) == 0
    error('Level set function is exactly zero at a node.')
end

if SignVect(1)==SignVect(2) & SignVect(1)==SignVect(3) & SignVect(1)==SignVect(4)
    error('Level set function does not cut through element, no int. points must be set.')
end

% Find the two intersection-points of the reference element with the level set function.
r = 10*ones(4, 1); s = 10*ones(4, 1); % (this value for r and s can never be reached, use 10 as marker)
if SignVect(1)~=SignVect(2) 
    r(1) = -1+2*ff(1)/(ff(1)-ff(2));
    s(1) = -1;
end
if SignVect(2)~=SignVect(3)
    r(2) = 1;
    s(2) = -1+2*ff(2)/(ff(2)-ff(3));
end
if SignVect(3)~=SignVect(4)
    r(3) = -1+2*ff(4)/(ff(4)-ff(3));
    s(3) = 1;
end
if SignVect(4)~=SignVect(1)
    r(4) = -1;
    s(4) = -1+2*ff(1)/(ff(1)-ff(4));
end

Pos = find(r~=10);
if length(Pos)~=2
    if length(Pos) == 4
        % This is a +-+- or -+-+ case.
        error('Special situation, disc. can not be uniquely determined from level set.')
    else
        error('Internal error.')
    end
end

ra = r(Pos(1)); sa = s(Pos(1));
rb = r(Pos(2)); sb = s(Pos(2));

% Compute normal vector in reference element (assumed to be constant).
dr = ra - rb; 
ds = sa - sb; 
ll = sqrt(dr*dr+ds*ds); 
nxRef = -ds/ll; nyRef = dr/ll;

% Compute position of intersection in real domain.
[xx1, yy1, wwIntDummy] = IntPoints2DRealElemQuad(xxElem, yyElem, ra, sa, 1, 1);
[xx2, yy2, wwIntDummy] = IntPoints2DRealElemQuad(xxElem, yyElem, rb, sb, 1, 1);

% Distribute integration points along the level set function in the reference 
% element. Furthermore, set helping points for the numerical evaluation of 
% the modification-factors of the integration weights.
[zz, ww] = IntPoints1DGauss(nQ);
xxIntRef = 0.5*(zz+1)*(rb-ra)+ra;
yyIntRef = 0.5*(zz+1)*(sb-sa)+sa;
wwIntRef = ww * sqrt((rb-ra)^2+(sb-sa)^2) / 2;

EPS = 1.e-6;
xxIntRefIncr = xxIntRef + EPS*(rb-ra);
yyIntRefIncr = yyIntRef + EPS*(sb-sa);

dxRefVect = xxIntRefIncr - xxIntRef;
dyRefVect = yyIntRefIncr - yyIntRef;

% Project all the points in the refernce element into the real element.
[xxInt, yyInt, wwIntDummy] = IntPoints2DRealElemQuad(xxElem, yyElem, xxIntRef, yyIntRef, ww, nQ);
[xxIntIncr, yyIntIncr, wwIntDummy] = IntPoints2DRealElemQuad(xxElem, yyElem, xxIntRefIncr, yyIntRefIncr, ww, nQ);

dxVect = xxIntIncr - xxInt;
dyVect = yyIntIncr - yyInt;

% Compute the modification-factors for the integration weights.
wwIntModFactors = sqrt(dxVect.^2+dyVect.^2) ./ sqrt(dxRefVect.^2+dyRefVect.^2);
wwInt = wwIntRef .* wwIntModFactors;

% Compute normal vectors in the real domain.
llVect = sqrt(dxVect.*dxVect+dyVect.*dyVect); 
nxVect = -dyVect./llVect; 
nyVect =  dxVect./llVect;

% The normal vector always points in the direction of the positive level set area.
if Pos(1)==1 & Pos(2)==2
    if SignVect(2)==1
        nxVect = -nxVect;
        nyVect = -nyVect;
    end
elseif Pos(1)==1 & Pos(2)==3
    if SignVect(2)==1
        nxVect = -nxVect;
        nyVect = -nyVect;
    end
elseif Pos(1)==1 & Pos(2)==4
    if SignVect(2)==1
        nxVect = -nxVect;
        nyVect = -nyVect;
    end
elseif Pos(1)==2 & Pos(2)==3
    if SignVect(2)==-1
        nxVect = -nxVect;
        nyVect = -nyVect;
    end
elseif Pos(1)==2 & Pos(2)==4
    if SignVect(2)==-1
        nxVect = -nxVect;
        nyVect = -nyVect;
    end
elseif Pos(1)==3 & Pos(2)==4
    if SignVect(2)==-1
        nxVect = -nxVect;
        nyVect = -nyVect;
    end
end

% % Plot situation.
% reset(cla), reset(clf)
% subplot(1, 2, 1), hold on
% patch([-1 1 1 -1], [-1 -1 1 1], 'y')
% plot(xxIntRef, yyIntRef, 'k*')
% plot(xxIntRefIncr, yyIntRefIncr, 'b*')
% for i = 1 : nQ
%     line([xxIntRef(i) xxIntRef(i)+nxRef*0.1], [yyIntRef(i) yyIntRef(i)+nyRef*0.1])
% end
% title('Reference element')
% axis equal
% axis([-1.1 1.1 -1.1 1.1])
% 
% subplot(1, 2, 2), hold on
% patch(xxElem, yyElem, 'y')
% plot(xxInt, yyInt, 'k*')
% plot(xxIntIncr, yyIntIncr, 'b*')
% for i = 1 : nQ
%     line([xxInt(i) xxInt(i)+nxVect(i)*0.1], [yyInt(i) yyInt(i)+nyVect(i)*0.1])
% end
% title('Real element')
% axis equal
% axis([min(xxElem) max(xxElem) min(yyElem) max(yyElem)])

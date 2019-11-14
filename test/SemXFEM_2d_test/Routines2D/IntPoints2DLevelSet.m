% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [xxIntTotal, yyIntTotal, wwIntTotal] = IntPoints2DLevelSet(ff, nQ);

% Set integration points and weigths in the 2D reference element, considering 
% special distributions for elements cut by the discontinuity. 
% Type 1-cases (++++, ----)
% Type 2-cases (-+++, +-++, ++-+, +++-, +---, -+--, --+-, ---+)
% Type 3-cases (+--+, ++--, +-+-, -++-, --++, -+-+)

% Example call:
% [xxIntTotal, yyIntTotal, wwIntTotal] = IntPoints2DLevelSet([-7  5  3  1], 5);
% [xxIntTotal, yyIntTotal, wwIntTotal] = IntPoints2DLevelSet([-7 -5  3  1], 5);
% [xxIntTotal, yyIntTotal, wwIntTotal] = IntPoints2DLevelSet([ 7  5  3  1], 5);

SignVect = sign(ff);
if isempty(find(SignVect==0)) == 0
    error('Level set function is exactly zero at a node.')
end

% Get integration points and weights in the reference element.
[xxIntRef, yyIntRef, wwIntRef] = IntPoints2DRefElemQuad(nQ);
nQnQ = nQ * nQ;

if SignVect(1)==SignVect(2) & SignVect(1)==SignVect(3) & SignVect(1)==SignVect(4)
    % Treatment of Type-1-elements (no disc., i.e. all level set values have the same sign).
    % (++++, ----)

    xxIntTotal = [xxIntRef];
    yyIntTotal = [yyIntRef];
    wwIntTotal = [wwIntRef];

    IntSide = SignVect(1) * ones(nQnQ, 1);

elseif SignVect(1)*SignVect(2)*SignVect(3)*SignVect(4) < 0
    % Treatment of Type-2-elements (one level function value is on the other side than the other three).
    % (-+++, +-++, ++-+, +++-, +---, -+--, --+-, ---+)
 
    % Find the one node on the other side.
    if SignVect(1)~=SignVect(2) & SignVect(1)~=SignVect(3) & SignVect(1)~=SignVect(4)
        Pos = 1;
        xA = -1+2*ff(4)/(ff(4)-ff(1));
        yB = -1+2*ff(1)/(ff(1)-ff(2));
        IntSide = -SignVect(1) * [ones(nQnQ, 1); ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1)];
    elseif SignVect(2)~=SignVect(1) & SignVect(2)~=SignVect(3) & SignVect(2)~=SignVect(4)
        Pos = 2;
        xA = -1+2*ff(1)/(ff(1)-ff(2));
        yB = -1+2*ff(2)/(ff(2)-ff(3));
        IntSide = -SignVect(2) * [ones(nQnQ, 1); ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1)];
    elseif SignVect(3)~=SignVect(1) & SignVect(3)~=SignVect(2) & SignVect(3)~=SignVect(4)
        Pos = 3;
        xA = -1+2*ff(2)/(ff(2)-ff(3));
        yB = -1+2*ff(3)/(ff(3)-ff(4));
        IntSide = -SignVect(3) * [ones(nQnQ, 1); ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1)];
    elseif SignVect(4)~=SignVect(1) & SignVect(4)~=SignVect(2) & SignVect(4)~=SignVect(3)
        Pos = 4;
        xA = -1+2*ff(3)/(ff(3)-ff(4));
        yB = -1+2*ff(4)/(ff(4)-ff(1));
        IntSide = -SignVect(4) * [ones(nQnQ, 1); ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1)];
    else
        error('Internal error.')
    end

    x1 =-1; y1 =-1;
    x2 = 1; y2 =-1;
    x3 = 1; y3 = 1;
    x4 =-1; y4 = 1;

    % Make sub-elements and int. points as if the SECOND node were on the other side (Pos=2).
    yA = -1;
    xB = 1;
    xC = 0.5*(xA+xB);
    yC = 0.5*(yA+yB);
    xD = 1/3+2/3*xC;
    yD = -1/3+2/3*yC;
    xE = 0.5*(1+xA);
    yE = -1;
    xF = 1;
    yF = 0.5*(-1+yB);

    % Project integration points into sub-elements.
    [xxInt1, yyInt1, wwInt1] = IntPoints2DRealElemQuad([x1 xA xC x4], [y1 yA yC y4], xxIntRef, yyIntRef, wwIntRef, nQnQ);
    [xxInt2, yyInt2, wwInt2] = IntPoints2DRealElemQuad([x4 xC xB x3], [y4 yC yB y3], xxIntRef, yyIntRef, wwIntRef, nQnQ);
    [xxInt3, yyInt3, wwInt3] = IntPoints2DRealElemQuad([xA xE xD xC], [yA yE yD yC], xxIntRef, yyIntRef, wwIntRef, nQnQ);
    [xxInt4, yyInt4, wwInt4] = IntPoints2DRealElemQuad([xB xC xD xF], [yB yC yD yF], xxIntRef, yyIntRef, wwIntRef, nQnQ);
    [xxInt5, yyInt5, wwInt5] = IntPoints2DRealElemQuad([x2 xF xD xE], [y2 yF yD yE], xxIntRef, yyIntRef, wwIntRef, nQnQ);

    xxIntTotal = [xxInt1, xxInt2, xxInt3, xxInt4, xxInt5];
    yyIntTotal = [yyInt1, yyInt2, yyInt3, yyInt4, yyInt5];
    wwIntTotal = [wwInt1, wwInt2, wwInt3, wwInt4, wwInt5];

%     patch([x1 xA xC x4], [y1 yA yC y4], 'b')
%     patch([x4 xC xB x3], [y4 yC yB y3], 'b')
%     patch([xA xE xD xC], [yA yE yD yC], 'b')
%     patch([xB xC xD xF], [yB yC yD yF], 'b')
%     patch([x2 xF xD xE], [y2 yF yD yE], 'b')
% 
%     plot(xxInt1, yyInt1, 'k*')
%     plot(xxInt2, yyInt2, 'k*')
%     plot(xxInt3, yyInt3, 'k*')
%     plot(xxInt4, yyInt4, 'k*')
%     plot(xxInt5, yyInt5, 'k*')

    % Rotate integration points according to which of the nodes is really on the other side.
    if Pos == 1
        T = [0 1; -1 0]; % Rotation -90 degree.
    elseif Pos == 2
        T = [1 0; 0 1]; % Nothing to do.
    elseif Pos == 3
        T = [0 -1; 1 0]; % Rotation 90 degree.
    elseif Pos == 4
        T = [-1 0; 0 -1]; % Rotation 180 degree.
    end

    HelpMat = T * [xxIntTotal; yyIntTotal];
    xxIntTotal = HelpMat(1,:);
    yyIntTotal = HelpMat(2,:);

else
    % Treatment of Type-3-elements (two level function values are on the other side than the other two).
    % (+--+, ++--, +-+-, -++-, --++, -+-+)

    if SignVect(1)==SignVect(2) & SignVect(3)==SignVect(4)
        Pos = 1;
        xA = -1+2*ff(4)/(ff(4)-ff(1));
        xB = -1+2*ff(3)/(ff(3)-ff(2));
        IntSide = SignVect(1) * [-ones(nQnQ, 1); ones(nQnQ, 1)];
    elseif SignVect(1)==SignVect(4) & SignVect(2)==SignVect(3)
        Pos = 2;
        xA = -1+2*ff(1)/(ff(1)-ff(2));
        xB = -1+2*ff(4)/(ff(4)-ff(3));
        IntSide = SignVect(1) * [ones(nQnQ, 1); -ones(nQnQ, 1)];
    elseif SignVect(1)==SignVect(3) & SignVect(2)==SignVect(4)
        error('Special situation, disc. can not be uniquely determined from level set.')
    else
        error('Internal error.')
    end

    x1 =-1; y1 =-1;
    x2 = 1; y2 =-1;
    x3 = 1; y3 = 1;
    x4 =-1; y4 = 1;

    % Make sub-elements and int. points as if the SECOND node were on the other side (Pos=2).
    yA =-1; yB = 1;

    % Project integration points into sub-elements.
    [xxInt1, yyInt1, wwInt1] = IntPoints2DRealElemQuad([x1 xA xB x4], [y1 yA yB y4], xxIntRef, yyIntRef, wwIntRef, nQnQ);
    [xxInt2, yyInt2, wwInt2] = IntPoints2DRealElemQuad([xB xA x2 x3], [yB yA y2 y3], xxIntRef, yyIntRef, wwIntRef, nQnQ);

    xxIntTotal = [xxInt1, xxInt2];
    yyIntTotal = [yyInt1, yyInt2];
    wwIntTotal = [wwInt1, wwInt2];

%     patch([x1 xA xB x4], [y1 yA yB y4], 'b')
%     patch([xB xA x2 x3], [yB yA y2 y3], 'b')
% 
%     plot(xxInt1, yyInt1, 'k*')
%     plot(xxInt2, yyInt2, 'k*')

    % Rotate integration points according to which of the nodes is really on the other side.
    if Pos == 1
        T = [0 1; -1 0]; % Rotation -90 degree.
    elseif Pos == 2
        T = [1 0; 0 1]; % Nothing to do.
    end

    HelpMat = T * [xxIntTotal; yyIntTotal];
    xxIntTotal = HelpMat(1,:);
    yyIntTotal = HelpMat(2,:);

end


% % Plot situation.
% reset(cla); reset(clf); hold on
% if exist('xA')==1
%     Cut1 = T * [xA; yA]; Cut2 = T * [xB; yB];
%     a = line([Cut1(1) Cut2(1)], [Cut1(2) Cut2(2)]);
%     set(a, 'LineWidth', 3, 'Color', 'b')
% end
% patch([-1 1 1 -1], [-1 -1 1 1], -0.01*[1 1 1 1], 'y');
% a = plot(xxIntTotal, yyIntTotal, 'ko');
% set(a, 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% axis equal
% axis([-1 1 -1 1])

% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function PlotLevelSet(Mesh, xx, yy, ff, ElemNum)

% Plot the discontinuity defined by the level-set function. Assume a
% straight line in the real element.

if isempty(find(sign(ff)==0)) == 0
    %error('Level-set function is exactly zero at a node!')
    disp('Level-set function is exactly zero at a node!')
end

for CurrElem = 1 : ElemNum

    Nodes = Mesh(CurrElem, :);

    xxElem = xx(Nodes);
    yyElem = yy(Nodes);
    ffElem = ff(Nodes);
    SignVect = sign(ffElem);

    if min(SignVect) == max(SignVect) % Element is not cut by disc.
        continue
    end

    % Find the element edges which are cut and compute intersection point.
    Count = 0;
    if SignVect(1) ~= SignVect(2)
        Count = Count + 1;
        xxS(Count) = Interpolate(xxElem(1), xxElem(2), ffElem(1), ffElem(2));
        yyS(Count) = Interpolate(yyElem(1), yyElem(2), ffElem(1), ffElem(2));
    end 
    if SignVect(2) ~= SignVect(3)
        Count = Count + 1;
        xxS(Count) = Interpolate(xxElem(2), xxElem(3), ffElem(2), ffElem(3));
        yyS(Count) = Interpolate(yyElem(2), yyElem(3), ffElem(2), ffElem(3));
    end 
    if SignVect(3) ~= SignVect(4)
        Count = Count + 1;
        xxS(Count) = Interpolate(xxElem(3), xxElem(4), ffElem(3), ffElem(4));
        yyS(Count) = Interpolate(yyElem(3), yyElem(4), ffElem(3), ffElem(4));
    end 
    if SignVect(4) ~= SignVect(1)
        Count = Count + 1;
        xxS(Count) = Interpolate(xxElem(4), xxElem(1), ffElem(4), ffElem(1));
        yyS(Count) = Interpolate(yyElem(4), yyElem(1), ffElem(4), ffElem(1));
    end 

    if Count == 2
        line([xxS(1) xxS(2)], [yyS(1) yyS(2)], -0.001*[1 1])
    else
        error('Internal error.')
    end

end

function [xStar] = Interpolate(x1, x2, f1, f2)

    xStar = x1 + (x2-x1) * f1 / (f1-f2);

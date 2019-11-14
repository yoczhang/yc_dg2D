% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function PlotLevelSet_yctest(Mesh, xx, yy, ffx, ffy, ElemNum)

% Plot the discontinuity defined by the level-set function. Assume a
% straight line in the real element.

% if isempty(find(sign(ff)==0)) == 0
%     error('Level-set function is exactly zero at a node!')
%     
% end

%------------------------------------
if isempty(find(sign(ffy)==0,1)) == 0
    disp('Level-set function is exactly zero at some nodes')
end 

m = @(x,singleNE) mod(x,singleNE)+(x==singleNE)*singleNE;
%------------------------------------

for CurrElem = 1 : ElemNum
    if iscell(Mesh)
        Nodes = Mesh{CurrElem};
    else
        Nodes = Mesh(CurrElem, :);
    end

    xxElem = xx(Nodes);
    yyElem = yy(Nodes);
    ffxElem = ffx(Nodes);
    ffyElem = ffy(Nodes);
    SignVect_x = sign(ffxElem);
    SignVect_y = sign(ffyElem);
    

%     if min(SignVect_y) == max(SignVect_y) % Element is not cut by disc.
%         continue
%     end
    
    if ( min(SignVect_y) == max(SignVect_y) ) 
        continue
    end
    
    
    if (max(SignVect_x) >= 0) && (min(SignVect_x) == max(SignVect_x))
        continue
    end
    
    % Find the element edges which are cut and compute intersection point.
    %-------------------------------
    xxS = zeros(2,1);
    yyS = zeros(2,1);
    singleNE = length(Nodes);
    Count = 0;
    for nn = 1:singleNE
        if SignVect_y(m(nn,singleNE)) ~= SignVect_y(m(nn+1,singleNE))
            Count = Count + 1;
            xxS(Count) = Interpolate(xxElem(m(nn,singleNE)), xxElem(m(nn+1,singleNE)), ...
                ffyElem(m(nn,singleNE)), ffyElem(m(nn+1,singleNE)));
            yyS(Count) = Interpolate(yyElem(m(nn,singleNE)), yyElem(m(nn+1,singleNE)), ...
                ffyElem(m(nn,singleNE)), ffyElem(m(nn+1,singleNE)));
            if xxS(Count) > 0.52
                xxS(Count) = 0.52;
            end
        end 
    end % for nn
    %-------------------------------
    
%     Count = 0;
%     if SignVect(1) ~= SignVect(2)
%         Count = Count + 1;
%         xxS(Count) = Interpolate(xxElem(1), xxElem(2), ffElem(1), ffElem(2));
%         yyS(Count) = Interpolate(yyElem(1), yyElem(2), ffElem(1), ffElem(2));
%     end 
%     if SignVect(2) ~= SignVect(3)
%         Count = Count + 1;
%         xxS(Count) = Interpolate(xxElem(2), xxElem(3), ffElem(2), ffElem(3));
%         yyS(Count) = Interpolate(yyElem(2), yyElem(3), ffElem(2), ffElem(3));
%     end 
%     if SignVect(3) ~= SignVect(4)
%         Count = Count + 1;
%         xxS(Count) = Interpolate(xxElem(3), xxElem(4), ffElem(3), ffElem(4));
%         yyS(Count) = Interpolate(yyElem(3), yyElem(4), ffElem(3), ffElem(4));
%     end 
%     if SignVect(4) ~= SignVect(1)
%         Count = Count + 1;
%         xxS(Count) = Interpolate(xxElem(4), xxElem(1), ffElem(4), ffElem(1));
%         yyS(Count) = Interpolate(yyElem(4), yyElem(1), ffElem(4), ffElem(1));
%     end 
    if Count == 0 
        continue
    elseif Count == 2
        line([xxS(1) xxS(2)], [yyS(1) yyS(2)], -0.001*[1 1])
    else
        [yyS1,ii1,jj1]=unique(yyS,'rows','stable');
        xxS1=xxS(ii1);
        line([xxS1(1) xxS1(2)], [yyS1(1) yyS1(2)], -0.001*[1 1]);
        %error('Internal error.')
    end

end
end % function

function [xStar] = Interpolate(x1, x2, f1, f2)

    xStar = x1 + (x2-x1) * f1 / (f1-f2);
end

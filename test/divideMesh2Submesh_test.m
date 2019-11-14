function [cutElem]=divideMesh2Submesh_test(node, elem)
%
%
%   We let Nnodes denote the number of nodes of Th,
%               Nelems denote the number of the elements of Th.
%
%   input:
%       node, [Nnodes x 2], the (x-coord, y-coord) of nodes of Th.
%       interfaceFunc, the 
%
%
%   YcZhang 15/9/2017
%
%   Last modified 23/9/2017
%

%% one build-in func
m = @(x,singleNE) mod(x,singleNE)+(x==singleNE)*singleNE;

%% Get level-set function.
aa=1; bb=-1; cc=0.5; % Vertical disc. at y=0.5.
ff = GetLevelSet(aa, bb, cc, node);

%% get the mesh coord and domain coord
xx = node(:,1);
yy = node(:,2);

xx_domain = [0, 0.5];
yy_domain = [-1, 1];

%% judge is there has the level-func on mesh node.
if isempty(find(sign(ff)==0,1)) == 0
    disp('Level-set function is exactly zero at some nodes')
end 

%%
Nelems = size(elem, 1);
cutElem = zeros(Nelems,1);
cutElemCount = 0;

for CurrElem = 1 : Nelems
    if iscell(elem)
        CurrNodes = elem{CurrElem};
    else
        CurrNodes = elem(CurrElem, :);
    end
    
    xxElem = xx(CurrNodes);
    yyElem = yy(CurrNodes);
    ffElem = ff(CurrNodes);
    SignVect = sign(ffElem);
    
    if min(SignVect) == max(SignVect) % Element is not cut by disc.
        continue
    end
    
    
    %% Find the element edges which are cut and compute intersection point.
    xxS = zeros(2,1);
    yyS = zeros(2,1);
    singleNE = length(CurrNodes);
    Count = 0;
    cutElemIndx_xcoord = 0;
    cutElemIndx_ycoord = 0;
    %x_count = 0;
    %y_count = 0;
    
    for nn = 1:singleNE
        if SignVect(m(nn,singleNE)) ~= SignVect(m(nn+1,singleNE))
            %% get the coord
            Count = Count + 1;
            xxS(Count) = Interpolate(xxElem(m(nn,singleNE)), xxElem(m(nn+1,singleNE)), ...
                ffElem(m(nn,singleNE)), ffElem(m(nn+1,singleNE)));
            yyS(Count) = Interpolate(yyElem(m(nn,singleNE)), yyElem(m(nn+1,singleNE)), ...
                ffElem(m(nn,singleNE)), ffElem(m(nn+1,singleNE)));
            
            %% amend the coord
            if xxS(Count) > max(xx_domain)
                xxS(Count) = max(xx_domain);
                
                cutElemIndx_xcoord = cutElemIndx_xcoord + 1;
                
                % here we also need to get the yyS according the level func.
                if bb~=0
                    yyS(Count) = -(aa/bb)*xxS(Count) - cc/bb;
                end
                %x_count = x_count +1;
            elseif xxS(Count) < min(xx_domain)
                xxS(Count) = min(xx_domain);
                
                cutElemIndx_xcoord = cutElemIndx_xcoord + 1;
                
                % here we also need to get the yyS according the level func.
                if bb~=0
                    yyS(Count) = -(aa/bb)*xxS(Count) - cc/bb;
                end
                %x_count = x_count +1;
            end
            
            if yyS(Count) > max(yy_domain)
                yyS(Count) = max(yy_domain);
                
                cutElemIndx_ycoord = cutElemIndx_ycoord + 1;
                
                % here we also need to get the xxS according the level func.
                if aa~=0
                    xxS(Count) = -(bb/aa)*yyS(Count) - cc/aa;
                end
            elseif yyS(Count) < min(yy_domain)
                yyS(Count) = min(yy_domain);
                
                cutElemIndx_ycoord = cutElemIndx_ycoord + 1;
                
                % here we also need to get the xxS according the level func.
                if aa~=0
                    xxS(Count) = -(bb/aa)*yyS(Count) - cc/aa;
                end
            end 
        end % if SignVect(m(nn,singleNE)) ~= SignVect(m(nn+1,singleNE))  
    end % for nn
    
    %% get the cut element index.
	if cutElemIndx_xcoord<2 && cutElemIndx_ycoord<2
        cutElemCount = cutElemCount + 1;
        cutElem(cutElemCount) = CurrElem;
	end % if
    
    %% plot the level-func
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
    
    
end % for CurrElem

cutElem = cutElem(cutElem>0);

end % function






%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%% sub function
function [ff] = GetLevelSet(a, b, c, node)

% Get the level-set function such that the discontinuity is 
% defined by the linear function yy = -(a/b)*xx - c/b.
%
% Some exmples:
% a = 1; b =-1; c = 0.5; => y = x + 0.5
% a = 0; b = 1; c =-0.5; => horizontal line at y = 0.5
% a = 1; b = 0; c =-0.5; => vertical line at x = 0.5

xx = node(:,1);
yy = node(:,2);

ff = (a*xx + b*yy + c) / (sqrt(a^2 + b^2));
end % function GetLevelSet


%% sub function
function [xStar] = Interpolate(x1, x2, f1, f2)

    xStar = x1 + (x2-x1) * f1 / (f1-f2);
end 

%--------------------------------------------------------------------------------------------------------------------%
%----------- end sub function ----------------------------------------------------------------------------------%

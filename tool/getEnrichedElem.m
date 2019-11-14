function [enrichedElem, cutNode] = getEnrichedElem(levelFunc, node, elem, plotLevelFunc)
%
%   [30/9/2017], add plotLevelFunc, one pars to jude whether to plot the
%                       level-func.
%
%   %-----------------------------------------------------
%       Just copy form divideMesh2Submesh_test.m (from test file)
%   %-----------------------------------------------------
%
%   input:
%       levelFunc, a structure-type, must contain the level-set func value: ff.
%       node, [Nnodes x 2], the (x-coord,y-coord) of mesh nodes.
%       elem, if elem is the cell-type, [Nelems x 1].
%                else [Nelems x singleNE].
%   output:
%       enrichedElem, [NenrichedElem x 1], the index of the enriched elems.
%       cutNode, [NcutNode x 3], the 1-th, 2-th column is the (x-coord, y-coord),
%           the 3-th column is the elemIdex which the cutNode belong to.
%
%
%   YcZhang 24/9/2017
%
%   Last modified 30/9/2017
%
%

%% ------ first we will need one build-in func
m = @(x,singleNE) mod(x,singleNE)+(x==singleNE)*singleNE;


%% ------ get the level-set func informations
ff = levelFunc.ff;
if isfield(levelFunc,'aa')
    aa = levelFunc.aa;
else
    aa = [];
end
if isfield(levelFunc,'bb')
    bb = levelFunc.bb;
else
    bb = [];
end
if isfield(levelFunc,'cc')
    cc = levelFunc.cc;
end
if isfield(levelFunc,'x_domain')
    x_domain = levelFunc.x_domain;
else
    x_domain = [];
end
if isfield(levelFunc,'y_domain')
    y_domain = levelFunc.y_domain;
else
    y_domain = [];
end


%% get the mesh coord and domain coord
xx = node(:,1);
yy = node(:,2);

%% judge is there has the level-func on mesh node.
if isempty(find(sign(ff)==0,1)) == 0
    disp(' ')
    disp(' ************************ warning ************************')
    disp(' in getEnrichedElem func:')
    disp('   Level-set function is exactly zero at some nodes')
    disp(' **********************************************************')
end 

%%
Nelems = size(elem, 1);
cutElem = zeros(Nelems,1);
cutElemCount = 0;
cutNode = [];

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
            if xxS(Count) > max(x_domain)
                xxS(Count) = max(x_domain);
                
                cutElemIndx_xcoord = cutElemIndx_xcoord + 1;
                
                % here we also need to get the yyS according the level func.
                if ~isempty(bb) && bb~=0
                    yyS(Count) = -(aa/bb)*xxS(Count) - cc/bb;
                end
                %x_count = x_count +1;
            elseif xxS(Count) < min(x_domain)
                xxS(Count) = min(x_domain);
                
                cutElemIndx_xcoord = cutElemIndx_xcoord + 1;
                
                % here we also need to get the yyS according the level func.
                if ~isempty(bb) && bb~=0
                    yyS(Count) = -(aa/bb)*xxS(Count) - cc/bb;
                end
                %x_count = x_count +1;
            end
            
            if yyS(Count) > max(y_domain)
                yyS(Count) = max(y_domain);
                
                cutElemIndx_ycoord = cutElemIndx_ycoord + 1;
                
                % here we also need to get the xxS according the level func.
                if ~isempty(aa) && aa~=0
                    xxS(Count) = -(bb/aa)*yyS(Count) - cc/aa;
                end
            elseif yyS(Count) < min(y_domain)
                yyS(Count) = min(y_domain);
                
                cutElemIndx_ycoord = cutElemIndx_ycoord + 1;
                
                % here we also need to get the xxS according the level func.
                if ~isempty(aa) && aa~=0
                    xxS(Count) = -(bb/aa)*yyS(Count) - cc/aa;
                end
            end 
        end % if SignVect(m(nn,singleNE)) ~= SignVect(m(nn+1,singleNE))  
    end % for nn
    
    
    %% plot the level-func
    if Count == 0 
        continue
    elseif Count == 2
        xxS1= xxS;
        yyS1 = yyS;
        if plotLevelFunc % the plotLevelFunc may take 1 or 0.
            line([xxS1(1) xxS1(2)], [yyS1(1) yyS1(2)], -0.001*[1 1])
        end
    else
        [yyS1,ii1,jj1]=unique(yyS,'rows','stable');
        xxS1=xxS(ii1);
        
        if plotLevelFunc
            if length(xxS1) == 1
                line(xxS1(1), yyS1(1), -0.001*1);
            elseif length(xxS1) == 2
                line([xxS1(1) xxS1(2)], [yyS1(1) yyS1(2)], -0.001*[1 1]);
            end
        end
        
        %error('Internal error.')
    end
    
    %% get the cut element index.
	if cutElemIndx_xcoord<2 && cutElemIndx_ycoord<2
        cutElemCount = cutElemCount + 1;
        cutElem(cutElemCount) = CurrElem;
        
        currt_cutnodes = [xxS1, yyS1, CurrElem*ones(size(xxS1))];
        cutNode_temp = [cutNode; currt_cutnodes];
        cutNode = cutNode_temp;
        
        %cutElemInfo.enrichedElemIndx = CurrElem;
        %cutNode = [xxS1,yyS1];
        %cutElemInfo.cutNode = cutNode;
	end % if
    
    
end % for CurrElem

enrichedElem = cutElem(cutElem>0);

end % function


%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%

%% sub function
function [xStar] = Interpolate(x1, x2, f1, f2)

    xStar = x1 + (x2-x1) * f1 / (f1-f2);
end 

%--------------------------------------------------------------------------------------------------------------------%
%----------- end sub function ----------------------------------------------------------------------------------%

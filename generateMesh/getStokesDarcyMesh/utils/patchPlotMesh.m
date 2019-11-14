function patchPlotMesh(node, elem)
%
%   %-----------------------------------------------------
%       Using the matlab func patch(..).
%       Just copy form PlotMesh.m (from SemXFEM_2d)
%   %-----------------------------------------------------
%
%
%
%   YcZhang 24/9/2017
%
%   Last modified 24/9/2017
%
%

% Plot the mesh.
xx = node(:,1);
yy = node(:,2);

ElemNum = size(elem,1);

figure
reset(cla), reset(clf), hold on

for CurrElem = 1 : ElemNum
    if iscell(elem)
        CurrNodes = elem{CurrElem};
    else
        CurrNodes = elem(CurrElem, :);
    end
    
    L = ones(1,length(CurrNodes));
    %patch(xx(CurrNodes), yy(CurrNodes), -0.005*L, 'w')
    patch(xx(CurrNodes), yy(CurrNodes), -0.005*L, [.95, .95, .95])
        %> also may use: patch(x, y, [r g b]) to control the color.
end

axis equal
axis([min(xx)-0.1 max(xx)+0.1 min(yy)-0.1 max(yy)+0.1])

end % function
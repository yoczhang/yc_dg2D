% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function PlotMesh_yctest(node, elem)

% Plot the mesh.
xx = node(:,1);
yy = node(:,2);

ElemNum = size(elem,1);


reset(cla), reset(clf), hold on

for CurrElem = 1 : ElemNum
    if iscell(elem)
        CurrNodes = elem{CurrElem};
    else
        CurrNodes = elem(CurrElem, :);
    end
    
    L = ones(1,length(CurrNodes));
    patch(xx(CurrNodes), yy(CurrNodes), -0.005*L, 'y')
end

axis equal
axis([min(xx)-0.1 max(xx)+0.1 min(yy)-0.1 max(yy)+0.1])

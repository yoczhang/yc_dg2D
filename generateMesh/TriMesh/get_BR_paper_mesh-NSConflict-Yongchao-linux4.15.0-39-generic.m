function [node, elem] = get_BR_paper_mesh(left,right,bottom,top,h_size)
%
%   First we get the quad mesh, then we get the BR's paper mesh.
%
%    ______
%   |          |
%   |______| quad mesh
%
%   ____
%   |_\/_|
%   |_/\_| BR's paper mesh
%
%
%   YcZhang 4/10/2017
%
%   Last modified 4/10/2017
%

%----------------- quad mesh ------------------------
% h_x = 1/20; h_y = h_x;
h_x = h_size; h_y = h_x;
h=[h_x,h_y];
% left=0;right=1;bottom=0;top=1;
[M,T]=generate_quad_P_T(left,right,bottom,top,h,1);
Qnode = M'; Qelem = T';
%---------------------------------------------------------

NQelems = size(Qelem,1);
NQnodes = size(Qnode,1);

Nelems = 4*NQelems;
Nnodes = NQnodes + NQelems;

node = zeros(Nnodes,2);
node(1:NQnodes,:) = Qnode;

elem = zeros(Nelems,3);

elemindex = 0;
for CurrElem = 1:NQelems
    CurrNodeindx = Qelem(CurrElem,:);
    CurrElemCoordx = Qnode(CurrNodeindx,1);
    CurrElemCoordy = Qnode(CurrNodeindx,2);
    
    baryElemx = sum(CurrElemCoordx)/4;
    baryElemy = sum(CurrElemCoordy)/4;
    
    %--- get node
    node(NQnodes+CurrElem,1) = baryElemx;
    node(NQnodes+CurrElem,2) = baryElemy;
    
    %--- get elem
    % 1-th elem
    elem(4*(CurrElem-1)+1,1) = CurrElem + NQnodes;
    elem(4*(CurrElem-1)+1,2) = CurrNodeindx(1);
    elem(4*(CurrElem-1)+1,3) = CurrNodeindx(2);
    % 2-th elem
    elem(4*(CurrElem-1)+2,1) = CurrElem + NQnodes;
    elem(4*(CurrElem-1)+2,2) = CurrNodeindx(2);
    elem(4*(CurrElem-1)+2,3) = CurrNodeindx(3);
    % 3-th elem
    elem(4*(CurrElem-1)+3,1) = CurrElem + NQnodes;
    elem(4*(CurrElem-1)+3,2) = CurrNodeindx(3);
    elem(4*(CurrElem-1)+3,3) = CurrNodeindx(4);
    % 4-th elem
    elem(4*(CurrElem-1)+4,1) = CurrElem + NQnodes;
    elem(4*(CurrElem-1)+4,2) = CurrNodeindx(4);
    elem(4*(CurrElem-1)+4,3) = CurrNodeindx(1);
    
    
%     for i = 1:4
%         elem(elemindex+i,1) = CurrElem + NQelems;
%         elem(elemindex+i,2) = CurrNodeindx(1);
%         
%         elemindex = elemindex + 1;
%     end % for i


end % for CurrElem


end % function

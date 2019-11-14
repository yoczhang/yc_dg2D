%------------------------------ PolyTop ----------------------------------%
% Ref: A Pereira, C Talischi, GH Paulino, IFM Menezes, MS Carvalho        %
%      "Fluid Flow Topology Optimization in PolyTop: Stability and        %
%      Computational Implementation", Struct Multidisc Optim,             %
%      DOI 10.1007/s00158-014-1182-z                                      %
%-------------------------------------------------------------------------%
function [BNode] = PolyBoundary(Node,Element)
NNode = size(Node,1); NElem = size(Element,1);
ElemNNode = cellfun(@length,Element); % number of nodes per element
NodeTri = zeros(NNode+NElem,2); NodeTri(1:NNode,:) = Node;
Tri = zeros(sum(ElemNNode),3); Count = 0;
for el=1:NElem
  NodeTri(NNode+el,:) = [mean(Node(Element{el},1)),...
                         mean(Node(Element{el},2))];
  nn = ElemNNode(el);
  Tri_e = zeros(nn,3);
  Tri_e(1:nn,1)=NNode + el;
  Tri_e(1:nn,2)=Element{el}(1:nn);
  Tri_e(1:(nn-1),3)=Element{el}(2:(nn-1)+1);
  Tri_e(nn,3)=Element{el}(1);
  Tri(Count+1:Count+nn,:) = Tri_e;
  Count = Count + nn;
end
TRep = TriRep(Tri,NodeTri(:,1),NodeTri(:,2));
BNode = unique(freeBoundary(TRep));
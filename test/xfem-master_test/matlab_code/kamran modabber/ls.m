function [shomare_node,ls] = ls(Qq,node,numnode)
shomare_node=zeros(numnode,2);
shomare_node(1:numnode,1)=1:numnode;
for j=1:size(Qq,1)-1
for i = 1 : numnode
    x = node(i,1);
    y = node(i,2);
    if x>=Qq(j,1) & x<=Qq(j+1,1)
       shomare_node(i,2)=j;
    end
end
end
for i = 1 : numnode
 
    if shomare_node(i,2)==0;

      shomare_node(i,2)=size(Qq,1)-1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%t   = 1/norm(seg)*seg;
for i = 1 : numnode
   x0=Qq(shomare_node(i,2),1);    y0=Qq(shomare_node(i,2),2);
   x1=Qq(shomare_node(i,2)+1,1);    y1=Qq(shomare_node(i,2)+1,2);
    
    x = node(i,1);
    y = node(i,2);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    ls(i,1) = phi/l;            % normal LS
   % ls(i,2) = ([x y]-xTip)*t';  % tangent LS
end


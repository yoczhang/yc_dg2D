function [DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT)


DGE=zeros(6,3*size(DGT,2));
DG_edge_flag=zeros(4,3*size(DGT,2));

for i=1:size(DGT,2)
    DGE(1:5,3*i-2)=[DGM(1,DGT(1,i)) DGM(2,DGT(1,i)) DGM(1,DGT(2,i)) DGM(2,DGT(2,i)) i];
    DGE(1:5,3*i-1)=[DGM(1,DGT(2,i)) DGM(2,DGT(2,i)) DGM(1,DGT(3,i)) DGM(2,DGT(3,i)) i];
    DGE(1:5,3*i)=[DGM(1,DGT(3,i)) DGM(2,DGT(3,i)) DGM(1,DGT(1,i)) DGM(2,DGT(1,i)) i];
end

edge_middle=[(DGE(1,:)+DGE(3,:))/2;(DGE(2,:)+DGE(4,:))/2];
for i=1:size(DGE,2)
    distance=sqrt((edge_middle(1,i)*ones(1,size(edge_middle,2))-edge_middle(1,:)).^2+(edge_middle(2,i)*ones(1,size(edge_middle,2))-edge_middle(2,:)).^2);
    [un_used,temp]=find(distance<0.000001);
    if length(temp)==2
        j=setdiff(temp,i);
        DGE(6,i)=j;
        ele_i=DGE(5,i);ele_j=DGE(5,j);
        DG_edge_flag(4,i)=sign(ele_i-ele_j);
        DG_edge_flag(2:3,i)=sign(ele_i-ele_j)*[DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]/norm([DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]);
    else
        DG_edge_flag(1,i)=1;
        DG_edge_flag(2:3,i)=[DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]/norm([DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]);
    end
end


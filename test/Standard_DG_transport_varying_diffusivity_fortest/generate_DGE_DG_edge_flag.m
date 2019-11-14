function [DGE,DG_edge_flag,F_i,F_in,F_out,F_HP_i]=generate_DGE_DG_edge_flag(DGM,DGT,eipsilon)


DGE=zeros(6,3*size(DGT,2));
DG_edge_flag=zeros(4,3*size(DGT,2));

for i=1:size(DGT,2)
    DGE(1:5,3*i-2)=[DGM(1,DGT(1,i)) DGM(2,DGT(1,i)) DGM(1,DGT(2,i)) DGM(2,DGT(2,i)) i];
    DGE(1:5,3*i-1)=[DGM(1,DGT(2,i)) DGM(2,DGT(2,i)) DGM(1,DGT(3,i)) DGM(2,DGT(3,i)) i];
    DGE(1:5,3*i)=[DGM(1,DGT(3,i)) DGM(2,DGT(3,i)) DGM(1,DGT(1,i)) DGM(2,DGT(1,i)) i];
end

edge_middle=[(DGE(1,:)+DGE(3,:))/2;(DGE(2,:)+DGE(4,:))/2];
F_HP_i=[];
for i=1:size(DGE,2)
    distance=sqrt((edge_middle(1,i)*ones(1,size(edge_middle,2))-edge_middle(1,:)).^2+(edge_middle(2,i)*ones(1,size(edge_middle,2))-edge_middle(2,:)).^2);
    [un_used,temp]=find(distance<0.000001);
    if length(temp)==2
        j=setdiff(temp,i);
        DGE(6,i)=j;
        
        eipsilon_i=eipsilon(1,DGE(5,i));eipsilon_j=eipsilon(1,DGE(5,j));
        if abs(eipsilon_i-eipsilon_j)>0
            DG_edge_flag(2:3,i)=sign(eipsilon_i-eipsilon_j)*[DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]'/norm([DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]);
            DG_edge_flag(4,i)=sign(eipsilon_i-eipsilon_j);
            if DG_edge_flag(2,i)<0  % means beta*n_F<0
                F_HP_i=[F_HP_i  i];
            end
        else
            ele_i=DGE(5,i);ele_j=DGE(5,j);
            DG_edge_flag(4,i)=sign(ele_i-ele_j);
            DG_edge_flag(2:3,i)=sign(ele_i-ele_j)*[DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]'/norm([DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]);
        end
    else
        DG_edge_flag(1,i)=1;
        DG_edge_flag(2:3,i)=[DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]'/norm([DGE(4,i)-DGE(2,i) DGE(1,i)-DGE(3,i)]);
    end
end

[un_used,interior_edges]=find(DGE(6,:));
F_i=setdiff(interior_edges,F_HP_i);

boundary_edges=setdiff(1:size(DGE,2),interior_edges);
[un_used,temp_1]=find(edge_middle(1,boundary_edges(:))==0);
F_in=boundary_edges(temp_1(:));

[un_used,temp_2]=find(edge_middle(1,boundary_edges(:))==2);
F_out=boundary_edges(temp_2(:));

% F_out=setdiff(boundary_edges,F_in);




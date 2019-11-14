function [new_DGM,new_DGT,new_eipsilon]=generate_refined_mesh_eipsilon(DGM,DGT,eipsilon,refinement)


%% Refinement
DGM_temp=DGM;DGT_temp=DGT;
for i=1:refinement
    [new_DGM,new_DGT]=uniformrefine_triangle(DGM_temp,DGT_temp);
    DGM_temp=new_DGM;DGT_temp=new_DGT;
end
clear DGM_temp DGT_temp


new_eipsilon=zeros(1,size(new_DGT,2));
for i=1:512
    new_eipsilon(1+(i-1)*4^refinement:i*4^refinement)=eipsilon(1,i)*ones(1,4^refinement);
end

% figure(3333);
% for i=1:size(new_DGT,2)
%     x=new_DGM(1,new_DGT(:,i));
%     y=new_DGM(2,new_DGT(:,i));
%     if new_eipsilon(1,i)==0
%        fill(x,y,'k');
%     elseif new_eipsilon(1,i)==1
%        fill(x,y,'w');
%     else
%        fill(x,y,[.5 .5 .5]);
%     end
%     hold on
% end
% axis equal, axis off
% % title('Refined Mesh');



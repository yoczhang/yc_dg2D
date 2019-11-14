function [DGM,DGT]=generate_DGM_DGT(M,T)
%DG_edge_flag(1,i)=0: the ith edge is an interior edge
%DG_edge_flag(1,i)=1: the ith edge is a Dirichelet boundary edge
%DG_edge_flag(1,i)=2: the ith edge is a Newmann boundary edge


%%
DGM=zeros(2,4*size(T,2));
DGT=zeros(4,size(T,2));

for i=1:size(T,2)
    DGT(1:4,i)=[4*i-3 4*i-2 4*i-1 4*i]';
end

for i=1:size(T,2)
    for j=1:4
        DGM(:,DGT(j,i))=M(:,T(j,i));
    end 
end

end % function



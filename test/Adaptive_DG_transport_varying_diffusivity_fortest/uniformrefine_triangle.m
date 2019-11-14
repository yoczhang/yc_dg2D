function [new_DGM,new_DGT]=uniformrefine_triangle(DGM,DGT)

new_DGT=zeros(3,4*size(DGT,2));
new_DGM=zeros(2,3*size(new_DGT,2));

for i=1:size(DGT,2)
    new_DGT(1:3,4*i-3)=[12*i-11,12*i-10,12*i-9]';
    new_DGT(1:3,4*i-2)=[12*i-8,12*i-7,12*i-6]';
    new_DGT(1:3,4*i-1)=[12*i-5,12*i-4,12*i-3]';
    new_DGT(1:3,4*i)=[12*i-2,12*i-1,12*i]';
    
    % 1st element
    new_DGM(1:2,12*i-11)=(DGM(1:2,DGT(1,i))+DGM(1:2,DGT(2,i)))/2;
    new_DGM(1:2,12*i-10)=DGM(1:2,DGT(2,i));
    new_DGM(1:2,12*i-9)=(DGM(1:2,DGT(2,i))+DGM(1:2,DGT(3,i)))/2;
    
    % 2nd element
    new_DGM(1:2,12*i-8)=(DGM(1:2,DGT(1,i))+DGM(1:2,DGT(2,i)))/2;
    new_DGM(1:2,12*i-7)=(DGM(1:2,DGT(2,i))+DGM(1:2,DGT(3,i)))/2;
    new_DGM(1:2,12*i-6)=DGM(1:2,DGT(1,i));
    
    % 3rd element
    new_DGM(1:2,12*i-5)=(DGM(1:2,DGT(1,i))+DGM(1:2,DGT(3,i)))/2;
    new_DGM(1:2,12*i-4)=DGM(1:2,DGT(1,i));
    new_DGM(1:2,12*i-3)=(DGM(1:2,DGT(2,i))+DGM(1:2,DGT(3,i)))/2;
    
    % 4th element
    new_DGM(1:2,12*i-2)=(DGM(1:2,DGT(1,i))+DGM(1:2,DGT(3,i)))/2;
    new_DGM(1:2,12*i-1)=(DGM(1:2,DGT(2,i))+DGM(1:2,DGT(3,i)))/2;
    new_DGM(1:2,12*i)=DGM(1:2,DGT(3,i));
end


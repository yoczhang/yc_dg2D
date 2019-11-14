function [new_DGM,new_DGT]=uniformrefine_quadrilateral(DGM,DGT)

new_DGT=zeros(4,4*size(DGT,2));
new_DGM=zeros(2,4*size(new_DGT,2));

for i=1:size(DGT,2)
    % Four new elements
    new_DGT(1:4,4*i-3)=[16*i-15,16*i-14,16*i-13,16*i-12]';
    new_DGT(1:4,4*i-2)=[16*i-11,16*i-10,16*i-9,16*i-8]';
    new_DGT(1:4,4*i-1)=[16*i-7,16*i-6,16*i-5,16*i-4]';
    new_DGT(1:4,4*i)=[16*i-3,16*i-2,16*i-1,16*i]';
    
    % 1st element
    new_DGM(1:2,16*i-15)=DGM(1:2,DGT(1,i));
    new_DGM(1:2,16*i-14)=(DGM(1:2,DGT(1,i))+DGM(1:2,DGT(2,i)))/2;
    new_DGM(1:2,16*i-13)=sum(DGM(1:2,DGT(:,i)),2)/4;
    new_DGM(1:2,16*i-12)=(DGM(1:2,DGT(1,i))+DGM(1:2,DGT(4,i)))/2;
    
    % 3rd element
    new_DGM(1:2,16*i-7)=(DGM(1:2,DGT(1,i))+DGM(1:2,DGT(2,i)))/2;
    new_DGM(1:2,16*i-6)=DGM(1:2,DGT(2,i));
    new_DGM(1:2,16*i-5)=(DGM(1:2,DGT(2,i))+DGM(1:2,DGT(3,i)))/2;
    new_DGM(1:2,16*i-4)=sum(DGM(1:2,DGT(:,i)),2)/4;
    
    % 2nd element
    new_DGM(1:2,16*i-11)=(DGM(1:2,DGT(1,i))+DGM(1:2,DGT(4,i)))/2;
    new_DGM(1:2,16*i-10)=sum(DGM(1:2,DGT(:,i)),2)/4;
    new_DGM(1:2,16*i-9)=(DGM(1:2,DGT(3,i))+DGM(1:2,DGT(4,i)))/2;
    new_DGM(1:2,16*i-8)=DGM(1:2,DGT(4,i));
    
    % 4th element
    new_DGM(1:2,16*i-3)=sum(DGM(1:2,DGT(:,i)),2)/4;
    new_DGM(1:2,16*i-2)=(DGM(1:2,DGT(2,i))+DGM(1:2,DGT(3,i)))/2;
    new_DGM(1:2,16*i-1)=DGM(1:2,DGT(3,i));
    new_DGM(1:2,16*i)=(DGM(1:2,DGT(3,i))+DGM(1:2,DGT(4,i)))/2;
end


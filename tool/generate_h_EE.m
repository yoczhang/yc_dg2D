function r=generate_h_EE(M,T)
center=zeros(2,size(T,2) );
for n=1:size(T,2) 
    T_boundary=cell2mat(T(:,n))'; 
    numbe_of_element_boundary=length(T_boundary);
    vertices=M(:,T_boundary(:,1));
    center(1,n)=sum( vertices(1,:))/numbe_of_element_boundary;
    center(2,n)=sum( vertices(2,:))/numbe_of_element_boundary;
end

r=[];
for i=1:size(center,2)
    s=repmat(center(:,i),size(cell2mat(T(:,i))))-M(:,cell2mat(T(:,i)));
    r=[r max(sqrt(sum(s.*s)))];
end

end % function

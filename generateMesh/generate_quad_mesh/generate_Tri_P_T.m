function [Pb,Tb]=generate_Tri_P_T(left,right,bottom,top,h)

% left = 0;
% right = 2;
% bottom = 0;
% top = 1;
% h = [1/64,1/64];
basis_type=1;

if basis_type==1

   N1=(right-left)/h(1);
   N2=(top-bottom)/h(2);
   tnp=(N1+1)*(N2+1);
   Pb=zeros(2,tnp);
   Tb=zeros(3,2*N1*N2);
   Q=zeros(N1+1,N2+1);

   for j=1:tnp
      if mod(j,N2+1)==0
         Pb(1,j)=left+(j/(N2+1)-1)*h(1);
         Pb(2,j)=top;
      else
         Pb(1,j)=left+fix(j/(N2+1))*h(1);
         Pb(2,j)=bottom+(mod(j,N2+1)-1)*h(2);
      end
   end

   for i=1:N1+1
      for j=1:N2+1
         Q(i,j)=(i-1)*(N2+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
   for n=1:N1*N2
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      Tb(1,2*n-1)=Q(column,row);
      Tb(2,2*n-1)=Q(column+1,row);
      Tb(3,2*n-1)=Q(column,row+1);  
  
      Tb(1,2*n)=Q(column,row+1);
      Tb(2,2*n)=Q(column+1,row);
      Tb(3,2*n)=Q(column+1,row+1);  
    
   end

elseif basis_type==2

   N1=(right-left)/h(1);
   N2=(top-bottom)/h(2);
   dh=h/2;
   dN1=N1*2;
   dN2=N2*2;
   tnp=(dN1+1)*(dN2+1);
   Pb=zeros(2,tnp);
   Tb=zeros(3,2*N1*N2);
   Q=zeros(dN1+1,dN2+1);

   for j=1:tnp
      if mod(j,dN2+1)==0
         Pb(1,j)=left+(j/(dN2+1)-1)*dh(1);
         Pb(2,j)=top;
      else
         Pb(1,j)=left+fix(j/(dN2+1))*dh(1);
         Pb(2,j)=bottom+(mod(j,dN2+1)-1)*dh(2);
      end
   end

   for i=1:dN1+1
      for j=1:dN2+1
         Q(i,j)=(i-1)*(dN2+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
   for n=1:N1*N2
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      Tb(1,2*n-1)=Q(2*column-1,2*row-1);
      Tb(2,2*n-1)=Q(2*column+1,2*row-1); 
      Tb(3,2*n-1)=Q(2*column-1,2*row+1);
      Tb(4,2*n-1)=Q(2*column,2*row-1);
      Tb(5,2*n-1)=Q(2*column,2*row);
      Tb(6,2*n-1)=Q(2*column-1,2*row);


      Tb(1,2*n)=Q(2*column-1,2*row+1);
      Tb(2,2*n)=Q(2*column+1,2*row-1);
      Tb(3,2*n)=Q(2*column+1,2*row+1);
      Tb(4,2*n)=Q(2*column,2*row);
      Tb(5,2*n)=Q(2*column+1,2*row);
      Tb(6,2*n)=Q(2*column,2*row+1); 

   end
end

Pb = Pb'; Tb = Tb';
end % function



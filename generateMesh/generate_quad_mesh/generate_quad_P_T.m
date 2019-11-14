function [M,T]=generate_quad_P_T(a,c,d,e,h_partition,basis_type)
%
%   %-------------------------------------------------
%       Just copy from Xiaoming He's code.
%   %-------------------------------------------------
%
%Generate the matrix M and T
%M stores the coordinates of all the grid points for the grid with step size h_partition
%M(i,j) is the ith coordinate of the jth grid point
%T stores the ordinal of the grid points of every element for the grid with step size h_partition
%T(i,j) is the ordinal of the ith grid point in th jth element
%The rules for the total ordinal for all the grid points and the local ordinal for the grid points in a element 
%are in my 'Standard 2-dimensional FE Tool boxes notes part 2'.
%note: Sometimes the grid points of the partition is not the same as the grid points of the basis.
%For example,Lagrange biquadratic FE has 9 grid points in every element but the partition still has only 4 grid points for every element,i.e, the 4 vertices.
%So,we have 2 kinds of M and T, which are M_partition,T_partition and M_basis, T_basis
%M_partition: store the coordinates of all the grid points of the partition,not the basis
%T_partition: store the ordinal of the grid points of every element for the partition,not the basis
%M_basis: store the coordinates of all the grid points of the basis,not the partition
%T_basis: store the ordinal of the grid points of every element for the basis,not the partition
%For bilinear FE,M_partition,T_partition are the same as M_basis,T_basis.
%For Langrange FE, they are different.The following 4 rows are what they do for Langrange FE
%M_partition only stores the coordinates of the 4 vertices for every element.
%M_basis stores the coordinates of all the 9 grid points for every element.
%T_partiton only stores the ordinal of the 4 vertices for every element
%T_basis stores the ordinal of the all the 9 grid points for every element
%note: the total ordinal for T_partiton is different the total ordinal for T_basis
%the total ordinal for T_partiton only counts the vertices of all the element
%the total ordinal for T_basis counts all the grid points for Lagrange biquadratic FE
%a: the left boundary of the retangular domain
%c: the right boundary of the retangular domain
%d: the lower boundary of the retangular domain
%e: the upper boundary of the retangular domain
%h_partition: the step size of the partition
%basis_type: the type of the basis
%basis_type=1:bilinear FE
%basis_type=2:Lagrange biquadratic FE

%N1 is the number of the intervals of the division in x-direction
%N2 is the number of the intervals of the division in y-direction
%tnp:total number of all the grid points,including the inner and boundary points
%The explanation of Q,row,column is in my Standard 2-dimension FE tool boxes notes part 2(6)

h=h_partition;
if basis_type==1
   N1=(c-a)/h(1);
   N2=(e-d)/h(2);
   tnp=(N1+1)*(N2+1);

   for j=1:tnp
      if mod(j,N2+1)==0
         M(1,j)=a+(fix(j/(N2+1))-1)*h(1);
         M(2,j)=e;
      else
         M(1,j)=a+fix(j/(N2+1))*h(1);
         M(2,j)=d+(mod(j,N2+1)-1)*h(2);
      end
   end

   for i=1:N1+1
      for j=1:N2+1
         Q(i,j)=(i-1)*(N2+1)+j;
      end
   end

   for n=1:N1*N2
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      T(1,n)=Q(column,row);
      T(2,n)=Q(column+1,row);
      T(3,n)=Q(column+1,row+1);    
      T(4,n)=Q(column,row+1);      
   end

elseif basis_type==2
   N1=(c-a)/h(1);
   N2=(e-d)/h(2);
   dh=h/2;
   dN1=N1*2;
   dN2=N2*2;
   tnp=(dN1+1)*(dN2+1);

   for j=1:tnp
      if mod(j,dN2+1)==0
         M(1,j)=a+fix(j/(dN2+1)-1)*dh(1);
         M(2,j)=e;
      else
         M(1,j)=a+fix(j/(dN2+1))*dh(1);
         M(2,j)=d+(mod(j,dN2+1)-1)*dh(2);
      end
   end

   for i=1:dN1+1
      for j=1:dN2+1
         Q(i,j)=(i-1)*(dN2+1)+j;
      end
   end

   for n=1:N1*N2
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      T(1,n)=Q(2*column-1,2*row-1);
      T(2,n)=Q(2*column+1,2*row-1); 
      T(3,n)=Q(2*column+1,2*row+1);
      T(4,n)=Q(2*column-1,2*row+1);
      T(5,n)=Q(2*column,2*row-1);
      T(6,n)=Q(2*column+1,2*row);
      T(7,n)=Q(2*column,2*row+1);
      T(8,n)=Q(2*column-1,2*row);
      T(9,n)=Q(2*column,2*row);
   end
end

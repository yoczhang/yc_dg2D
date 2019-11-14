function [Gauss_coefficient_local,Gauss_point_local]=generate_Gauss_local(Gauss_coefficient_reference,Gauss_point_reference,left_lower_point,h_partition,vertices)
%Generate the Gauss coefficients and Gauss points on the local retangular element by using affine tranformation
%Gauss_coefficient_local,Gauss_point_local:the Gauss coefficients and Gauss points on the local retangular element
%Some details are in my standard 2-dimension FE Tool boxes notes part 6
%Gauss_coefficient_reference,Gauss_point_reference: the Gauss coefficients and Gauss points on the reference square 
%left_lower_point: The coordinates of the left-lower vertice of the current element whose local basis we are evaluating
%h_partition: the step size of the partition

% Jacobi=h_partition(1)*h_partition(2)/4;
% Gauss_coefficient_local=Gauss_coefficient_reference*Jacobi;
% Gauss_point_local(:,1)=h_partition(1)*Gauss_point_reference(:,1)/2+left_lower_point(1)+h_partition(1)/2;
% Gauss_point_local(:,2)=h_partition(2)*Gauss_point_reference(:,2)/2+left_lower_point(2)+h_partition(2)/2;
x_1=vertices(1,1);x_2=vertices(1,2);x_3=vertices(1,3);x_4=vertices(1,4);
y_1=vertices(2,1);y_2=vertices(2,2);y_3=vertices(2,3);y_4=vertices(2,4);
for i=1:length(Gauss_coefficient_reference)
    Gauss_point_local(i,1)=vertices(1,1)*(1-Gauss_point_reference(i,1))*(1-Gauss_point_reference(i,2))/4....
                          +vertices(1,2)*(1+Gauss_point_reference(i,1))*(1-Gauss_point_reference(i,2))/4....
                          +vertices(1,3)*(1+Gauss_point_reference(i,1))*(1+Gauss_point_reference(i,2))/4....
                          +vertices(1,4)*(1-Gauss_point_reference(i,1))*(1+Gauss_point_reference(i,2))/4;
    Gauss_point_local(i,2)=vertices(2,1)*(1-Gauss_point_reference(i,1))*(1-Gauss_point_reference(i,2))/4....
                          +vertices(2,2)*(1+Gauss_point_reference(i,1))*(1-Gauss_point_reference(i,2))/4....
                          +vertices(2,3)*(1+Gauss_point_reference(i,1))*(1+Gauss_point_reference(i,2))/4....
                          +vertices(2,4)*(1-Gauss_point_reference(i,1))*(1+Gauss_point_reference(i,2))/4;
    x_bar=Gauss_point_reference(i,1);
    y_bar=Gauss_point_reference(i,2);
    Gauss_coefficient_local(i)=Gauss_coefficient_reference(i)*abs(-(x_2*y_1 - x_1*y_2 + x_1*y_4 - x_2*y_3 + x_3*y_2 - x_4*y_1...
                              - x_3*y_4 + x_4*y_3 + x_1*x_bar*y_3 - x_3*x_bar*y_1 - x_1*x_bar*y_4 - x_2*x_bar*y_3 + x_3*x_bar*y_2 ...
                              + x_4*x_bar*y_1 + x_2*x_bar*y_4 - x_4*x_bar*y_2 + x_1*y_2*y_bar - x_2*y_1*y_bar - x_1*y_3*y_bar ...
                              + x_3*y_1*y_bar + x_2*y_4*y_bar - x_4*y_2*y_bar - x_3*y_4*y_bar + x_4*y_3*y_bar)/8);
end

% temp_1=[vertices(1,2)-vertices(1,1) vertices(1,3)-vertices(1,1);vertices(2,2)-vertices(2,1) vertices(2,3)-vertices(2,1)];
% temp_2=[vertices(1,4)-vertices(1,1) vertices(1,3)-vertices(1,1);vertices(2,4)-vertices(2,1) vertices(2,3)-vertices(2,1)];
% area=abs(det(temp_1))/2+abs(det(temp_2))/2;
% Gauss_coefficient_local=Gauss_coefficient_reference;


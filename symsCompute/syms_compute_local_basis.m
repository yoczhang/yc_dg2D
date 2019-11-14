function syms_compute_local_basis
clc
clear

%% 2D
% syms x_b xT y_b yT hT
% x = (x_b-xT)/hT;
% y = (y_b-yT)/hT;
% 
% Pb1 = [1, x, y]
% Pb1_x = diff(Pb1,x_b)
% Pb1_y = diff(Pb1,y_b)
% 
% 
% Pb2 = [1 x y x.^2 x.*y y.^2]
% Pb2_x = diff(Pb2,x_b)
% Pb2_y = diff(Pb2,y_b)
% 
% Pb3 = [1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3]
% Pb3_x = diff(Pb3,x_b)
% Pb3_y = diff(Pb3,y_b)
% 
% Pb4 = [1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3	y.^4]
% Pb4_x = diff(Pb4,x_b)
% Pb4_y = diff(Pb4,y_b)
% 
% Pb5 = [...
%     1, ...
%     x, y, ...
%     x^2, x*y,  y^2, ...
%     x^3, x^2*y, x*y^2, y^3, ...
%     x^4, x^3*y, x^2*y^2, x*y^3, y^4, ...
%     x^5, x^4*y, x^3*y^2, x^2*y^3, x*y^4, y^5]
% Pb5_x = diff(Pb5,x_b)
% Pb5_y = diff(Pb5,y_b)
% 
% Pb6 = [...
%     1, ...
%     x, y, ...
%     x^2, x*y,  y^2, ...
%     x^3, x^2*y, x*y^2, y^3, ...
%     x^4, x^3*y, x^2*y^2, x*y^3, y^4, ...
%     x^5, x^4*y, x^3*y^2, x^2*y^3, x*y^4, y^5, ...
%     x^6, x^5*y, x^4*y^2, x^3*y^3, x^2*y^4, x*y^5, y^6]
% Pb6_x = diff(Pb6,x_b)
% Pb6_y = diff(Pb6,y_b)
% 
% Pb7 = [...
%     1, ...
%     x, y, ...
%     x^2, x*y,  y^2, ...
%     x^3, x^2*y, x*y^2, y^3, ...
%     x^4, x^3*y, x^2*y^2, x*y^3, y^4, ...
%     x^5, x^4*y, x^3*y^2, x^2*y^3, x*y^4, y^5, ...
%     x^6, x^5*y, x^4*y^2, x^3*y^3, x^2*y^4, x*y^5, y^6, ...
%     x^7, x^6*y, x^5*y^2, x^4*y^3, x^3*y^4, x^2*y^5, x*y^6, y^7]
% Pb7_x = diff(Pb7,x_b)
% Pb7_y = diff(Pb7,y_b)



%% 1D
syms x xE hE
xb = (x-xE)/hE;

Pb1 = [1, xb]
Pb1_x = diff(Pb1,x)

Pb2 = [1, xb, xb^2]
Pb2_x = diff(Pb2,x)

Pb3 = [1, xb, xb^2, xb^3]
Pb3_x = diff(Pb3,x)

Pb4 = [1, xb, xb^2, xb^3, xb^4]
Pb4_x = diff(Pb4,x)

Pb5 = [1, xb, xb^2, xb^3, xb^4, xb^5]
Pb5_x = diff(Pb5,x)

Pb6 = [1, xb, xb^2, xb^3, xb^4, xb^5, xb^6]
Pb6_x = diff(Pb6,x)

Pb7 = [1, xb, xb^2, xb^3, xb^4, xb^5, xb^6, xb^7]
Pb7_x = diff(Pb7,x)


end % function
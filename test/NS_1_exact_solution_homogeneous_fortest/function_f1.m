function r=function_f1(x,y)

nu=function_nu(0,0);
r=2*(2*y - 1)*(- 3*nu*x^4 + 6*nu*x^3 - 6*nu*x^2*y^2 + 6*nu*x^2*y - 3*nu*x^2 + 6*nu*x*y^2 - 6*nu*x*y - nu*y^2 + nu*y + 1)...
    +function_u1(x,y)*function_u1_x(x,y)+function_u2(x,y)*function_u1_y(x,y);
function r=function_f2(x,y)

nu=function_nu(0,0);
r=2*(2*x - 1)*(6*nu*x^2*y^2 - 6*nu*x^2*y + nu*x^2 - 6*nu*x*y^2 + 6*nu*x*y - nu*x + 3*nu*y^4 - 6*nu*y^3 + 3*nu*y^2 + 1)...
    +function_u1(x,y)*function_u2_x(x,y)+function_u2(x,y)*function_u2_y(x,y);
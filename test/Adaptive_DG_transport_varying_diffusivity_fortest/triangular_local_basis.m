function r=triangular_local_basis(x,y,basis_type,derivative_degree_x,derivative_degree_y)

if basis_type==1
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1 x y];
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0 1 0];
    elseif derivative_degree_x==0&&derivative_degree_y==1
        r=[0 0 1];   
    end
elseif basis_type==2
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1 x y x.^2 x.*y y.^2];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0 1 0 2.*x y 0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1
        r=[0 0 1 0 x 2.*y];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0 0 0 0 1 0];
    end
elseif basis_type==3
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0];
    end
elseif basis_type==4
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3 y.^4];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0      4*x.^3   3*x^2*y    2*x*y^2      y^3     0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2 0        x^3        2*x^2*y      3*x*y^2 4*y^3];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0      0        3*x^2      4*x*y        3*y^2   0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6*x     2*y      0        0      12*x^2   6*x*y      2*y^2        0       0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2*x      6*y    0        0          2*x^2        6*x*y   12*y^2];
    end
end
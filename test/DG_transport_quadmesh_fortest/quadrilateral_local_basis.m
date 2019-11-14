function r=quadrilateral_local_basis(x,y,basis_type,derivative_degree_x,derivative_degree_y)

if basis_type==1
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1 x y];
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0 1 0];
    elseif derivative_degree_x==0&&derivative_degree_y==1
        r=[0 0 1];   
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0 0 0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0 0 0]; 
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
elseif basis_type==5
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3 y.^4     x^5      x^4*y      x^3*y^2      x^2*y^3      x*y^4       y^5];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0      4*x.^3   3*x^2*y    2*x*y^2      y^3     0        5*x^4    4*x^3*y    3*x^2*y^2    2*x*y^3      y^4         0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1 
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2 0        x^3        2*x^2*y      3*x*y^2 4*y^3    0        x^4        2*x^3*y      3*x^2*y^2    4*x*y^3     5*y^4];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0      0        3*x^2      4*x*y        3*y^2   0        0        4*x^3      6*x^2*y      6*x*y^2      4*y^3       0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6*x     2*y      0        0      12*x^2   6*x*y      2*y^2        0       0        20*x^3   12*x^2*y   6*x*y^2      2*y^3        0           0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2*x      6*y    0        0          2*x^2        6*x*y   12*y^2   0        0          2*x^3        6*x^2*y     12*x*y^2     20*y^3];
    end
elseif basis_type==6
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3 y.^4     x^5      x^4*y      x^3*y^2      x^2*y^3      x*y^4       y^5      x^6        x^5*y        x^4*y^2       x^3*y^3       x^2*y^4        x*y^5        y^6];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0      4*x.^3   3*x^2*y    2*x*y^2      y^3     0        5*x^4    4*x^3*y    3*x^2*y^2    2*x*y^3      y^4         0        6*x^5      5*x^4*y      4*x^3*y^2     3*x^2*y^3     2*x*y^4        y^5          0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1 
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2 0        x^3        2*x^2*y      3*x*y^2 4*y^3    0        x^4        2*x^3*y      3*x^2*y^2    4*x*y^3     5*y^4    0          x^5          2*x^4*y       3*x^3*y^2     4*x^2*y^3      5*x*y^4      6*y^5];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0      0        3*x^2      4*x*y        3*y^2   0        0        4*x^3      6*x^2*y      6*x*y^2      4*y^3       0        0          5*x^4        8*x^3*y       9*x^2*y^2     8*x*y^3        5*y^4        0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6*x     2*y      0        0      12*x^2   6*x*y      2*y^2        0       0        20*x^3   12*x^2*y   6*x*y^2      2*y^3        0           0        30*x^4     20*x^3*y     12*x^2*y^2    6*x*y^3       2*y^4          0            0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2*x      6*y    0        0          2*x^2        6*x*y   12*y^2   0        0          2*x^3        6*x^2*y     12*x*y^2     20*y^3   0          0            2*x^4         6*x^3*y       12*x^2*y^2     20*x*y^3     30*y^4];
    end
elseif basis_type==7
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3 y.^4     x^5      x^4*y      x^3*y^2      x^2*y^3      x*y^4       y^5      x^6        x^5*y        x^4*y^2       x^3*y^3       x^2*y^4        x*y^5        y^6      x^7     x^6*y     x^5*y^2     x^4*y^3       x^3*y^4     x^2*y^5    x*y^6     y^7];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0      4*x.^3   3*x^2*y    2*x*y^2      y^3     0        5*x^4    4*x^3*y    3*x^2*y^2    2*x*y^3      y^4         0        6*x^5      5*x^4*y      4*x^3*y^2     3*x^2*y^3     2*x*y^4        y^5          0        7*x^6   6*x^5*y   5*x^4*y^2   4*x^3*y^3     3*x^2*y^4   2*x*y^5    y^6       0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1 
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2 0        x^3        2*x^2*y      3*x*y^2 4*y^3    0        x^4        2*x^3*y      3*x^2*y^2    4*x*y^3     5*y^4    0          x^5          2*x^4*y       3*x^3*y^2     4*x^2*y^3      5*x*y^4      6*y^5    0       x^6       2*x^5*y     3*x^4*y^2     4*x^3*y^3   5*x^2*y^4  6*x*y^5   7*y^6];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0      0        3*x^2      4*x*y        3*y^2   0        0        4*x^3      6*x^2*y      6*x*y^2      4*y^3       0        0          5*x^4        8*x^3*y       9*x^2*y^2     8*x*y^3        5*y^4        0        0       6*x^5     10*x^4*y    12*x^3*y^2    12*x^2*y^3  10*x*y^4   6*y^5     0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6*x     2*y      0        0      12*x^2   6*x*y      2*y^2        0       0        20*x^3   12*x^2*y   6*x*y^2      2*y^3        0           0        30*x^4     20*x^3*y     12*x^2*y^2    6*x*y^3       2*y^4          0            0        42*x^5  30*x^4*y  20*x^3*y^2  12*x^2*y^3    6*x*y^4     2*y^5      0         0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2*x      6*y    0        0          2*x^2        6*x*y   12*y^2   0        0          2*x^3        6*x^2*y     12*x*y^2     20*y^3   0          0            2*x^4         6*x^3*y       12*x^2*y^2     20*x*y^3     30*y^4   0       0         2*x^5       6*x^4*y       12*x^3*y^2  20*x^2*y^3 30*x*y^4  42*y^5];
    end
elseif basis_type==8
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[1  x  y  x.^2  x.*y  y.^2   x.^3    x.^2.*y  x.*y.^2  y.^3   x.^4     x.^3.*y    x.^2.*y.^2   x.*y.^3 y.^4     x^5      x^4*y      x^3*y^2      x^2*y^3      x*y^4       y^5      x^6        x^5*y        x^4*y^2       x^3*y^3       x^2*y^4        x*y^5        y^6      x^7     x^6*y     x^5*y^2     x^4*y^3       x^3*y^4     x^2*y^5    x*y^6     y^7      x^8     x^7*y      x^6*y^2      x^5*y^3       x^4*y^4      x^3*y^5      x^2*y^6      x*y^7     y^8];         
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[0  1  0  2.*x  y     0      3*x.^2  2*x.*y   y.^2     0      4*x.^3   3*x^2*y    2*x*y^2      y^3     0        5*x^4    4*x^3*y    3*x^2*y^2    2*x*y^3      y^4         0        6*x^5      5*x^4*y      4*x^3*y^2     3*x^2*y^3     2*x*y^4        y^5          0        7*x^6   6*x^5*y   5*x^4*y^2   4*x^3*y^3     3*x^2*y^4   2*x*y^5    y^6       0        8*x^7   7*x^6*y    6*x^5*y^2    5*x^4*y^3     4*x^3*y^4    3*x^2*y^5    2*x*y^6      y^7       0];          
    elseif derivative_degree_x==0&&derivative_degree_y==1 
        r=[0  0  1  0     x     2.*y   0       x.^2     2*x.*y   3*y.^2 0        x^3        2*x^2*y      3*x*y^2 4*y^3    0        x^4        2*x^3*y      3*x^2*y^2    4*x*y^3     5*y^4    0          x^5          2*x^4*y       3*x^3*y^2     4*x^2*y^3      5*x*y^4      6*y^5    0       x^6       2*x^5*y     3*x^4*y^2     4*x^3*y^3   5*x^2*y^4  6*x*y^5   7*y^6    0       x^7        2*x^6*y      3*x^5*y^2     4*x^4*y^3    5*x^3*y^4    6*x^2*y^5    7*x*y^6   8*y^7];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[0  0  0  0     1     0      0       2*x      2*y      0      0        3*x^2      4*x*y        3*y^2   0        0        4*x^3      6*x^2*y      6*x*y^2      4*y^3       0        0          5*x^4        8*x^3*y       9*x^2*y^2     8*x*y^3        5*y^4        0        0       6*x^5     10*x^4*y    12*x^3*y^2    12*x^2*y^3  10*x*y^4   6*y^5     0        0       7*x^6      12*x^5*y     15*x^4*y^2    16*x^3*y^3   15*x^2*y^4   12*x*y^5     7*y^6     0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[0  0  0  2     0     0      6*x     2*y      0        0      12*x^2   6*x*y      2*y^2        0       0        20*x^3   12*x^2*y   6*x*y^2      2*y^3        0           0        30*x^4     20*x^3*y     12*x^2*y^2    6*x*y^3       2*y^4          0            0        42*x^5  30*x^4*y  20*x^3*y^2  12*x^2*y^3    6*x*y^4     2*y^5      0         0        56*x^6  42*x^5*y   30*x^4*y^2   20*x^3*y^3    12*x^2*y^4   6*x*y^5      2*y^6        0         0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[0  0  0  0     0     2      0       0        2*x      6*y    0        0          2*x^2        6*x*y   12*y^2   0        0          2*x^3        6*x^2*y     12*x*y^2     20*y^3   0          0            2*x^4         6*x^3*y       12*x^2*y^2     20*x*y^3     30*y^4   0       0         2*x^5       6*x^4*y       12*x^3*y^2  20*x^2*y^3 30*x*y^4  42*y^5   0       0          2*x^6        6*x^5*y       12*x^4*y^2   20*x^3*y^3   30*x^2*y^4   42*x*y^5  56*y^6];
    end
elseif basis_type==9
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[ 1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3, x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^5, x^4*y, x^3*y^2, x^2*y^3, x*y^4, y^5, x^6, x^5*y, x^4*y^2, x^3*y^3, x^2*y^4, x*y^5, y^6, x^7, x^6*y, x^5*y^2, x^4*y^3, x^3*y^4, x^2*y^5, x*y^6, y^7, x^8, x^7*y, x^6*y^2, x^5*y^3, x^4*y^4, x^3*y^5, x^2*y^6, x*y^7, y^8, x^9, x^8*y, x^7*y^2, x^6*y^3, x^5*y^4, x^4*y^5, x^3*y^6, x^2*y^7, x*y^8, y^9];
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[ 0, 1, 0, 2*x, y, 0, 3*x^2, 2*x*y, y^2, 0, 4*x^3, 3*x^2*y, 2*x*y^2, y^3, 0, 5*x^4, 4*x^3*y, 3*x^2*y^2, 2*x*y^3, y^4, 0, 6*x^5, 5*x^4*y, 4*x^3*y^2, 3*x^2*y^3, 2*x*y^4, y^5, 0, 7*x^6, 6*x^5*y, 5*x^4*y^2, 4*x^3*y^3, 3*x^2*y^4, 2*x*y^5, y^6, 0, 8*x^7, 7*x^6*y, 6*x^5*y^2, 5*x^4*y^3, 4*x^3*y^4, 3*x^2*y^5, 2*x*y^6, y^7, 0, 9*x^8, 8*x^7*y, 7*x^6*y^2, 6*x^5*y^3, 5*x^4*y^4, 4*x^3*y^5, 3*x^2*y^6, 2*x*y^7, y^8, 0];
    elseif derivative_degree_x==0&&derivative_degree_y==1 
        r=[ 0, 0, 1, 0, x, 2*y, 0, x^2, 2*x*y, 3*y^2, 0, x^3, 2*x^2*y, 3*x*y^2, 4*y^3, 0, x^4, 2*x^3*y, 3*x^2*y^2, 4*x*y^3, 5*y^4, 0, x^5, 2*x^4*y, 3*x^3*y^2, 4*x^2*y^3, 5*x*y^4, 6*y^5, 0, x^6, 2*x^5*y, 3*x^4*y^2, 4*x^3*y^3, 5*x^2*y^4, 6*x*y^5, 7*y^6, 0, x^7, 2*x^6*y, 3*x^5*y^2, 4*x^4*y^3, 5*x^3*y^4, 6*x^2*y^5, 7*x*y^6, 8*y^7, 0, x^8, 2*x^7*y, 3*x^6*y^2, 4*x^5*y^3, 5*x^4*y^4, 6*x^3*y^5, 7*x^2*y^6, 8*x*y^7, 9*y^8];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[ 0, 0, 0, 0, 1, 0, 0, 2*x, 2*y, 0, 0, 3*x^2, 4*x*y, 3*y^2, 0, 0, 4*x^3, 6*x^2*y, 6*x*y^2, 4*y^3, 0, 0, 5*x^4, 8*x^3*y, 9*x^2*y^2, 8*x*y^3, 5*y^4, 0, 0, 6*x^5, 10*x^4*y, 12*x^3*y^2, 12*x^2*y^3, 10*x*y^4, 6*y^5, 0, 0, 7*x^6, 12*x^5*y, 15*x^4*y^2, 16*x^3*y^3, 15*x^2*y^4, 12*x*y^5, 7*y^6, 0, 0, 8*x^7, 14*x^6*y, 18*x^5*y^2, 20*x^4*y^3, 20*x^3*y^4, 18*x^2*y^5, 14*x*y^6, 8*y^7, 0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[ 0, 0, 0, 2, 0, 0, 6*x, 2*y, 0, 0, 12*x^2, 6*x*y, 2*y^2, 0, 0, 20*x^3, 12*x^2*y, 6*x*y^2, 2*y^3, 0, 0, 30*x^4, 20*x^3*y, 12*x^2*y^2, 6*x*y^3, 2*y^4, 0, 0, 42*x^5, 30*x^4*y, 20*x^3*y^2, 12*x^2*y^3, 6*x*y^4, 2*y^5, 0, 0, 56*x^6, 42*x^5*y, 30*x^4*y^2, 20*x^3*y^3, 12*x^2*y^4, 6*x*y^5, 2*y^6, 0, 0, 72*x^7, 56*x^6*y, 42*x^5*y^2, 30*x^4*y^3, 20*x^3*y^4, 12*x^2*y^5, 6*x*y^6, 2*y^7, 0, 0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[ 0, 0, 0, 0, 0, 2, 0, 0, 2*x, 6*y, 0, 0, 2*x^2, 6*x*y, 12*y^2, 0, 0, 2*x^3, 6*x^2*y, 12*x*y^2, 20*y^3, 0, 0, 2*x^4, 6*x^3*y, 12*x^2*y^2, 20*x*y^3, 30*y^4, 0, 0, 2*x^5, 6*x^4*y, 12*x^3*y^2, 20*x^2*y^3, 30*x*y^4, 42*y^5, 0, 0, 2*x^6, 6*x^5*y, 12*x^4*y^2, 20*x^3*y^3, 30*x^2*y^4, 42*x*y^5, 56*y^6, 0, 0, 2*x^7, 6*x^6*y, 12*x^5*y^2, 20*x^4*y^3, 30*x^3*y^4, 42*x^2*y^5, 56*x*y^6, 72*y^7];
    end
elseif basis_type==10
    if derivative_degree_x==0&&derivative_degree_y==0
        r=[ 1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3, x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^5, x^4*y, x^3*y^2, x^2*y^3, x*y^4, y^5, x^6, x^5*y, x^4*y^2, x^3*y^3, x^2*y^4, x*y^5, y^6, x^7, x^6*y, x^5*y^2, x^4*y^3, x^3*y^4, x^2*y^5, x*y^6, y^7, x^8, x^7*y, x^6*y^2, x^5*y^3, x^4*y^4, x^3*y^5, x^2*y^6, x*y^7, y^8, x^9, x^8*y, x^7*y^2, x^6*y^3, x^5*y^4, x^4*y^5, x^3*y^6, x^2*y^7, x*y^8, y^9, x^10, x^9*y, x^8*y^2, x^7*y^3, x^6*y^4, x^5*y^5, x^4*y^6, x^3*y^7, x^2*y^8, x*y^9, y^10];
    elseif derivative_degree_x==1&&derivative_degree_y==0
        r=[ 0, 1, 0, 2*x, y, 0, 3*x^2, 2*x*y, y^2, 0, 4*x^3, 3*x^2*y, 2*x*y^2, y^3, 0, 5*x^4, 4*x^3*y, 3*x^2*y^2, 2*x*y^3, y^4, 0, 6*x^5, 5*x^4*y, 4*x^3*y^2, 3*x^2*y^3, 2*x*y^4, y^5, 0, 7*x^6, 6*x^5*y, 5*x^4*y^2, 4*x^3*y^3, 3*x^2*y^4, 2*x*y^5, y^6, 0, 8*x^7, 7*x^6*y, 6*x^5*y^2, 5*x^4*y^3, 4*x^3*y^4, 3*x^2*y^5, 2*x*y^6, y^7, 0, 9*x^8, 8*x^7*y, 7*x^6*y^2, 6*x^5*y^3, 5*x^4*y^4, 4*x^3*y^5, 3*x^2*y^6, 2*x*y^7, y^8, 0, 10*x^9, 9*x^8*y, 8*x^7*y^2, 7*x^6*y^3, 6*x^5*y^4, 5*x^4*y^5, 4*x^3*y^6, 3*x^2*y^7, 2*x*y^8, y^9, 0];
    elseif derivative_degree_x==0&&derivative_degree_y==1 
        r=[ 0, 0, 1, 0, x, 2*y, 0, x^2, 2*x*y, 3*y^2, 0, x^3, 2*x^2*y, 3*x*y^2, 4*y^3, 0, x^4, 2*x^3*y, 3*x^2*y^2, 4*x*y^3, 5*y^4, 0, x^5, 2*x^4*y, 3*x^3*y^2, 4*x^2*y^3, 5*x*y^4, 6*y^5, 0, x^6, 2*x^5*y, 3*x^4*y^2, 4*x^3*y^3, 5*x^2*y^4, 6*x*y^5, 7*y^6, 0, x^7, 2*x^6*y, 3*x^5*y^2, 4*x^4*y^3, 5*x^3*y^4, 6*x^2*y^5, 7*x*y^6, 8*y^7, 0, x^8, 2*x^7*y, 3*x^6*y^2, 4*x^5*y^3, 5*x^4*y^4, 6*x^3*y^5, 7*x^2*y^6, 8*x*y^7, 9*y^8, 0, x^9, 2*x^8*y, 3*x^7*y^2, 4*x^6*y^3, 5*x^5*y^4, 6*x^4*y^5, 7*x^3*y^6, 8*x^2*y^7, 9*x*y^8, 10*y^9];
    elseif derivative_degree_x==1&&derivative_degree_y==1 
        r=[ 0, 0, 0, 0, 1, 0, 0, 2*x, 2*y, 0, 0, 3*x^2, 4*x*y, 3*y^2, 0, 0, 4*x^3, 6*x^2*y, 6*x*y^2, 4*y^3, 0, 0, 5*x^4, 8*x^3*y, 9*x^2*y^2, 8*x*y^3, 5*y^4, 0, 0, 6*x^5, 10*x^4*y, 12*x^3*y^2, 12*x^2*y^3, 10*x*y^4, 6*y^5, 0, 0, 7*x^6, 12*x^5*y, 15*x^4*y^2, 16*x^3*y^3, 15*x^2*y^4, 12*x*y^5, 7*y^6, 0, 0, 8*x^7, 14*x^6*y, 18*x^5*y^2, 20*x^4*y^3, 20*x^3*y^4, 18*x^2*y^5, 14*x*y^6, 8*y^7, 0, 0, 9*x^8, 16*x^7*y, 21*x^6*y^2, 24*x^5*y^3, 25*x^4*y^4, 24*x^3*y^5, 21*x^2*y^6, 16*x*y^7, 9*y^8, 0];
    elseif derivative_degree_x==2&&derivative_degree_y==0
        r=[ 0, 0, 0, 2, 0, 0, 6*x, 2*y, 0, 0, 12*x^2, 6*x*y, 2*y^2, 0, 0, 20*x^3, 12*x^2*y, 6*x*y^2, 2*y^3, 0, 0, 30*x^4, 20*x^3*y, 12*x^2*y^2, 6*x*y^3, 2*y^4, 0, 0, 42*x^5, 30*x^4*y, 20*x^3*y^2, 12*x^2*y^3, 6*x*y^4, 2*y^5, 0, 0, 56*x^6, 42*x^5*y, 30*x^4*y^2, 20*x^3*y^3, 12*x^2*y^4, 6*x*y^5, 2*y^6, 0, 0, 72*x^7, 56*x^6*y, 42*x^5*y^2, 30*x^4*y^3, 20*x^3*y^4, 12*x^2*y^5, 6*x*y^6, 2*y^7, 0, 0, 90*x^8, 72*x^7*y, 56*x^6*y^2, 42*x^5*y^3, 30*x^4*y^4, 20*x^3*y^5, 12*x^2*y^6, 6*x*y^7, 2*y^8, 0, 0];
    elseif derivative_degree_x==0&&derivative_degree_y==2
        r=[ 0, 0, 0, 0, 0, 2, 0, 0, 2*x, 6*y, 0, 0, 2*x^2, 6*x*y, 12*y^2, 0, 0, 2*x^3, 6*x^2*y, 12*x*y^2, 20*y^3, 0, 0, 2*x^4, 6*x^3*y, 12*x^2*y^2, 20*x*y^3, 30*y^4, 0, 0, 2*x^5, 6*x^4*y, 12*x^3*y^2, 20*x^2*y^3, 30*x*y^4, 42*y^5, 0, 0, 2*x^6, 6*x^5*y, 12*x^4*y^2, 20*x^3*y^3, 30*x^2*y^4, 42*x*y^5, 56*y^6, 0, 0, 2*x^7, 6*x^6*y, 12*x^5*y^2, 20*x^4*y^3, 30*x^3*y^4, 42*x^2*y^5, 56*x*y^6, 72*y^7, 0, 0, 2*x^8, 6*x^7*y, 12*x^6*y^2, 20*x^5*y^3, 30*x^4*y^4, 42*x^3*y^5, 56*x^2*y^6, 72*x*y^7, 90*y^8];
    end
end

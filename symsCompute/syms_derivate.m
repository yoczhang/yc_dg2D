function syms_derivate()
clc;
clearvars;
syms X1 X2 PI x y muu K k t

u_case = 1;
switch u_case
    case 1        
        phi = @(t,x,y) (exp(y)-exp(-y)).*sin(x)*exp(t);
        u1 = @(t,x,y) (K/PI*sin(2*PI*y).*cos(x))*exp(t);
        u2 = @(t,x,y) ((-2*K + K/PI^2*(sin(PI*y)).^2).*sin(x))*exp(t);
        p = @(t,x,y) 0.*x;
        
    case 2
        phi = @(x,y) y.*sin(PI*x);
        u1 = @(x,y) -sin(PI*y).*cos(PI*x);
        u2 = @(x,y) sin(PI*x).*cos(PI*y);
        p = @(x,y) sin(PI*x);
        
    case 3 
        phi = @(x,y) 0.*x;
        
        u1 = @(x,y) x.^2.*(x-1).^2.*y.*(y-1).*(2.*y-1);
        u2 = @(x,y) -x.*(x-1).*(2.*x-1).*y.^2.*(y-1).^2;
        p = @(x,y) (2.*x-1).*(2.*y-1);
        
    case 4
        u1 = @(x,y) cos(2*pi*x)*sin(2*pi*y);
        u2 = @(x,y) -sin(2*pi*x)*cos(2*pi*y);
        p = @(x,y) 0.*x;
        
    case 5
        u1 = @(x,y) 2*pi*sin(pi*x)^2*sin(pi*y)*cos(pi*y);
        u2 = @(x,y) -2*pi*sin(pi*x)*sin(pi*y)^2*cos(pi*x);
        p = @(x,y) cos(pi*x)*cos(pi*y);
        
    case 6
        u1 = @(x,y) -x^2*y*(x-1)*(3*y-2);
        u2 = @(x,y) x*y^2*(y-1)*(3*x-2);
        p = @(x,y) (2*x-1)*(2*y-1);
    
    case 7
        phi = @(x,y) -pi/4*cos(pi*x/2).*y;
        u1 = @(x,y) -cos(pi*y/2).^2.*sin(pi*x/2);
        u2 = @(x,y) 1/4*cos(pi*x/2).*(sin(pi*y)+pi*y);
        p = @(x,y) -pi/4*cos(pi*x/2).*(y-2*cos(pi*y/2).^2);
    
    case 8
        phi = @(x,y) (2-pi*sin(pi*x)).*(1-y-cos(pi*y));
        u1 = @(x,y) x.^2.*(y-1).^2+y;
        u2 = @(x,y) -2/3*x.*(y-1).^3+2-pi*sin(pi*x);
        p = @(x,y) (2-pi*sin(pi*x)).*sin(0.5*pi*y);
        
    case 9
        phi = @(x,y) 2*muu*x+(x.*(1-x).*(y-1)+1/3*y.^3-y.^2+y);
        u1 = @(x,y) y.^2-2*y+1;
        u2 = @(x,y) x.^2-x;
        p = @(x,y) 2*muu*(x+y-1)+1/(3*K);
        
end % switch

phit = diff(phi,t)
u1t = diff(u1,t)
u2t = diff(u2,t)

phix=diff(phi,x)
phiy=diff(phi,y)
phixy=diff(phix,y)
phixx=diff(phix,x)
phiyy=diff(phiy,y)

uf1 = func2str(u1)
u1x=diff(u1,x)
u1y=diff(u1,y)
u1xy=diff(u1x,y)
u1xx = diff(u1x,x)
u1yy = diff(u1y,y)


uf2 = func2str(u2)
u2x=diff(u2,x)
u2y=diff(u2,y)
u2xy=diff(u2x,y)
u2xx = diff(u2x,x)
u2yy = diff(u2y,y)

pf = func2str(p)
px=diff(p,x)
py=diff(p,y)

end % function
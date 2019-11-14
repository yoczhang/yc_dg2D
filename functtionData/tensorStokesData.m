function pde = tensorStokesData(coeff_case, u_case)
%
%   stokesData data for stokes equation
%
%
%	YcZhang 13/8/2017
%
%   Last modified 15/8/2017
%

switch coeff_case
    case 1
        mu = 1e-0;
        k_11 = @(x,y) 1+0.*x;
%         k_11x = @(x,y) 0.*x; % \partial_x kappa
%         k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1+0.*x;
%         k_22x = @(x,y) 0.*x; % \partial_x kappa
%         k_22y = @(x,y) 0.*x; % \partial_y kappa
    case 0
        mu = 1;
        k_11 = @(x,y) 1+0.*x;
%         k_11x = @(x,y) 0.*x; % \partial_x kappa
%         k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1+0.*x;
%         k_22x = @(x,y) 0.*x; % \partial_x kappa
%         k_22y = @(x,y) 0.*x; % \partial_y kappa
end

switch u_case
    case 1
        u1 = @(x,y) 2*x+3*y;
        u1x = @(x,y) 2+0.*x;
        u1y = @(x,y) 3+0.*x;
        u1xx = @(x,y) 0.*x;
        u1yy = @(x,y) 0.*x;
    case 2
        u1 = @(x,y) x.^2+y.*x;
        u1x = @(x,y) 2*x + y;
        u1y = @(x,y) x + 0.*x;
        u1xx = @(x,y) 2+0.*x;
        u1yy = @(x,y) 0.*x;
    case 3
        u1 = @(x,y) 2*x.^2.*y;
        u1x = @(x,y) 4*x.*y;
        u1y = @(x,y) 2*x.^2;
        u1xx = @(x,y) 4*y;
        u1yy = @(x,y) 0.*x;
    case 4
        u1 = @(x,y) sin(x.*y);
        u1x = @(x,y) y.*cos(x.*y);
        u1y = @(x,y) x.*cos(x.*y);
        u1xx = @(x,y) -y.^2.*sin(x.*y);
        u1yy = @(x,y) -x.^2.*sin(x.*y);
    case 5
        u1 = @(x,y) 10*x.^2.*(1-x).^2.*y.^2.*(1-y).^2;
        u1x = @(x,y) 20.*x.*y.^2.*(x - 1).^2.*(y - 1).^2 + 10.*x.^2.*y.^2.*(2.*x - 2).*(y - 1).^2;
        u1y = @(x,y) 20.*x.^2.*y.*(x - 1).^2.*(y - 1).^2 + 10.*x.^2.*y.^2.*(2.*y - 2).*(x - 1).^2;
        u1xx = @(x,y) 20.*x.^2.*y.^2.*(y - 1).^2 + 20.*y.^2.*(x - 1).^2.*(y - 1).^2 ...
            + 40.*x.*y.^2.*(2.*x - 2).*(y - 1).^2;
        u1yy = @(x,y) 20.*x.^2.*y.^2.*(x - 1).^2 + 20.*x.^2.*(x - 1).^2.*(y - 1).^2 ...
            + 40.*x.^2.*y.*(2.*y - 2).*(x - 1).^2;
    case 6
        u1 = @(x,y) exp(-x-y.^2);
        u1x = @(x,y) -exp(- y.^2 - x);
        u1y = @(x,y) -2*y.*exp(- y.^2 - x);
        u1xx = @(x,y) exp(- y.^2 - x);
        u1yy = @(x,y) 4*y.^2.*exp(- y.^2 - x) - 2*exp(- y.^2 - x);
    case 7
        u1 = @(x,y) x.^2.*(x-1).^2.*y.*(y-1).*(2.*y-1);
        u1x = @(x,y) 2.*x.*y.*(2.*y - 1).*(x - 1).^2.*(y - 1) + x.^2.*y.*(2.*x - 2).*(2.*y - 1).*(y - 1);
        u1y = @(x,y) x.^2.*y.*(2.*y - 1).*(x - 1).^2 + x.^2.*(2.*y - 1).*(x - 1).^2.*(y - 1) + 2.*x.^2.*y.*(x - 1).^2.*(y - 1);
        u1xy = @(x,y) x.^2.*y.*(2.*x - 2).*(2.*y - 1) + 2.*x.*(2.*y - 1).*(x - 1).^2.*(y - 1) + 4.*x.*y.*(x - 1).^2.*(y - 1) + x.^2.*(2.*x - 2).*(2.*y - 1).*(y - 1) + 2.*x.*y.*(2.*y - 1).*(x - 1).^2 + 2.*x.^2.*y.*(2.*x - 2).*(y - 1);
        u1xx = @(x,y) 2.*y.*(2.*y - 1).*(x - 1).^2.*(y - 1) + 2.*x.^2.*y.*(2.*y - 1).*(y - 1) + 4.*x.*y.*(2.*x - 2).*(2.*y - 1).*(y - 1);
        u1yy = @(x,y) 4.*x.^2.*(x - 1).^2.*(y - 1) + 2.*x.^2.*(2.*y - 1).*(x - 1).^2 + 4.*x.^2.*y.*(x - 1).^2;
        
        u2 = @(x,y) -x.*(x-1).*(2.*x-1).*y.^2.*(y-1).^2;
        u2x = @(x,y) - x.*y.^2.*(2.*x - 1).*(y - 1).^2 - y.^2.*(2.*x - 1).*(x - 1).*(y - 1).^2 - 2.*x.*y.^2.*(x - 1).*(y - 1).^2;
        u2y = @(x,y) - 2.*x.*y.*(2.*x - 1).*(x - 1).*(y - 1).^2 - x.*y.^2.*(2.*x - 1).*(2.*y - 2).*(x - 1);
        u2xy = @(x,y) - x.*y.^2.*(2.*x - 1).*(2.*y - 2) - 2.*y.*(2.*x - 1).*(x - 1).*(y - 1).^2 - 4.*x.*y.*(x - 1).*(y - 1).^2 - y.^2.*(2.*x - 1).*(2.*y - 2).*(x - 1) - 2.*x.*y.*(2.*x - 1).*(y - 1).^2 - 2.*x.*y.^2.*(2.*y - 2).*(x - 1);
        u2xx = @(x,y) - 4.*y.^2.*(x - 1).*(y - 1).^2 - 2.*y.^2.*(2.*x - 1).*(y - 1).^2 - 4.*x.*y.^2.*(y - 1).^2;
        u2yy = @(x,y) - 2.*x.*(2.*x - 1).*(x - 1).*(y - 1).^2 - 2.*x.*y.^2.*(2.*x - 1).*(x - 1) - 4.*x.*y.*(2.*x - 1).*(2.*y - 2).*(x - 1);
        
        p = @(x,y) (2.*x-1).*(2.*y-1);
        px = @(x,y) 2*(2*y-1);
        py = @(x,y) 2*(2*x-1);
        
    case 8
        u1 = @(x,y) -sin(pi*y).*cos(pi*x);
        u1x = @(x,y) pi.*sin(pi.*x).*sin(pi.*y);
        u1y = @(x,y) -pi.*cos(pi.*x).*cos(pi.*y);
        u1xy = @(x,y) pi.^2.*cos(pi.*y).*sin(pi.*x);
        u1xx = @(x,y) pi.^2.*cos(pi.*x).*sin(pi.*y);
        u1yy = @(x,y) pi.^2.*cos(pi.*x).*sin(pi.*y);
        
        u2 = @(x,y) sin(pi*x).*cos(pi*y);
        u2x = @(x,y) pi.*cos(pi.*x).*cos(pi.*y);
        u2y = @(x,y) -pi.*sin(pi.*x).*sin(pi.*y);
        u2xy = @(x,y) -pi.^2.*cos(pi.*x).*sin(pi.*y);
        u2xx = @(x,y) -pi.^2.*cos(pi.*y).*sin(pi.*x);
        u2yy = @(x,y) -pi.^2.*cos(pi.*y).*sin(pi.*x);
        
        p = @(x,y) sin(pi*x)-2/pi;
        px = @(x,y) pi.*cos(pi.*x);
        py = @(x,y) 0.*x;
        
    case 9
        KK = 1;
        u1 = @(x,y) KK/pi.*sin(2.*pi.*y).*cos(x);
        u1x = @(x,y) -(KK.*sin(2.*pi.*y).*sin(x))/pi;
        u1y = @(x,y) 2.*KK.*cos(2.*pi.*y).*cos(x);
        u1xy = @(x,y) -2.*KK.*cos(2.*pi.*y).*sin(x);
        u1xx = @(x,y) -(KK.*sin(2.*pi.*y).*cos(x))/pi;
        u1yy = @(x,y) -4.*KK.*pi.*sin(2.*pi.*y).*cos(x);

        u2 = @(x,y) (-2.*KK+KK/pi.^2.*sin(pi.*y).*sin(pi.*y)).*sin(x);
        u2x = @(x,y) -cos(x).*(2.*KK - (KK.*sin(pi.*y).^2)/pi.^2);
        u2y = @(x,y) (2.*KK.*cos(pi.*y).*sin(pi.*y).*sin(x))/pi;
        u2xy = @(x,y) (2.*KK.*cos(pi.*y).*sin(pi.*y).*cos(x))/pi;
        u2xx = @(x,y) sin(x).*(2.*KK - (KK.*sin(pi.*y).^2)/pi.^2);
        u2yy = @(x,y) 2.*KK.*cos(pi.*y).^2.*sin(x) - 2.*KK.*sin(pi.*y).^2.*sin(x);

        p = @(x,y) 0.*x;
        px = @(x,y) 0.*x; 
        py = @(x,y) 0.*x; 
end


f1 = @(x,y) -mu.*(2*u1xx(x,y)+u1yy(x,y)+u2xy(x,y))+px(x,y);
f2 = @(x,y) -mu.*(u1xy(x,y)+u2xx(x,y)+2*u2yy(x,y))+py(x,y);

gD1 = u1;
gD2 = u2;

funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1 + 0.*x;

pde = struct(...
    'mu', mu,  ...
    'k_11', k_11, ...
    'k_22', k_22, ...
    'u1',u1, ...
    'u1x',u1x, ...
    'u1y',u1y, ...
    'u2',u2, ...
    'u2x',u2x, ...
    'u2y',u2y, ...
    'p', p, ...
    'px', px, ...
    'py', py, ...
    'f1',f1, ...
    'f2',f2, ...
    'gD1',gD1, ... % Dirichlet function
    'gD2',gD2, ... % Dirichlet function
    'funcZero',funcZero, ...
    'funcOne',funcOne ...
    );

end % pde function

function pde = poissonData(coeff_case, u_case)
%
%   poissonData data for poisson equation
%
%
%	YcZhang 13/5/2017
%
%   Last modified 14/5/2017
%

switch coeff_case
    case 1
        k_11 = @(x,y) 2+sin(x).*sin(y);
        k_11x = @(x,y) cos(x).*sin(y); % \partial_x kappa
        k_11y = @(x,y) sin(x).*cos(y); % \partial_y kappa
        k_22 = @(x,y) 2+sin(x).*sin(y);
        k_22x = @(x,y) cos(x).*sin(y); % \partial_x kappa
        k_22y = @(x,y) sin(x).*cos(y); % \partial_y kappa
    case 0
        k_11 = @(x,y) 1+0.*x;
        k_11x = @(x,y) 0.*x; % \partial_x kappa
        k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1+0.*x;
        k_22x = @(x,y) 0.*x; % \partial_x kappa
        k_22y = @(x,y) 0.*x; % \partial_y kappa
end

switch u_case
    case 1
        u = @(x,y) 2*x+3*y;
        ux = @(x,y) 2+0.*x;
        uy = @(x,y) 3+0.*x;
        uxx = @(x,y) 0.*x;
        uyy = @(x,y) 0.*x;
    case 2
        u = @(x,y) x.^2+y.*x;
        ux = @(x,y) 2*x + y;
        uy = @(x,y) x + 0.*x;
        uxx = @(x,y) 2+0.*x;
        uyy = @(x,y) 0.*x;
    case 3
        u = @(x,y) 2*x.^2.*y;
        ux = @(x,y) 4*x.*y;
        uy = @(x,y) 2*x.^2;
        uxx = @(x,y) 4*y;
        uyy = @(x,y) 0.*x;
    case 4
        u = @(x,y) sin(x.*y);
        ux = @(x,y) y.*cos(x.*y);
        uy = @(x,y) x.*cos(x.*y);
        uxx = @(x,y) -y.^2.*sin(x.*y);
        uyy = @(x,y) -x.^2.*sin(x.*y);
    case 5
        u = @(x,y) 10*x.^2.*(1-x).^2.*y.^2.*(1-y).^2;
        ux = @(x,y) 20.*x.*y.^2.*(x - 1).^2.*(y - 1).^2 + 10.*x.^2.*y.^2.*(2.*x - 2).*(y - 1).^2;
        uy = @(x,y) 20.*x.^2.*y.*(x - 1).^2.*(y - 1).^2 + 10.*x.^2.*y.^2.*(2.*y - 2).*(x - 1).^2;
        uxx = @(x,y) 20.*x.^2.*y.^2.*(y - 1).^2 + 20.*y.^2.*(x - 1).^2.*(y - 1).^2 ...
            + 40.*x.*y.^2.*(2.*x - 2).*(y - 1).^2;
        uyy = @(x,y) 20.*x.^2.*y.^2.*(x - 1).^2 + 20.*x.^2.*(x - 1).^2.*(y - 1).^2 ...
            + 40.*x.^2.*y.*(2.*y - 2).*(x - 1).^2;
    case 6
        u = @(x,y) exp(-x-y.^2);
        ux = @(x,y) -exp(- y.^2 - x);
        uy = @(x,y) -2*y.*exp(- y.^2 - x);
        uxx = @(x,y) exp(- y.^2 - x);
        uyy = @(x,y) 4*y.^2.*exp(- y.^2 - x) - 2*exp(- y.^2 - x);
end


f = @(x,y) -(k_11x(x,y).*ux(x,y)+k_11(x,y).*uxx(x,y))...
             -(k_22y(x,y).*uy(x,y)+k_22(x,y).*uyy(x,y));

gD = u;

func0 = @(x,y) 0.*x;

pde = struct(...
    'k_11',k_11, 'k_22', k_22, ...
    'u',u, ...
    'ux',ux, ...
    'uy',uy, ...
    'f',f, ...
    'gD',gD, ... % Dirichlet function
    'func0',func0 ...
    );

end % pde function

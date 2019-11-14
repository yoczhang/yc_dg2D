function pde = parabolicData(coeff_case, u_case)
%
%   parabolicData data for poisson equation
%
%
%	YcZhang 3/9/2017
%
%   Last modified 3/9/2017
%

switch coeff_case
    case 0
        k_11 = @(x,y) 1+0.*x;
        k_11x = @(x,y) 0.*x; % \partial_x kappa
        k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1+0.*x;
        k_22x = @(x,y) 0.*x; % \partial_x kappa
        k_22y = @(x,y) 0.*x; % \partial_y kappa
    case 1
        k_11 = @(x,y) 2+sin(x).*sin(y);
        k_11x = @(x,y) cos(x).*sin(y); % \partial_x kappa
        k_11y = @(x,y) sin(x).*cos(y); % \partial_y kappa
        k_22 = @(x,y) 2+sin(x).*sin(y);
        k_22x = @(x,y) cos(x).*sin(y); % \partial_x kappa
        k_22y = @(x,y) sin(x).*cos(y); % \partial_y kappa
end


switch u_case
    case 1
        u = @(t,x,y) x.*(1-x).*y.*(1-y).*exp(t);
        ut = @(t,x,y) x.*y.*exp(t).*(x - 1).*(y - 1);
        ux = @(t,x,y) y.*exp(t).*(2.*x - 1).*(y - 1);
        uy = @(t,x,y) x.*exp(t).*(2.*y - 1).*(x - 1);
        uxx = @(t,x,y) 2*y.*exp(t).*(y - 1);
        uyy = @(t,x,y) 2*x.*exp(t).*(x - 1);
        
    case 2
        u = @(t,x,y) (exp(y)-exp(-y)).*sin(x).*exp(t);
        ut = @(t,x,y) -exp(t).*sin(x).*(exp(-y) - exp(y));
        ux = @(t,x,y) -exp(t).*cos(x).*(exp(-y) - exp(y));
        uy = @(t,x,y) exp(t).*sin(x).*(exp(-y) + exp(y));
        uxy = @(t,x,y) exp(t).*cos(x).*(exp(-y) + exp(y));
        uxx = @(t,x,y) exp(t).*sin(x).*(exp(-y) - exp(y));
        uyy = @(t,x,y) -exp(t).*sin(x).*(exp(-y) - exp(y));
        
    case 3
        u = @(t,x,y) x.*(1-x).*y.*(1-y).*t;
        ut = @(t,x,y) x.*y.*(x - 1).*(y - 1);
        ux = @(t,x,y) t.*x.*y.*(y - 1) + t.*y.*(x - 1).*(y - 1);
        uy = @(t,x,y) t.*x.*y.*(x - 1) + t.*x.*(x - 1).*(y - 1);
        uxy = @(t,x,y) t.*x.*(y - 1) + t.*y.*(x - 1) + t.*(x - 1).*(y - 1) + t.*x.*y;
        uxx = @(t,x,y) 2.*t.*y.*(y - 1);
        uyy = @(t,x,y) 2.*t.*x.*(x - 1);
        
    case 4
        alpha = 3;
        beta = 1.2;
        
        u = @(t,x,y) 1+x.^2 + alpha*y.^2 + beta*t;
        ut = @(t,x,y) beta + 0.*x;
        ux = @(t,x,y) 2*x;
        uy = @(t,x,y) 2*alpha*y;
        uxy = @(t,x,y) 0.*x;
        uxx = @(t,x,y) 2 + 0.*x;
        uyy = @(t,x,y) 2*alpha + 0.*x;
        
end 


f = @(t,x,y) ut(t,x,y) -(k_11x(x,y).*ux(t,x,y)+k_11(x,y).*uxx(t,x,y)) ...
    -(k_22y(x,y).*uy(t,x,y)+k_22(x,y).*uyy(t,x,y));


gD = u;
u0 =@(x,y) u(0,x,y);

funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1 + 0.*x;


pde = struct(...
    'k_11',k_11, 'k_22', k_22, ...
    'u',u, ...
    'ux',ux, ...
    'uy',uy, ...
    'f',f, ...
    'gD',gD, ... % Dirichlet function
    'u0',u0, ...
    'funcZero',funcZero, ...
    'funcOne',funcOne ...
    );


end % function
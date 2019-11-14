function pde = navierstokesData1(coeff_case, u_case)
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       Using the Saclar-Stokes, 
%       the dgNavierStokes2 using the Tensor-Stokes.
%   %---------------------------------------------
%
%   navierstokesData data for navier-stokes equation
%
%
%	YcZhang 17/10/2017
%
%   Last modified 17/10/2017
%

switch coeff_case
    case 1
        mu = @(x,y) 1e-3+0.*x;
                k_11 = @(x,y) 1+0.*x;
%         k_11x = @(x,y) 0.*x; % \partial_x kappa
%         k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1+0.*x;
%         k_22x = @(x,y) 0.*x; % \partial_x kappa
%         k_22y = @(x,y) 0.*x; % \partial_y kappa
    case 0
        mu = @(x,y) 1+0.*x;
        k_11 = @(x,y) 1+0.*x;
%         k_11x = @(x,y) 0.*x; % \partial_x kappa
%         k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1+0.*x;
%         k_22x = @(x,y) 0.*x; % \partial_x kappa
%         k_22y = @(x,y) 0.*x; % \partial_y kappa
end

switch u_case
    case 1
        u1 = @(x,y) x.^2.*(x-1).^2.*y.*(y-1).*(2.*y-1);
        u1x = @(x,y) 2.*x.*y.*(2.*y - 1).*(x - 1).^2.*(y - 1) + x.^2.*y.*(2.*x - 2).*(2.*y - 1).*(y - 1);
        u1y = @(x,y) x.^2.*y.*(2.*y - 1).*(x - 1).^2 + x.^2.*(2.*y - 1).*(x - 1).^2.*(y - 1) + 2.*x.^2.*y.*(x - 1).^2.*(y - 1);
        u1xx = @(x,y) 2.*y.*(2.*y - 1).*(x - 1).^2.*(y - 1) + 2.*x.^2.*y.*(2.*y - 1).*(y - 1) + 4.*x.*y.*(2.*x - 2).*(2.*y - 1).*(y - 1);
        u1yy = @(x,y) 4.*x.^2.*(x - 1).^2.*(y - 1) + 2.*x.^2.*(2.*y - 1).*(x - 1).^2 + 4.*x.^2.*y.*(x - 1).^2;
        
        u2 = @(x,y) -x.*(x-1).*(2.*x-1).*y.^2.*(y-1).^2;
        u2x = @(x,y) - x.*y.^2.*(2.*x - 1).*(y - 1).^2 - y.^2.*(2.*x - 1).*(x - 1).*(y - 1).^2 - 2.*x.*y.^2.*(x - 1).*(y - 1).^2;
        u2y = @(x,y) - 2.*x.*y.*(2.*x - 1).*(x - 1).*(y - 1).^2 - x.*y.^2.*(2.*x - 1).*(2.*y - 2).*(x - 1);
        u2xx = @(x,y) - 4.*y.^2.*(x - 1).*(y - 1).^2 - 2.*y.^2.*(2.*x - 1).*(y - 1).^2 - 4.*x.*y.^2.*(y - 1).^2;
        u2yy = @(x,y) - 2.*x.*(2.*x - 1).*(x - 1).*(y - 1).^2 - 2.*x.*y.^2.*(2.*x - 1).*(x - 1) - 4.*x.*y.*(2.*x - 1).*(2.*y - 2).*(x - 1);
        
        p = @(x,y) (2.*x-1).*(2.*y-1);
        px = @(x,y) 2*(2*y-1);
        py = @(x,y) 2*(2*x-1);
        
        f1 = @(x,y) -mu(x,y).*(u1xx(x,y)+u1yy(x,y))+(u1(x,y).*u1x(x,y)+u2(x,y).*u1y(x,y))+px(x,y);
        f2 = @(x,y) -mu(x,y).*(u2xx(x,y)+u2yy(x,y))+(u1(x,y).*u2x(x,y)+u2(x,y).*u2y(x,y))+py(x,y);

        
    case 2 
        u1 = @(x,y) cos(2.*pi.*x).*sin(2.*pi.*y);
        u1x = @(x,y) -2.*pi.*sin(2.*pi.*x).*sin(2.*pi.*y);
        u1y = @(x,y) 2.*pi.*cos(2.*pi.*x).*cos(2.*pi.*y);
        u1xx = @(x,y) -4.*pi.^2.*cos(2.*pi.*x).*sin(2.*pi.*y);
        u1yy = @(x,y) -4.*pi.^2.*cos(2.*pi.*x).*sin(2.*pi.*y);
        
        u2 = @(x,y) -sin(2.*pi.*x).*cos(2.*pi.*y);
        u2x = @(x,y) -2.*pi.*cos(2.*pi.*x).*cos(2.*pi.*y);
        u2y = @(x,y) 2.*pi.*sin(2.*pi.*x).*sin(2.*pi.*y);
        u2xx = @(x,y) 4.*pi.^2.*cos(2.*pi.*y).*sin(2.*pi.*x);
        u2yy = @(x,y) 4.*pi.^2.*cos(2.*pi.*y).*sin(2.*pi.*x);
        
        p = @(x,y) 0.*x;
        px = @(x,y) 0.*x;
        py = @(x,y) 0.*x;
        
        f1 = @(x,y) -mu(x,y).*(u1xx(x,y)+u1yy(x,y))+(u1(x,y).*u1x(x,y)+u2(x,y).*u1y(x,y))+px(x,y);
        f2 = @(x,y) -mu(x,y).*(u2xx(x,y)+u2yy(x,y))+(u1(x,y).*u2x(x,y)+u2(x,y).*u2y(x,y))+py(x,y);
        
    case 3
        u1 = @(x,y) 2.*pi.*sin(pi.*x).^2.*sin(pi.*y).*cos(pi.*y);
        u1x = @(x,y) 4.*pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*x).*sin(pi.*y);
        u1y = @(x,y) 2.*pi.^2.*cos(pi.*y).^2.*sin(pi.*x).^2 - 2.*pi.^2.*sin(pi.*x).^2.*sin(pi.*y).^2;
        u1xx = @(x,y) 4.*pi.^3.*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y) - 4.*pi.^3.*cos(pi.*y).*sin(pi.*x).^2.*sin(pi.*y);
        u1yy = @(x,y) -8.*pi.^3.*cos(pi.*y).*sin(pi.*x).^2.*sin(pi.*y);
        
        u2 = @(x,y) -2.*pi.*sin(pi.*x).*sin(pi.*y).^2.*cos(pi.*x);
        u2x = @(x,y) 2.*pi.^2.*sin(pi.*x).^2.*sin(pi.*y).^2 - 2.*pi.^2.*cos(pi.*x).^2.*sin(pi.*y).^2;
        u2y = @(x,y) -4.*pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*x).*sin(pi.*y);
        u2xx = @(x,y) 8.*pi.^3.*cos(pi.*x).*sin(pi.*x).*sin(pi.*y).^2;
        u2yy = @(x,y) 4.*pi.^3.*cos(pi.*x).*sin(pi.*x).*sin(pi.*y).^2 - 4.*pi.^3.*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x);
        
        p = @(x,y) cos(pi.*x).*cos(pi.*y);
        px = @(x,y) -pi.*cos(pi.*y).*sin(pi.*x);
        py = @(x,y) -pi.*cos(pi.*x).*sin(pi.*y);
        
        f1 = @(x,y) -mu(x,y).*(u1xx(x,y)+u1yy(x,y))+(u1(x,y).*u1x(x,y)+u2(x,y).*u1y(x,y))+px(x,y);
        f2 = @(x,y) -mu(x,y).*(u2xx(x,y)+u2yy(x,y))+(u1(x,y).*u2x(x,y)+u2(x,y).*u2y(x,y))+py(x,y);
        
    case 4 % NS Step Channel Problem
        u1 = @(x,y) velocity_u1(x,y);
        u2 = @(x,y) velocity_u2(x,y);
        p = @(x,y) 0.*x;

        f1 = @(x,y) 0.*x;
        f2 = @(x,y) 0.*x;
        
        %---- the following is just to get the structure: pde
        u1x = @(x,y) 0.*x;
        u1y = @(x,y) 0.*x;
        u2x = @(x,y) 0.*x;
        u2y = @(x,y) 0.*x;
        px = @(x,y) 0.*x;
        py = @(x,y) 0.*x;
        
    case 5
        u1 = @(x,y) -x.^2.*y.*(x-1).*(3.*y-2);
        u1x = @(x,y) - x.^2.*y.*(3.*y - 2) - 2.*x.*y.*(3.*y - 2).*(x - 1);
        u1y = @(x,y) - x.^2.*(3.*y - 2).*(x - 1) - 3.*x.^2.*y.*(x - 1);
        u1xx = @(x,y) - 4.*x.*y.*(3.*y - 2) - 2.*y.*(3.*y - 2).*(x - 1);
        u1yy = @(x,y) -6.*x.^2.*(x - 1);
        
        u2 = @(x,y) x.*y.^2.*(y-1).*(3.*x-2);
        u2x = @(x,y) y.^2.*(3.*x - 2).*(y - 1) + 3.*x.*y.^2.*(y - 1);
        u2y = @(x,y) x.*y.^2.*(3.*x - 2) + 2.*x.*y.*(3.*x - 2).*(y - 1);
        u2xx = @(x,y) 6.*y.^2.*(y - 1);
        u2yy = @(x,y) 4.*x.*y.*(3.*x - 2) + 2.*x.*(3.*x - 2).*(y - 1);
        
        p = @(x,y) (2.*x-1).*(2.*y-1);
        px = @(x,y) 4.*y - 2;
        py = @(x,y) 4.*x - 2;
        
        f1 = @(x,y) -mu(x,y).*(u1xx(x,y)+u1yy(x,y))+(u1(x,y).*u1x(x,y)+u2(x,y).*u1y(x,y))+px(x,y);
        f2 = @(x,y) -mu(x,y).*(u2xx(x,y)+u2yy(x,y))+(u1(x,y).*u2x(x,y)+u2(x,y).*u2y(x,y))+py(x,y);

end



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

%%
%%>> -- Begin sub function  -------------------------------------------------------------
function r = velocity_u1(x,y)

if abs(x-0)<=5e-8
    r=y.*(1-y);
elseif abs(x-10)<=5e-8
    r=y.*(1-y);
else
    r=0.*x;
end

end % function


function r = velocity_u2(x,y)

r=0.*x;

end % function

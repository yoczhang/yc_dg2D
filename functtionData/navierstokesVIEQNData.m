function pde = navierstokesVIEQNData(coeff_case, u_case)
%
%   navierstokesData data for navier-stokes equation
%
%
%	YcZhang 17/10/2017
%
%   Last modified 17/10/2017
%

switch coeff_case
    case 0
        mu = 1;
        k_11 = @(x,y) 1+0.*x;
%         k_11x = @(x,y) 0.*x; % \partial_x kappa
%         k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1+0.*x;
%         k_22x = @(x,y) 0.*x; % \partial_x kappa
%         k_22y = @(x,y) 0.*x; % \partial_y kappa

    case 1
        mu = 1e-3;
        k_11 = @(x,y) 1+0.*x;
%         k_11x = @(x,y) 0.*x; % \partial_x kappa
%         k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1+0.*x;
%         k_22x = @(x,y) 0.*x; % \partial_x kappa
%         k_22y = @(x,y) 0.*x; % \partial_y kappa
end

switch u_case
    case 1
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
        
        f1 = @(x,y) -mu.*(u1xx(x,y)+u1yy(x,y))+(u1(x,y).*u1x(x,y)+u2(x,y).*u1y(x,y))+px(x,y);
        f2 = @(x,y) -mu.*(u2xx(x,y)+u2yy(x,y))+(u1(x,y).*u2x(x,y)+u2(x,y).*u2y(x,y))+py(x,y);

        
        
end % swtich



gD1 = u1;
gD2 = u2;
fric_g1 = @(x,y) 4*mu^2.*x.^2.*(1-x);
fric_g2 = @(x,y) 4*mu^2.*y.^2.*(1-y);



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
    'fric_g1',fric_g1, ...
    'fric_g2',fric_g2, ...
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

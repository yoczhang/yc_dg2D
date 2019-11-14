function pde = stokesData(coeff_case, u_case)
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
        pdecase = 'case1';
        u1 = @(x,y) y.^2-2*y+1;
        u1x = @(x,y) 0.*x;
        u1y = @(x,y) 2.*y - 2;
        u1xy = @(x,y) 0.*x;
        u1xx = @(x,y) 0.*x;
        u1yy = @(x,y) 2+0.*x;
        
        u2 = @(x,y) x.^2-x;
        u2x = @(x,y) 2.*x - 1;
        u2y = @(x,y) 0.*x;
        u2xy = @(x,y) 0.*x;
        u2xx = @(x,y) 2+0.*x;
        u2yy = @(x,y) 0.*x;
        
        p = @(x,y) 2*(x+y-1)+1/(3*1);
        px = @(x,y) 2+0.*x;
        py = @(x,y) 2+0.*x;
    case 2
        pdecase = 'case2';
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
        
    case 3
        pdecase = 'case3';
        
        u1 = @(x,y) -cos(pi*y/2).^2.*sin(pi*x/2);
        u1x = @(x,y) -(pi.*cos((pi.*x)/2).*cos((pi.*y)/2).^2)/2;
        u1y = @(x,y) pi.*cos((pi.*y)/2).*sin((pi.*x)/2).*sin((pi.*y)/2);
        u1xy = @(x,y) (pi.^2.*cos((pi.*x)/2).*cos((pi.*y)/2).*sin((pi.*y)/2))/2;
        u1xx = @(x,y) (pi.^2.*cos((pi.*y)/2).^2.*sin((pi.*x)/2))/4;
        u1yy = @(x,y) (pi.^2.*cos((pi.*y)/2).^2.*sin((pi.*x)/2))/2 - (pi.^2.*sin((pi.*x)/2).*sin((pi.*y)/2).^2)/2;
        
        u2 = @(x,y) 1/4*cos(pi*x/2).*(sin(pi*y)+pi*y);
        u2x = @(x,y) -(pi.*sin((pi.*x)/2).*(sin(pi.*y) + pi.*y))/8;
        u2y = @(x,y) (cos((pi.*x)/2).*(pi + pi.*cos(pi.*y)))/4;
        u2xy = @(x,y) -(pi.*sin((pi.*x)/2).*(pi + pi.*cos(pi.*y)))/8;
        u2xx = @(x,y) -(pi.^2.*cos((pi.*x)/2).*(sin(pi.*y) + pi.*y))/16;
        u2yy = @(x,y) -(pi.^2.*cos((pi.*x)/2).*sin(pi.*y))/4;
        
        p = @(x,y) -pi/4*cos(pi*x/2).*(y-2*cos(pi*y/2).^2)+1/4;
        px = @(x,y) (pi.^2.*sin((pi.*x)/2).*(y - 2.*cos((pi.*y)/2).^2))/8;
        py = @(x,y) -(pi.*cos((pi.*x)/2).*(2.*pi.*cos((pi.*y)/2).*sin((pi.*y)/2) + 1))/4;
end


f1 = @(x,y) -mu(x,y).*(u1xx(x,y)+u1yy(x,y))+px(x,y);
f2 = @(x,y) -mu(x,y).*(u2xx(x,y)+u2yy(x,y))+py(x,y);

gD1 = u1;
gD2 = u2;

funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1 + 0.*x;

pde = struct(...
    'mu', mu,  ...
    'k_11', k_11, ...
    'k_22', k_22, ...
    'pdecase', pdecase, ...
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

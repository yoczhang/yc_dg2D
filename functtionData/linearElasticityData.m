function pde = linearElasticityData(coeff_case, u_case)
%
%   stokesData data for linearElasticity equation
%
%
%	YcZhang 22/9/2017
%
%   Last modified 22/9/2017
%

% switch coeff_case
%     case 1
%         mu = @(x,y) 1e-3+0.*x;
%                 k_11 = @(x,y) 1+0.*x;
% %         k_11x = @(x,y) 0.*x; % \partial_x kappa
% %         k_11y = @(x,y) 0.*x; % \partial_y kappa
%         k_22 = @(x,y) 1+0.*x;
% %         k_22x = @(x,y) 0.*x; % \partial_x kappa
% %         k_22y = @(x,y) 0.*x; % \partial_y kappa
%     case 0
%         mu = @(x,y) 1+0.*x;
%         k_11 = @(x,y) 1+0.*x;
% %         k_11x = @(x,y) 0.*x; % \partial_x kappa
% %         k_11y = @(x,y) 0.*x; % \partial_y kappa
%         k_22 = @(x,y) 1+0.*x;
% %         k_22x = @(x,y) 0.*x; % \partial_x kappa
% %         k_22y = @(x,y) 0.*x; % \partial_y kappa
% end

switch coeff_case
    case 1
        lambda = 1; 
        mu = 1;    
    
end 

switch u_case
    case 1
        u1 = @(x,y) sin(pi*x).*sin(pi*y);
        u1x = @(x,y) pi.*cos(pi.*x).*sin(pi.*y);
        u1y = @(x,y) pi.*cos(pi.*y).*sin(pi.*x);
        u1xy = @(x,y) pi^2*cos(pi*x).*cos(pi*y);
        u1xx = @(x,y) -pi.^2.*sin(pi.*x).*sin(pi.*y);
        u1yy = @(x,y) -pi.^2.*sin(pi.*x).*sin(pi.*y);
        
        u2 = @(x,y) x.*(x-1).*y.*(y-1);
        u2x = @(x,y) x.*y.*(y - 1) + y.*(x - 1).*(y - 1);
        u2y = @(x,y) x.*y.*(x - 1) + x.*(x - 1).*(y - 1);
        u2xy = @(x,y) (x - 1).*(y - 1) + x.*y + x.*(y - 1) + y.*(x - 1);
        u2xx = @(x,y) 2.*y.*(y - 1);
        u2yy = @(x,y) 2.*x.*(x - 1);
        
end


f1 = @(x,y) -(lambda*( u1xx(x,y)+u2xy(x,y) )+ 2*mu*u1xx(x,y) + mu*(u1yy(x,y)+u2xy(x,y)));
f2 = @(x,y) -(mu*(u1xy(x,y) + u2xx(x,y) ) + lambda*(u1xy(x,y) + u2yy(x,y)) +2*mu*u2yy(x,y) );

gD1 = u1;
gD2 = u2;

funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1 + 0.*x;

pde = struct(...
    'lambda', lambda, ...
    'mu', mu,  ...
    'u1',u1, ...
    'u1x',u1x, ...
    'u1y',u1y, ...
    'u2',u2, ...
    'u2x',u2x, ...
    'u2y',u2y, ...
    'f1',f1, ...
    'f2',f2, ...
    'gD1',gD1, ... % Dirichlet function
    'gD2',gD2, ... % Dirichlet function
    'funcZero',funcZero, ...
    'funcOne',funcOne ...
    );

end % pde function

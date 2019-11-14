function pde = StokesDarcyData(u_case)
%
%   StokesDarcyData data for stokes equation
%
%
%	YcZhang 27/8/2017
%
%   Last modified 27/8/2017
%

switch u_case
    case 1
        %---- Compute domain:
        % Stokes: [0,1]x[0,1]; Darcy: [0,1]x[-1,0];
        %
        mu = 1e-2; % fluid kinematic viscosity
        K = 1e-2; % Darcy permeability matrix
        g = 1; % gravitational acceleration
        Kappa = 1; % this is the coefficient of (P(uh),P(vh))_interface
        pdecase = 'case1';
        
        phi = @(x,y) (exp(y)-exp(-y)).*sin(x);
        phix = @(x,y) -cos(x).*(exp(-y) - exp(y));
        phiy = @(x,y) sin(x).*(exp(-y) + exp(y));
        phixy = @(x,y) cos(x).*(exp(-y) + exp(y));
        phixx = @(x,y) sin(x).*(exp(-y) - exp(y));
        phiyy = @(x,y) -sin(x).*(exp(-y) - exp(y));

        u1 = @(x,y) K/pi.*sin(2.*pi.*y).*cos(x);
        u1x = @(x,y) -(K.*sin(2.*pi.*y).*sin(x))/pi;
        u1y = @(x,y) 2.*K.*cos(2.*pi.*y).*cos(x);
        u1xy = @(x,y) -2.*K.*cos(2.*pi.*y).*sin(x);
        u1xx = @(x,y) -(K.*sin(2.*pi.*y).*cos(x))/pi;
        u1yy = @(x,y) -4.*K.*pi.*sin(2.*pi.*y).*cos(x);

        u2 = @(x,y) (-2.*K+K/pi.^2.*(sin(pi.*y)).^2).*sin(x);
        u2x = @(x,y) -cos(x).*(2.*K - (K.*sin(pi.*y).^2)/pi.^2);
        u2y = @(x,y) (2.*K.*cos(pi.*y).*sin(pi.*y).*sin(x))/pi;
        u2xy = @(x,y) (2.*K.*cos(pi.*y).*sin(pi.*y).*cos(x))/pi;
        u2xx = @(x,y) sin(x).*(2.*K - (K.*sin(pi.*y).^2)/pi.^2);
        u2yy = @(x,y) 2.*K.*cos(pi.*y).^2.*sin(x) - 2.*K.*sin(pi.*y).^2.*sin(x);

        p = @(x,y) 0.*x;
        px = @(x,y) 0.*x; 
        py = @(x,y) 0.*x; 
        
    case 2
        %---- Compute domain:
        % Stokes: [0,1]x[0,1]; Darcy: [0,1]x[1,2];
        %
        mu = 1; % fluid kinematic viscosity
        K = 1; % Darcy permeability matrix
        g = 1; % gravitational acceleration
        Kappa = 1; % this is the coefficient of (P(uh),P(vh))_interface
        pdecase = 'case2';
        
        phi = @(x,y) y.*sin(pi.*x);
        phix = @(t,x,y) pi.*y.*cos(pi.*x);
        phiy = @(t,x,y) sin(pi.*x);
        phixy = @(t,x,y) pi.*cos(pi.*x);
        phixx = @(t,x,y) -pi.^2.*y.*sin(pi.*x);
        phiyy = @(t,x,y) 0.*x; 


        u1 = @(x,y) -sin(pi.*y).*cos(pi.*x);
        u1x = @(x,y) pi.*sin(pi.*x).*sin(pi.*y);
        u1y = @(x,y) -pi.*cos(pi.*x).*cos(pi.*y);
        u1xy = @(x,y) pi.^2.*cos(pi.*y).*sin(pi.*x);
        u1xx = @(x,y) pi.^2.*cos(pi.*x).*sin(pi.*y);
        u1yy = @(x,y) pi.^2.*cos(pi.*x).*sin(pi.*y);


        u2 = @(x,y) sin(pi.*x).*cos(pi.*y);
        u2x = @(x,y) pi.*cos(pi.*x).*cos(pi.*y);
        u2y = @(x,y) -pi.*sin(pi.*x).*sin(pi.*y);
        u2xy = @(x,y) -pi.^2.*cos(pi.*x).*sin(pi.*y);
        u2xx = @(x,y) -pi.^2.*cos(pi.*y).*sin(pi.*x);
        u2yy = @(x,y) -pi.^2.*cos(pi.*y).*sin(pi.*x);


        p = @(x,y) sin(pi.*x);
        px = @(x,y) pi.*cos(pi.*x);
        py = @(x,y) 0.*x;
        
    case 3
        K = @(x,y) k_function(x,y);

end

fp = @(x,y) -K*(phixx(x,y) + phiyy(x,y));
f1 = @(x,y) -mu*(2*u1xx(x,y)+u1yy(x,y)+u2xy(x,y))+px(x,y);
f2 = @(x,y) -mu*(u1xy(x,y)+u2xx(x,y)+2*u2yy(x,y))+py(x,y);

gDp = phi;
gD1 = u1;
gD2 = u2;

funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1 + 0.*x;

pde = struct(...
    'mu', mu,  ...
    'K', K, ...
    'g', g, ...
    'Kappa', Kappa, ...
    'pdecase', pdecase, ...
    'phi', phi, ...
    'phix', phix, ...
    'phiy', phiy, ...
    'u1',u1, ...
    'u1x',u1x, ...
    'u1y',u1y, ...
    'u2',u2, ...
    'u2x',u2x, ...
    'u2y',u2y, ...
    'p', p, ...
    'px', px, ...
    'py', py, ...
    'fp', fp, ...
    'f1',f1, ...
    'f2',f2, ...
    'gDp', gDp, ...
    'gD1',gD1, ... % Dirichlet function
    'gD2',gD2, ... % Dirichlet function
    'funcZero',funcZero, ...
    'funcOne',funcOne ...
    );

end % pde function

%%
%%>> -- Begin sub function  -------------------------------------------------------------
function r = k_function(x,y)

if 0<=x && x<=0.75 && 0.5<=y && y<=0.75
    r = 1e-4;
elseif 1.25<=x && x<=2 && 0.25<=y && y<=0.5
    r = 1e-6;
else
    r = 1;
    
end

end % function k_function


function pde = timeStokesDarcyData(u_case)
%
%   Non-stationary StokesDarcyData data for Stokes-Darcy equation
%
%
%	YcZhang 7/11/2017
%
%   Last modified 7/11/2017
%

switch u_case
    case 1
        %---- Compute domain:
        % Stokes: [0,1]x[0,1]; Darcy: [0,1]x[-1,0];
        %
        mu = 1e-2; % fluid kinematic viscosity
        K = 1e-2; % Darcy permeability matrix
        g = 1; % gravitational acceleration
        Kappa = 1; % this is the coefficient of (P_t(uh),P_t(vh))_interface
        pdecase = 'case1';
        
        phi = @(t,x,y) (exp(y)-exp(-y)).*sin(x).*exp(t);
        phit = @(t,x,y) -exp(t).*sin(x).*(exp(-y) - exp(y));
        phix = @(t,x,y) -exp(t).*cos(x).*(exp(-y) - exp(y));
        phiy = @(t,x,y) exp(t).*sin(x).*(exp(-y) + exp(y));
        phixy = @(t,x,y) exp(t).*cos(x).*(exp(-y) + exp(y));
        phixx = @(t,x,y) exp(t).*sin(x).*(exp(-y) - exp(y));
        phiyy = @(t,x,y) -exp(t).*sin(x).*(exp(-y) - exp(y));
        
        u1 = @(t,x,y) (K/pi.*sin(2.*pi.*y).*cos(x)).*exp(t);
        u1t = @(t,x,y) (K.*exp(t).*sin(2.*pi.*y).*cos(x))/pi;
        u1x = @(t,x,y) -(K.*exp(t).*sin(2.*pi.*y).*sin(x))/pi;
        u1y = @(t,x,y) 2.*K.*exp(t).*cos(2.*pi.*y).*cos(x);
        u1xy = @(t,x,y) -2.*K.*exp(t).*cos(2.*pi.*y).*sin(x);
        u1xx = @(t,x,y) -(K.*exp(t).*sin(2.*pi.*y).*cos(x))/pi;
        u1yy = @(t,x,y) -4.*K.*pi.*exp(t).*sin(2.*pi.*y).*cos(x);
        
        u2 = @(t,x,y) ((-2.*K+K/pi.^2.*(sin(pi.*y)).^2).*sin(x)).*exp(t);
        u2t = @(t,x,y) -exp(t).*sin(x).*(2.*K - (K.*sin(pi.*y).^2)/pi.^2);
        u2x = @(t,x,y) -exp(t).*cos(x).*(2.*K - (K.*sin(pi.*y).^2)/pi.^2);
        u2y = @(t,x,y) (2.*K.*exp(t).*cos(pi.*y).*sin(pi.*y).*sin(x))/pi;
        u2xy = @(t,x,y) (2.*K.*exp(t).*cos(pi.*y).*sin(pi.*y).*cos(x))/pi;
        u2xx = @(t,x,y) exp(t).*sin(x).*(2.*K - (K.*sin(pi.*y).^2)/pi.^2);
        u2yy = @(t,x,y) 2.*K.*exp(t).*cos(pi.*y).^2.*sin(x) - 2.*K.*exp(t).*sin(pi.*y).^2.*sin(x);

        p = @(t,x,y) 0.*x;
        px = @(t,x,y) 0.*x;
        py = @(t,x,y) 0.*x;
        
    case 2
        %---- Compute domain:
        % Stokes: [0,1]x[0,1]; Darcy: [0,1]x[1,2];
        %
        mu = 1; % fluid kinematic viscosity
        K = 1; % Darcy permeability matrix
        g = 1; % gravitational acceleration
        Kappa = 1; % this is the coefficient of (P(uh),P(vh))_interface
        pdecase = 'case2';
        
        phi = @(t,x,y) (y).*sin(pi.*x).*cos(t);
        phit = @(t,x,y) -y.*sin(pi.*x).*sin(t);
        phix = @(t,x,y) pi.*y.*cos(pi.*x).*cos(t);
        phiy = @(t,x,y) sin(pi.*x).*cos(t);
        phixy = @(t,x,y) pi.*cos(pi.*x).*cos(t);
        phixx = @(t,x,y) -pi.^2.*y.*sin(pi.*x).*cos(t);
        phiyy = @(t,x,y) 0.*x; 
        
        u1 = @(t,x,y) -sin(pi.*y).*cos(pi.*x).*cos(t);
        u1t = @(t,x,y) cos(pi.*x).*sin(pi.*y).*sin(t);
        u1x = @(t,x,y) pi.*sin(pi.*x).*sin(pi.*y).*cos(t);
        u1y = @(t,x,y) -pi.*cos(pi.*x).*cos(pi.*y).*cos(t);
        u1xy = @(t,x,y) pi.^2.*cos(pi.*y).*sin(pi.*x).*cos(t);
        u1xx = @(t,x,y) pi.^2.*cos(pi.*x).*sin(pi.*y).*cos(t);
        u1yy = @(t,x,y) pi.^2.*cos(pi.*x).*sin(pi.*y).*cos(t);

        u2 = @(t,x,y) sin(pi.*x).*cos(pi.*y).*cos(t);
        u2t = @(t,x,y) -cos(pi.*y).*sin(pi.*x).*sin(t);
        u2x = @(t,x,y) pi.*cos(pi.*x).*cos(pi.*y).*cos(t);
        u2y = @(t,x,y) -pi.*sin(pi.*x).*sin(pi.*y).*cos(t);
        u2xy = @(t,x,y) -pi.^2.*cos(pi.*x).*sin(pi.*y).*cos(t);
        u2xx = @(t,x,y) -pi.^2.*cos(pi.*y).*sin(pi.*x).*cos(t);
        u2yy = @(t,x,y) -pi.^2.*cos(pi.*y).*sin(pi.*x).*cos(t);

        p = @(t,x,y) sin(pi.*x).*cos(t);
        px = @(t,x,y) pi.*cos(pi.*x).*cos(t);
        py = @(t,x,y) 0.*x; 
        
        
    case 3
        K = @(x,y) k_function(x,y);

end

%---- initiation condition -------------------
phi0 = @(x,y) phi(0,x,y);
u10 = @(x,y) u1(0,x,y);
u20 = @(x,y) u2(0,x,y);
p0 = @(x,y) p(0,x,y);


%---- rhs ----------------------------------------
fp = @(t,x,y) phit(t,x,y)-K*(phixx(t,x,y) + phiyy(t,x,y)); % porous, the Darcy eqn's RHS f.
f1 = @(t,x,y) u1t(t,x,y)-mu*(2*u1xx(t,x,y)+u1yy(t,x,y)+u2xy(t,x,y))+px(t,x,y);
f2 = @(t,x,y) u2t(t,x,y)-mu*(u1xy(t,x,y)+u2xx(t,x,y)+2*u2yy(t,x,y))+py(t,x,y);

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
    'phi0', phi0, ...
    'u10', u10, ...
    'u20', u20, ...
    'p0', p0, ...
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


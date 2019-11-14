function pde = StokesDarcyDataUsingScalarStokes(u_case)
%
%   StokesDarcyDataUsingScalarStokes data for StokesDarcy using the Scalar Stokes equation
%
%
%	YcZhang 31/10/2017
%
%   Last modified 31/10/2017
%

switch u_case
    case 3
        % this case example is from: 
        % LiRui, 2017 A Weak Galerkin Finite Element Method for a Coupled Stokes-Darcy Problem
        % example 2.
        mu = 1;
        K = 1; % Darcy permeability matrix
        %g = 1;
        Kappa = 1; % this is the coefficient of (P(uh),P(vh))_interface
        pdecase = 'case3';
        
        phi = @(x,y) -pi/4*cos(pi*x/2).*y;
        phix = @(x,y) (y.*pi.^2.*sin((pi.*x)/2))/8;
        phiy = @(x,y) -(pi.*cos((pi.*x)/2))/4;
        phixy = @(x,y) (pi.^2.*sin((pi.*x)/2))/8;
        phixx = @(x,y) (y.*pi.^3.*cos((pi.*x)/2))/16; 
        phiyy = @(x,y) 0.*x;
        
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
        
        p = @(x,y) -pi/4*cos(pi*x/2).*(y-2*cos(pi*y/2).^2);
        px = @(x,y) (pi.^2.*sin((pi.*x)/2).*(y - 2.*cos((pi.*y)/2).^2))/8;
        py = @(x,y) -(pi.*cos((pi.*x)/2).*(2.*pi.*cos((pi.*y)/2).*sin((pi.*y)/2) + 1))/4;
        
        
    case 4
        %--------------------------------------
        % this case is WRONG!!!
        %--------------------------------------
        % this case example is from: 
        % DANAIL VASSILEV, 2009 Coupling Stokes–darcy flow with transport.
        % example 1.
        pdecase = 'case4';
        mu = 0.1;
        K = 1; % Darcy permeability matrix
        alpha = 0.5;
        G = sqrt(mu*K)/alpha;
        omega = 1.05;
        Kappa = mu*alpha/sqrt(K); % this is the coefficient of (P(uh),P(vh))_interface
        
        % Stokes Darcy data 
        phi = @(x,y) G/K*cos(x/G+omega).*exp(y/G);
        phix = @(x,y) -(exp(y/G).*sin(omega + x/G))/K;
        phiy = @(x,y) (exp(y/G).*cos(omega + x/G))/K;
        phixy = @(x,y) -(exp(y/G).*sin(omega + x/G))/(G.*K);
        phixx = @(x,y) -(exp(y/G).*cos(omega + x/G))/(G.*K); 
        phiyy = @(x,y) (exp(y/G).*cos(omega + x/G))/(G.*K);
        
        u1 = @(x,y) sin(x/G+omega).*exp(y/G);
        u1x = @(x,y) (exp(y/G).*cos(omega + x/G))/G;
        u1y = @(x,y) (exp(y/G).*sin(omega + x/G))/G;
        u1xy = @(x,y) (exp(y/G).*cos(omega + x/G))/G^2;
        u1xx = @(x,y) -(exp(y/G).*sin(omega + x/G))/G^2;
        u1yy = @(x,y) (exp(y/G).*sin(omega + x/G))/G^2;
        
        u2 = @(x,y) -cos(x/G+omega).*exp(y/G);
        u2x = @(x,y) (exp(y/G).*sin(omega + x/G))/G;
        u2y = @(x,y) -(exp(y/G).*cos(omega + x/G))/G;
        u2xy = @(x,y) (exp(y/G).*sin(omega + x/G))/G^2;
        u2xx = @(x,y) (exp(y/G).*cos(omega + x/G))/G^2;
        u2yy = @(x,y) -(exp(y/G).*cos(omega + x/G))/G^2;
        
        p = @(x,y) (G/K-mu/G)*cos(x/G+omega).*exp(1/(2*G))+y-0.5;
        px = @(x,y) -(exp(1/(2.*G)).*sin(omega + x/G).*(G/K - mu/G))/G;
        py = @(x,y) 1+0.*x;
    
    case 5
        %--------------------------------------
        % this case is WRONG!!!
        %--------------------------------------
        % this case example is from: 
        % DANAIL VASSILEV, 2009 Coupling Stokes–darcy flow with transport.
        % example 2. 
        pdecase = 'case5';
        Kappa = 1; % this is the coefficient of (P(uh),P(vh))_interface
        mu = 0.1;
        K = 1; % Darcy permeability matrix
        alpha = 0.5;
        G = sqrt(mu*K)/alpha;
        omega = 1.05;
        xi = (1-G)/(2+2*G);
        kappa = (-30*xi-17)/48;
        
        % Stokes Darcy data 
        phi = @(x,y) 1/K*(x.^2/2-2.*x).*(0.5-xi)+kappa/K*(-(y.^2+y)/2-1);
        phix = @(x,y) -((x - 2).*(xi - 1/2))/K;
        phiy = @(x,y) -(kappa.*(y + 1/2))/K;
        phixy = @(x,y) 0.*x;
        phixx = @(x,y) -(xi - 1/2)/K; 
        phiyy = @(x,y) -kappa/K;
        
        u1 = @(x,y) (2-x).*(1.5-y).*(y-xi);
        u1x = @(x,y) -(xi - y).*(y - 3/2);
        u1y = @(x,y) (x - 2).*(y - 3/2) - (xi - y).*(x - 2);
        u1xy = @(x,y) 2.*y - xi - 3/2;
        u1xx = @(x,y) 0.*x;
        u1yy = @(x,y) 2.*x - 4;
        
        u2 = @(x,y) -y.^3/3+y.^2/2.*(xi+1.5)-1.5*xi.*y-0.5;
        u2x = @(x,y) 0.*x;
        u2y = @(x,y) - y.^2 + (xi + 3/2).*y - (3.*xi)/2;
        u2xy = @(x,y) 0.*x;
        u2xx = @(x,y) 0.*x;
        u2yy = @(x,y) xi - 2.*y + 3/2;
        
        p = @(x,y) 1/K*(x.^2/2-2.*x).*(0.5-xi)-11*xi/(8*K)+mu*(0.5-xi)+y-0.5;
        px = @(x,y) -((x - 2).*(xi - 1/2))/K;
        py = @(x,y) 1+0.*x;
        
    case 6
        % this case example is from: 
        % LiRui, 2017 A Weak Galerkin Finite Element Method for a Coupled Stokes-Darcy Problem
        % example 1.
        mu = 1;
        K = 1; % Darcy permeability matrix
        %g = 1;
        Kappa = 1; % this is the coefficient of (P(uh),P(vh))_interface
        pdecase = 'case6';
        
        phi = @(x,y) (2-pi*sin(pi*x)).*(1-y-cos(pi*y));
        phix = @(x,y) pi.^2.*cos(pi.*x).*(y + cos(pi.*y) - 1);
        phiy = @(x,y) -(pi.*sin(pi.*x) - 2).*(pi.*sin(pi.*y) - 1);
        phixy = @(x,y) -pi.^2.*cos(pi.*x).*(pi.*sin(pi.*y) - 1);
        phixx = @(x,y) -pi.^3.*sin(pi.*x).*(y + cos(pi.*y) - 1); 
        phiyy = @(x,y) -pi.^2.*cos(pi.*y).*(pi.*sin(pi.*x) - 2);
        
        u1 = @(x,y) x.^2.*(y-1).^2+y;
        u1x = @(x,y) 2.*x.*(y - 1).^2;
        u1y = @(x,y) (2.*y - 2).*x.^2 + 1;
        u1xy = @(x,y) 2.*x.*(2.*y - 2);
        u1xx = @(x,y) 2.*(y - 1).^2;
        u1yy = @(x,y) 2.*x.^2;
        
        u2 = @(x,y) -2/3*x.*(y-1).^3+2-pi*sin(pi*x);
        u2x = @(x,y) - (2.*(y - 1).^3)/3 - pi.^2.*cos(pi.*x);
        u2y = @(x,y) -2.*x.*(y - 1).^2;
        u2xy = @(x,y) -2.*(y - 1).^2;
        u2xx = @(x,y) pi.^3.*sin(pi.*x);
        u2yy = @(x,y) -2.*x.*(2.*y - 2);
        
        p = @(x,y) (2-pi*sin(pi*x)).*sin(0.5*pi*y);
        px = @(x,y) -pi.^2.*cos(pi.*x).*sin((pi.*y)/2);
        py = @(x,y) -(pi.*cos((pi.*y)/2).*(pi.*sin(pi.*x) - 2))/2;
        
    case 7
        %----------------------------------
        %   This case, the true solution is wrong! 
        %   See, case 8.
        %----------------------------------
        % this case example is from: 
        % LiRui, 2017 A Weak Galerkin Finite Element Method for a Coupled Stokes-Darcy Problem
        % example 3.
        mu = 1;
        K = 1e-1; % Darcy permeability matrix
        %g = 1;
        Kappa = 1; % this is the coefficient of (P(uh),P(vh))_interface
        pdecase = 'case7';
        
        phi = @(x,y) 2*mu*x+(x.*(1-x).*(y-1)+1/3*y.^3-y.^2+y);
        phix = @(x,y) 2.*mu - (x - 1).*(y - 1) - x.*(y - 1);
        phiy = @(x,y) y.^2 - 2.*y - x.*(x - 1) + 1;
        phixy = @(x,y) 1 - 2.*x;
        phixx = @(x,y) 2 - 2.*y; 
        phiyy = @(x,y) 2.*y - 2;
        
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
        
        p = @(x,y) 2*mu*(x+y-1)+1/(3*K);
        px = @(x,y) 2*mu + 0.*x;
        py = @(x,y) 2*mu + 0.*x;
        
    case 8
        mu = 1;
        K = 1e-1; % Darcy permeability matrix
        %g = 1;
        Kappa = 1; % this is the coefficient of (P(uh),P(vh))_interface
        pdecase = 'case8';
        
        phi = @(x,y) 2.*mu.*x+(y-y.^2+y.^3/3-x.*(x-1).*(y-1))/K;
        phix = @(x,y) 2.*mu - ((x - 1).*(y - 1) + x.*(y - 1))/K;
        phiy = @(x,y) -(2.*y + x.*(x - 1) - y.^2 - 1)/K;
        phixy = @(x,y) -(2.*x - 1)/K;
        phixx = @(x,y) -(2.*y - 2)/K;
        phiyy = @(x,y) (2.*y - 2)/K;


        u1 = @(x,y) y.^2-2.*y+1;
        u1x = @(x,y) 0.*x; 
        u1y = @(x,y) 2.*y - 2;
        u1xy = @(x,y) 0.*x; 
        u1xx = @(x,y) 0.*x; 
        u1yy = @(x,y) 2;


        u2 = @(x,y) x.^2-x;
        u2x = @(x,y) 2.*x - 1;
        u2y = @(x,y) 0.*x; 
        u2xy = @(x,y) 0.*x; 
        u2xx = @(x,y) 2;
        u2yy = @(x,y) 0.*x; 


        p = @(x,y) 2.*mu.*(x+y-1)+1/(3.*K);
        px = @(x,y) 2.*mu;
        py = @(x,y) 2.*mu;
        
        
end

fp = @(x,y) -K*(phixx(x,y) + phiyy(x,y));
f1 = @(x,y) -mu*(u1xx(x,y)+u1yy(x,y))+px(x,y);
f2 = @(x,y) -mu*(u2xx(x,y)+u2yy(x,y))+py(x,y);
if strcmpi(pdecase,'case4') || strcmpi(pdecase,'case5') 
    fp = @(x,y) u1x(x,y) + u2y(x,y);
end 

gDp = phi;
gD1 = u1;
gD2 = u2;

funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1 + 0.*x;

pde = struct(...
    'mu', mu,  ...
    'K', K, ...
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


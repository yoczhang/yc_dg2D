function pde = TimeDecoupledSDTData2(coeff_case, u_case)
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   %--------------------------------------------
%       dgTimeDecoupledSDT2, 
%       for the transport eqn, have inflow and outflow boundaryEdges.
%   %---------------------------------------------
%
%
%
%	YcZhang 18/11/2017
%
%   Last modified 18/11/2017
%

DiffusivityCoeff = 1e-5;

switch coeff_case
    case 0
%         K = 1.0; % permeability matrix
%         K_11 = @(x,y) 0.*x;
%         K_11x = @(x,y) 0.*x; % \partial_x kappa
%         K_11y = @(x,y) 0.*x; % \partial_y kappa
%         K_22 = @(x,y) 0.*x;
%         K_22x = @(x,y) 0.*x; % \partial_x kappa
%         K_22y = @(x,y) 0.*x; % \partial_y kappa
        
        C_11 = @(x,y) DiffusivityCoeff+0.*x; % the coefficient parameter for transport_c function
        C_11x = @(x,y) 0.*x; % \partial_x 
        C_11y = @(x,y) 0.*x; % \partial_y 
        C_22 = @(x,y) DiffusivityCoeff+0.*x; % the coefficient parameter for transport_c function
        C_22x = @(x,y) 0.*x; % \partial_x 
        C_22y = @(x,y) 0.*x; % \partial_y 
end

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
        %c_darcy = K/mu; % this is the coefficient of 
        
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
        
        % Transport data
        c = @(t,x,y) t*(cos(pi*x) + cos(pi*y))/pi;
        ct = @(t,x,y) (cos(pi*x) + cos(pi*y))/pi;
        cx = @(t,x,y) -t*sin(pi*x);
        cy = @(t,x,y) -t*sin(pi*y);
        cxx = @(t,x,y) -t*pi*cos(pi*x);
        cyy = @(t,x,y) -t*pi*cos(pi*y);
        
    case 2
        K = @(x,y) K_function(x,y);
        
end

%---- initiation condition -------------------
phi0 = @(x,y) phi(0,x,y);
u10 = @(x,y) u1(0,x,y);
u20 = @(x,y) u2(0,x,y);
p0 = @(x,y) p(0,x,y);
c0 = @(t,x,y) c(0,x,y);

%---- rhs ----------------------------------------
fp = @(t,x,y) phit(t,x,y)-K*(phixx(t,x,y) + phiyy(t,x,y)); % porous, the Darcy eqn's RHS f.
f1 = @(t,x,y) u1t(t,x,y)-mu*(2*u1xx(t,x,y)+u1yy(t,x,y)+u2xy(t,x,y))+px(t,x,y);
f2 = @(t,x,y) u2t(t,x,y)-mu*(u1xy(t,x,y)+u2xx(t,x,y)+2*u2yy(t,x,y))+py(t,x,y);

fc_S = @(t,x,y) ct(t,x,y) + u1(t,x,y).*cx(t,x,y) + u2(t,x,y).*cy(t,x,y) ...
    - (C_11x(x,y).*cx(t,x,y)+C_11(x,y).*cxx(t,x,y)) - (C_22y(x,y).*cy(t,x,y)+C_22(x,y).*cyy(t,x,y));
fc_D = @(t,x,y) ct(t,x,y) + (-K)*phix(t,x,y).*cx(t,x,y) + (-K)*phiy(t,x,y).*cy(t,x,y) ...
    - (C_11x(x,y).*cx(t,x,y)+C_11(x,y).*cxx(t,x,y)) - (C_22y(x,y).*cy(t,x,y)+C_22(x,y).*cyy(t,x,y));
f_inflow = c;

%---- Dirichlet B.C. -------------------------
gDp = phi;
gD1 = u1;
gD2 = u2;
gDc = c;

%---- other auxiliary function ------------
funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1 + 0.*x;

pde = struct(...
    'DiffusivityCoeff', DiffusivityCoeff, ...
    'pdecase', pdecase, ...
    'Kappa', Kappa, ...
    'mu', mu,  ...
    'K', K, ...
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
    'c', c, ...
    'ct',ct, ...
    'cx', cx, ...
    'cxx', cxx, ...
    'cyy', cyy, ...
    'cy', cy, ...
    'phi0', phi0, ...
    'u10', u10, ...
    'u20', u20, ...
    'p0', p0, ...
    'c0', c0, ...
    'fp', fp, ...
    'f1',f1, ...
    'f2',f2, ...
    'fc_S',fc_S, ...
    'fc_D',fc_D, ...
    'f_inflow',f_inflow, ...
    'gDp', gDp, ...
    'gD1',gD1, ... % Dirichlet function
    'gD2',gD2, ... % Dirichlet function
    'gDc',gDc, ... % Dirichlet function
    'funcZero',funcZero, ...
    'funcOne',funcOne ...
    );

end % pde function

%%
%%>> -- Begin sub function  -------------------------------------------------------------
function r = K_function(x,y)

if 0<=x && x<=0.75 && 0.5<=y && y<=0.75
    r = 1e-4;
elseif 1.25<=x && x<=2 && 0.25<=y && y<=0.5
    r = 1e-6;
else
    r = 1;
    
end

end % function k_function



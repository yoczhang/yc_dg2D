function pde = SDTransportData(coeff_case, u_case)
%
%   %------------------------------------------------------
%       These examples is from the paper:
%       'coupling stokes-darcy flow with transport'--2009,SIAM, Dannil Vassilev.
%   %------------------------------------------------------
%
%   SDTransportData data for Stokes-Darcy-Transport equation
%
%
%	YcZhang 10/9/2017
%
%   Last modified 10/9/2017
%

switch coeff_case
    case 0
        k = 1.0;
        k_11 = @(x,y) 0.*x;
        k_11x = @(x,y) 0.*x; % \partial_x kappa
        k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 0.*x;
        k_22x = @(x,y) 0.*x; % \partial_x kappa
        k_22y = @(x,y) 0.*x; % \partial_y kappa
        
        C_11 = @(x,y) 1+0.*x; % the coefficient parameter for transport_c function
        C_11x = @(x,y) 0.*x; % \partial_x 
        C_11y = @(x,y) 0.*x; % \partial_y 
        C_22 = @(x,y) 1+0.*x; % the coefficient parameter for transport_c function
        C_22x = @(x,y) 0.*x; % \partial_x 
        C_22y = @(x,y) 0.*x; % \partial_y 
    case 1
        k = 1.0;
        k_11 = @(x,y) 1e-3+0.*x;
        k_11x = @(x,y) 0.*x; % \partial_x kappa
        k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1e-3+0.*x;
        k_22x = @(x,y) 0.*x; % \partial_x kappa
        k_22y = @(x,y) 0.*x; % \partial_y kappa

        C_11 = @(x,y) 1e-3+0.*x; % the coefficient parameter for transport_c function
        C_11x = @(x,y) 0.*x; % \partial_x 
        C_11y = @(x,y) 0.*x; % \partial_y 
        C_22 = @(x,y) 1e-3+0.*x; % the coefficient parameter for transport_c function
        C_22x = @(x,y) 0.*x; % \partial_x 
        C_22y = @(x,y) 0.*x; % \partial_y 
end

switch u_case
    case 1
        mu = 0.1;
        K = 1;
        alpha = 0.5;
        G = sqrt(mu*K)/alpha;
        omega = 1.05;
        
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
        
        % Transport data
        c = @(t,x,y) t*(cos(pi*x) + cos(pi*y))/pi;
        ct = @(t,x,y) (cos(pi*x) + cos(pi*y))/pi;
        cx = @(t,x,y) -t*sin(pi*x);
        cy = @(t,x,y) -t*sin(pi*y);
        cxx = @(t,x,y) -t*pi*cos(pi*x);
        cyy = @(t,x,y) -t*pi*cos(pi*y);
        
    case 2
        pdecase = 'case2';
        c_tau = 1; % this is the coefficient of (P(uh),P(vh))_interface
        mu = 0.1;
        K = 1;
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
        
        % Transport data
        c = @(t,x,y) t*(cos(pi*x) + cos(pi*y))/pi;
        ct = @(t,x,y) (cos(pi*x) + cos(pi*y))/pi;
        cx = @(t,x,y) -t*sin(pi*x);
        cy = @(t,x,y) -t*sin(pi*y);
        cxx = @(t,x,y) -t*pi*cos(pi*x);
        cyy = @(t,x,y) -t*pi*cos(pi*y);
        
        %--- this is just for test
        f_inflow = @(t,x,y) 0+0.*x;
        
    case 3 
        mu = 0.1;
        K = 1;
        alpha = 0.5;
        G = sqrt(mu*K)/alpha;
        omega = 6;
        xi = (1-G)/(2+2*G);
        kappa = (-30*xi-17)/48;
        
        % StokesDarcy data
        phi = @(x,y) -kappa/K*(y+0.5).^2/2-sin(omega*x).*y/K;
        phix = @(x,y) -(omega*y.*cos(omega*x))/K;
        phiy = @(x,y) - sin(omega*x)/K - (kappa*(2*y + 1))/(2*K);
        phixy = @(x,y) -(omega*cos(omega*x))/K;
        phixx = @(x,y) (omega.^2.*y.*sin(omega*x))/K; 
        phiyy = @(x,y) -kappa/K;
        
        u1 = @(x,y) (2-x).*(1.5-y).*(y-xi);
        u1x = @(x,y) -(xi - y).*(y - 3/2);
        u1y = @(x,y) (x - 2).*(y - 3/2) - (xi - y).*(x - 2);
        u1xy = @(x,y) 2.*y - xi - 3/2;
        u1xx = @(x,y) 0.*x;
        u1yy = @(x,y) 2.*x - 4;
        
        u2 = @(x,y) -y.^3/3+y.^2/2.*(xi+1.5)-1.5*xi.*y-0.5+sin(omega*x);
        u2x = @(x,y) omega*cos(omega*x);
        u2y = @(x,y) - y.^2 + (xi + 3/2).*y - (3.*xi)/2;
        u2xy = @(x,y) 0.*x;
        u2xx = @(x,y) -omega.^2.*sin(omega.*x);
        u2yy = @(x,y) xi - 2.*y + 3/2;
        
        p = @(x,y) -(sin(omega*x)+kappa)/(2*K)+mu*(0.5-xi)+cos(pi*y);
        px = @(x,y) -(omega*cos(omega*x))/(2*K);
        py = @(x,y) -pi*sin(pi*y);
        
        % Transport data
        c = @(t,x,y) t*(cos(pi*x) + cos(pi*y))/pi;
        ct = @(t,x,y) (cos(pi*x) + cos(pi*y))/pi;
        cx = @(t,x,y) -t*sin(pi*x);
        cy = @(t,x,y) -t*sin(pi*y);
        cxx = @(t,x,y) -t*pi*cos(pi*x);
        cyy = @(t,x,y) -t*pi*cos(pi*y);
        
    case 4
        k = @(x,y) k_function(x,y);
        
end

fp = @(x,y) -k*(phixx(x,y) + phiyy(x,y));
f1 = @(x,y) -mu*(2*u1xx(x,y)+u1yy(x,y)+u2xy(x,y))+px(x,y);
f2 = @(x,y) -mu*(u1xy(x,y)+u2xx(x,y)+2*u2yy(x,y))+py(x,y);
fc = @(t,x,y,flag) ct(t,x,y) + vector_c_1(x,y,flag).*cx(t,x,y) + vector_c_2(x,y,flag).*cy(t,x,y) ...
    - (C_11x(x,y).*cx(t,x,y)+C_11(x,y).*cxx(t,x,y)) - (C_22y(x,y).*cy(t,x,y)+C_22(x,y).*cyy(t,x,y));
fc_inflow = @(t,x,y) 1/2.^(t/10).*(cos((pi.*x)/5).*cos((pi.*y)/5) + 1) - (1/2.^(t/10).*pi.*cos((pi.*y)/5).*sin((pi.*x)/5))/5;

gDp = phi;
gD1 = u1;
gD2 = u2;
gDc = c;

c0 = @(x,y) c(0,x,y);

funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1 + 0.*x;

pde = struct(...
    'pdecase', pdecase, ...
    'c_tau', c_tau, ...
    'mu', mu,  ...
    'k', k, ...
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
    'c0', c0, ...
    'f_inflow', f_inflow, ...
    'fp', fp, ...
    'f1',f1, ...
    'f2',f2, ...
    'fc',fc, ...
    'fc_inflow',fc_inflow, ...
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
function r = k_function(x,y)

if 0<=x && x<=0.75 && 0.5<=y && y<=0.75
    r = 1e-4;
elseif 1.25<=x && x<=2 && 0.25<=y && y<=0.5
    r = 1e-6;
else
    r = 1;
    
end

end % function k_function



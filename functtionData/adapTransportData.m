function pde = adapTransportData(coeff_case, u_case)
%
%   transportData data for transport equation
%
%
%	YcZhang 23/9/2017
%
%   Last modified 23/9/2017
%

diffusivityCoeff = 1e-5;

switch coeff_case
    case 0
        k_11 = @(x,y) diffusivityCoeff+0.*x; % the coefficient parameter for Poisson function
        k_11x = @(x,y) 0.*x; % \partial_x kappa
        k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) diffusivityCoeff+0.*x;
        k_22x = @(x,y) 0.*x; % \partial_x kappa
        k_22y = @(x,y) 0.*x; % \partial_y kappa
        vector_u_1 = @(x,y,flag) external_velocity(x,y,flag,1);
        vector_u_2 = @(x,y,flag) external_velocity(x,y,flag,2);
end


switch u_case
    case 1
        
        c = @(t,x,y) t*(cos(pi*x) + cos(pi*y))/pi;
        ct = @(t,x,y) (cos(pi*x) + cos(pi*y))/pi;
        cx = @(t,x,y) -t*sin(pi*x);
        cy = @(t,x,y) -t*sin(pi*y);
        cxx = @(t,x,y) -t*pi*cos(pi*x);
        cyy = @(t,x,y) -t*pi*cos(pi*y);
        
    case 2
        c = @(t,x,y) (1+cos(pi*x/5).*cos(pi*y/5))*2^(-t/10);
        ct = @(t,x,y) -(1/2.^(t/10).*log(2).*(cos((pi.*x)/5).*cos((pi.*y)/5) + 1))/10;
        cx = @(t,x,y) -(1/2.^(t/10).*pi.*cos((pi.*y)/5).*sin((pi.*x)/5))/5;
        cy = @(t,x,y) -(1/2.^(t/10).*pi.*cos((pi.*x)/5).*sin((pi.*y)/5))/5;
        cxx = @(t,x,y) -(1/2.^(t/10).*pi.^2.*cos((pi.*x)/5).*cos((pi.*y)/5))/25;
        cyy = @(t,x,y) -(1/2.^(t/10).*pi.^2.*cos((pi.*x)/5).*cos((pi.*y)/5))/25;
        
    case 3
        c = @(t,x,y) t^4*(cos(pi*x) + cos(pi*y))/pi;
        ct = @(t,x,y) 4*t^3*(cos(pi*x) + cos(pi*y))/pi;
        cx = @(t,x,y) -t^4*sin(pi*x);
        cy = @(t,x,y) -t^4*sin(pi*y);
        cxx = @(t,x,y) -t^4*pi*cos(pi*x);
        cyy = @(t,x,y) -t^4*pi*cos(pi*y);
        
    case 4
        c = @(t,x,y) exp(-t)*(cos(pi*x) + cos(pi*y))/pi;
        ct = @(t,x,y) -exp(-t)*(cos(pi*x) + cos(pi*y))/pi;
        cx = @(t,x,y) -exp(-t)*sin(pi*x);
        cy = @(t,x,y) -exp(-t)*sin(pi*y);
        cxx = @(t,x,y) -exp(-t)*pi*cos(pi*x);
        cyy = @(t,x,y) -exp(-t)*pi*cos(pi*y);
end 


f_rhs = @(t,x,y,flag) ct(t,x,y) + vector_u_1(x,y,flag).*cx(t,x,y) + vector_u_2(x,y,flag).*cy(t,x,y) ...
    - (k_11x(x,y).*cx(t,x,y)+k_11(x,y).*cxx(t,x,y)) - (k_22y(x,y).*cy(t,x,y)+k_22(x,y).*cyy(t,x,y));
% f_rhs = @(t,x,y,flag) 0.*x;
f_rhs_test = @(t,x,y,flag) (1/2.^(t/10).*(8.*pi.^2.*cos((pi.*x)/5).*cos((pi.*y)/5) - log(2) + 8.*pi.*cos((pi.*x)/5).*sin((pi.*y)/5) + 20.*pi.*cos((pi.*y)/5).*sin((pi.*x)/5) - log(2).*cos((pi.*x)/5).*cos((pi.*y)/5)))/100;
f_inflow = @(t,x,y) 1/2.^(t/10).*(cos((pi.*x)/5).*cos((pi.*y)/5) + 1) - (1/2.^(t/10).*pi.*cos((pi.*y)/5).*sin((pi.*x)/5))/5;

gD = c;
u0 =@(t,x,y,flag) c(0,x,y);

funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1 + 0.*x;


pde = struct(...
    'diffusivityCoeff', diffusivityCoeff, ...
    'k_11',k_11, 'k_22', k_22, ...
    'vector_u_1',vector_u_1,  ...
    'vector_u_2',vector_u_2,  ...
    'u',c, ...
    'ux',cx, ...
    'uy',cy, ...
    'f_rhs',f_rhs, ...
    'f_rhs_test',f_rhs_test, ...
    'f_inflow', f_inflow, ...
    'gD',gD, ... % Dirichlet function
    'u0',u0, ...
    'funcZero',funcZero, ...
    'funcOne',funcOne ...
    );


end % function



%%
%%>> -- Begin sub function  -------------------------------------------------------------
function r = external_velocity(x,y,flag,component)

%% case 1

%--- 
% if flag == 1 && component == 1
% %     r = 1+0.*x;
% %     r = -1+0.*x;
%     r = 1e-0+0.*x;
% elseif flag == 1 && component == 2
% %     r = 0.5+0.*y;
% %     r = -0.4+0.*y;
%     r = 0.*x;
% elseif flag == 2 && component == 1
%     r = 0+0.*x;
% elseif flag == 2 && component == 2
%     r = 0+0.*y;
% end

%% case 2

%--- Stokes-Darcy velocity setting
mu = 1e-2;
K = 1e-2;
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

if flag == 1 && component == 1 % Stokes, veloctiy1
    r = u1(x,y);
elseif flag == 1 && component == 2
    r = u2(x,y);
elseif flag == -1 && component == 1
    r = -K*phix(x,y);
elseif flag == -1 && component == 2
    r = -K*phiy(x,y);
end 

end % function k_function


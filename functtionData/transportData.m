function pde = transportData(coeff_case, u_case)
%
%   transportData data for transport equation
%
%
%	YcZhang 3/9/2017
%
%   Last modified 3/9/2017
%

switch coeff_case
    case 0
        k_11 = @(x,y) 1+0.*x; % the coefficient parameter for Poisson function
        k_11x = @(x,y) 0.*x; % \partial_x kappa
        k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1+0.*x;
        k_22x = @(x,y) 0.*x; % \partial_x kappa
        k_22y = @(x,y) 0.*x; % \partial_y kappa
        vector_u_1 = @(x,y,flag) external_velocity(x,y,flag,1);
        vector_u_2 = @(x,y,flag) external_velocity(x,y,flag,2);
    case 1
        k_11 = @(x,y) 1e-3+0.*x;
        k_11x = @(x,y) 0.*x; % \partial_x kappa
        k_11y = @(x,y) 0.*x; % \partial_y kappa
        k_22 = @(x,y) 1e-3+0.*x;
        k_22x = @(x,y) 0.*x; % \partial_x kappa
        k_22y = @(x,y) 0.*x; % \partial_y kappa
        vector_u_1 = @(x,y,flag) external_velocity(x,y,flag,1);
        vector_u_2 = @(x,y,flag) external_velocity(x,y,flag,2);
    case 2
        k_11 = @(x,y) 2+sin(x).*sin(y);
        k_11x = @(x,y) cos(x).*sin(y); % \partial_x kappa
        k_11y = @(x,y) sin(x).*cos(y); % \partial_y kappa
        k_22 = @(x,y) 2+sin(x).*sin(y);
        k_22x = @(x,y) cos(x).*sin(y); % \partial_x kappa
        k_22y = @(x,y) sin(x).*cos(y); % \partial_y kappa
end


switch u_case
    case 1
        u = @(t,x,y) t*(cos(pi*x) + cos(pi*y))/pi;
        ut = @(t,x,y) (cos(pi*x) + cos(pi*y))/pi;
        ux = @(t,x,y) -t*sin(pi*x);
        uy = @(t,x,y) -t*sin(pi*y);
        uxx = @(t,x,y) -t*pi*cos(pi*x);
        uyy = @(t,x,y) -t*pi*cos(pi*y);
        
    case 2
        u = @(t,x,y) (1+cos(pi*x/5).*cos(pi*y/5))*2^(-t/10);
        ut = @(t,x,y) -(1/2.^(t/10).*log(2).*(cos((pi.*x)/5).*cos((pi.*y)/5) + 1))/10;
        ux = @(t,x,y) -(1/2.^(t/10).*pi.*cos((pi.*y)/5).*sin((pi.*x)/5))/5;
        uy = @(t,x,y) -(1/2.^(t/10).*pi.*cos((pi.*x)/5).*sin((pi.*y)/5))/5;
        uxx = @(t,x,y) -(1/2.^(t/10).*pi.^2.*cos((pi.*x)/5).*cos((pi.*y)/5))/25;
        uyy = @(t,x,y) -(1/2.^(t/10).*pi.^2.*cos((pi.*x)/5).*cos((pi.*y)/5))/25;
        
    case 3
        u = @(t,x,y) t^4*(cos(pi*x) + cos(pi*y))/pi;
        ut = @(t,x,y) 4*t^3*(cos(pi*x) + cos(pi*y))/pi;
        ux = @(t,x,y) -t^4*sin(pi*x);
        uy = @(t,x,y) -t^4*sin(pi*y);
        uxx = @(t,x,y) -t^4*pi*cos(pi*x);
        uyy = @(t,x,y) -t^4*pi*cos(pi*y);
        
    case 4
        u = @(t,x,y) exp(-t)*(cos(pi*x) + cos(pi*y))/pi;
        ut = @(t,x,y) -exp(-t)*(cos(pi*x) + cos(pi*y))/pi;
        ux = @(t,x,y) -exp(-t)*sin(pi*x);
        uy = @(t,x,y) -exp(-t)*sin(pi*y);
        uxx = @(t,x,y) -exp(-t)*pi*cos(pi*x);
        uyy = @(t,x,y) -exp(-t)*pi*cos(pi*y);
end 


f_rhs = @(t,x,y,flag) ut(t,x,y) + vector_u_1(x,y,flag).*ux(t,x,y) + vector_u_2(x,y,flag).*uy(t,x,y) ...
    - (k_11x(x,y).*ux(t,x,y)+k_11(x,y).*uxx(t,x,y)) - (k_22y(x,y).*uy(t,x,y)+k_22(x,y).*uyy(t,x,y));
f_rhs_test = @(t,x,y,flag) (1/2.^(t/10).*(8.*pi.^2.*cos((pi.*x)/5).*cos((pi.*y)/5) - log(2) + 8.*pi.*cos((pi.*x)/5).*sin((pi.*y)/5) + 20.*pi.*cos((pi.*y)/5).*sin((pi.*x)/5) - log(2).*cos((pi.*x)/5).*cos((pi.*y)/5)))/100;
f_inflow = @(t,x,y) 1/2.^(t/10).*(cos((pi.*x)/5).*cos((pi.*y)/5) + 1) - (1/2.^(t/10).*pi.*cos((pi.*y)/5).*sin((pi.*x)/5))/5;

gD = u;
u0 =@(x,y) u(0,x,y);

funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1 + 0.*x;


pde = struct(...
    'k_11',k_11, 'k_22', k_22, ...
    'vector_u_1',vector_u_1,  ...
    'vector_u_2',vector_u_2,  ...
    'u',u, ...
    'ux',ux, ...
    'uy',uy, ...
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

if flag == 1 && component == 1
%     r = 1+0.*x;
    r = -1+0.*x;
elseif flag == 1 && component == 2
%     r = 0.5+0.*y;
    r = -0.4+0.*y;
elseif flag == 2 && component == 1
    r = 0+0.*x;
elseif flag == 2 && component == 2
    r = 0+0.*y;
end

end % function k_function


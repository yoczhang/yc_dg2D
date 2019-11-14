function r=function_f(x,y,t,flag_ele)

if flag_ele==1
    r=(1/2^(t/10)*(8*pi^2*cos((pi*x)/5)*cos((pi*y)/5) - log(2) + 8*pi*cos((pi*x)/5)*sin((pi*y)/5) + 20*pi*cos((pi*y)/5)*sin((pi*x)/5) - log(2)*cos((pi*x)/5)*cos((pi*y)/5)))/100;
else
    r=(2/2^(t/10)*pi^2*cos((pi*x)/5)*cos((pi*y)/5))/25 - (1/2^(t/10)*log(2)*(cos((pi*x)/5)*cos((pi*y)/5) + 1))/100;
end
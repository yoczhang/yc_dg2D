function r=function_inflow(x,y,t)

r=1/2^(t/10)*(cos((pi*x)/5)*cos((pi*y)/5) + 1) - (1/2^(t/10)*pi*cos((pi*y)/5)*sin((pi*x)/5))/5;
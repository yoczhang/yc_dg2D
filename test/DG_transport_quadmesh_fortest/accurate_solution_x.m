function r=accurate_solution_x(x,y,t)

r=-(1/2^(t/10)*pi*cos((pi*y)/5)*sin((pi*x)/5))/5;
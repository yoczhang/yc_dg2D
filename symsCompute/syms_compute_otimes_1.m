
clearvars

syms A11 A12 A21 A22 a1 a2 b1 b2

A=[A11 A12;A21 A22];
a = [a1;a2];
b = [b1;b2];

c1=dot(A*a,b)
c2=